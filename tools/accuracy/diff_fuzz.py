"""Find the inputs where the working tree (or `--head`) is less accurate than a base revision.

Usage (from the repository root; or `just diff-fuzz`):

    uv run --project tools tools/accuracy/diff_fuzz.py [--base REV] [--head REV]
        [--only NAME ...] [--cases N] [--seed S] [--max-refs N] [--jobs N]

(e.g. `--cases 10 --max-refs 10` for a quick run of all specs).

Steps:

1. Build the dump program of `accuracy.py` for both revisions; a revision other than the working
   tree is exported with `git archive` and built in a temporary directory, and its program is
   cached in `.cache/diff-fuzz/`.
2. Evaluate both on the same `--cases` seeded inputs per spec of `specs.py` (or of `--only`, by
   spec or function name): samples of its domains, with the parameters scaled log-uniformly or
   rounded to (half-)integers, single arguments moved to log-uniform, (half-)integer, boundary
   or special values, a parameter or real argument set to (nearly) a signed sum of others plus
   an integer, and the imaginary parts of complex arguments set to 0 or tiny values of either
   sign. Skip the specs whose function a revision does not have.
3. Flag the inputs where the results differ by more than `RTOL` relative to the larger, in a part
   of a complex result, or where one is NaN, infinite, zero, a panic or a timeout and the other is
   not, or where the signs of zeros differ.
4. Count a panic or timeout of one revision as a regression (or an improvement). For a share of
   the other flagged inputs (at most `--max-refs` in all; those that differ in kind first),
   compute references in `REF_TIMEOUT` seconds, with Arb where `specs.ARB` or `ARB` has the
   function (NaN if it is not real for a real function), else with the reference of the spec.
   Classify each by the largest relative error of the parts that differ, from the reference
   rounded to binary64: regressed (more than `WORSE` times the base error, and more than
   `FLOOR`), improved (vice versa), both wrong, or both right.
5. For `NEAR_ROUNDS` rounds, repeat 3 and 4 for `--cases` mutants of the regressions and the
   flagged inputs so far (or, if none are flagged yet, of new inputs as in 2), with the
   parameters scaled (keeping that sum), two arguments shifted together, or one moved slightly
   or as in 2.
6. Print the counts, write the regressions to `.cache/diff-fuzz.tsv`, and exit with 1 if there
   are any.
"""

import argparse
import concurrent.futures
import contextlib
import io
import math
import multiprocessing
import os
import random
import shutil
import struct
import subprocess
import sys
import tarfile
import tempfile
import time
from collections import Counter
from collections.abc import Callable
from enum import StrEnum
from multiprocessing.connection import Connection, wait
from pathlib import Path
from typing import Final, NamedTuple

import flint
import mpmath as mp
import specs

import accuracy

RTOL: Final = 1e-12
WORSE: Final = 2.0
FLOOR: Final = 4e-16
# a part of a reference below 2^-TINY_BITS of its modulus counts as 0, so that its error is
# relative to the modulus
TINY_BITS: Final = 200
# the bits of the Arb references, relative to their modulus, so that the other parts have 64
# (more than the 53 of binary64)
ARB_PREC: Final = 64 + TINY_BITS
# the magnitude of a signed zero imaginary argument for Arb, for the side of a branch cut: far
# below what its working precision resolves
ARB_ZERO: Final = flint.arb(2) ** -(2 * specs.ARB_MAX_PREC)
# seconds per reference, so that a few slow ones cannot stall a run
REF_TIMEOUT: Final = 10
CACHE: Final = accuracy.CACHE / "diff-fuzz"
TSV: Final = accuracy.CACHE / "diff-fuzz.tsv"
# lines per run of a dump program; a failed run is bisected down to the failing inputs
CHUNK: Final = 1000
# seconds per run of a dump program, and per line in it
RUN_TIMEOUT: Final = 10
LINE_TIMEOUT: Final = 0.02
# seconds per spec and revision that timeouts of runs may cost, after which the lines of a timed
# out run are not bisected further (a library call can hang)
TIMEOUT_BUDGET: Final = 120
# references per worker process, since mpmath caches grow
WORKER_TASKS: Final = 500
NEAR_ROUNDS: Final = 3


def _hyp1f1_c(a: flint.arb, b: flint.arb, re: flint.arb, im: flint.arb) -> flint.acb:
    return flint.acb(re, im).hypgeom_1f1(a, b)


def _hyp2f1_c(
    a: flint.arb, b: flint.arb, c: flint.arb, re: flint.arb, im: flint.arb
) -> flint.acb:
    return flint.acb(re, im).hypgeom_2f1(a, b, c)


# with the complex functions that `specs.ARB` lacks, since their specs use mpmath
ARB: Final = specs.ARB | {"hyp1f1#c": _hyp1f1_c, "hyp2f1#c": _hyp2f1_c}
IMAG: Final = {"im", "xi", "yi"}
SPECIAL: Final = (0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 5e-324, 1e-300, 1e300)
SPECIAL_NONFINITE: Final = (math.inf, -math.inf, math.nan)

type Args = accuracy.Args
# a result, or "unsupported", "panic", "timeout" or "skipped" (not run, after too many timeouts)
type Value = complex | str
type RefParts = list[str] | str  # the parts of a reference, or why there is none


class Kind(StrEnum):
    """The verdict on a flagged input."""

    REGRESSED = "regressed"
    IMPROVED = "improved"
    BOTH_WRONG = "both wrong"
    BOTH_RIGHT = "both right"
    NO_REF = "no ref"


COLUMNS: Final = ["function", "cases", "flagged", *Kind, "unchecked"]


class Budget:
    """The seconds that timeouts may still cost for one revision of a spec."""

    def __init__(self) -> None:
        self.left: float = TIMEOUT_BUDGET


def program(rev: str | None) -> str:
    """The dump program of the working tree (`rev` None), or of a commit (cached)."""
    if rev is None:
        return accuracy.build_rust()
    git = ["git", "-C", str(accuracy.ROOT)]
    sha = subprocess.check_output(
        [*git, "rev-parse", f"{rev}^{{commit}}"], text=True
    ).strip()
    exe = CACHE / f"xsf-accuracy-{sha}"
    if not exe.exists():
        accuracy.log(f"building {rev} ({sha[:12]})")
        # (outside the working tree, so that cargo takes no config from it, and published by a
        # rename, for concurrent runs)
        with tempfile.TemporaryDirectory() as tmp:
            tar = subprocess.check_output([*git, "archive", sha])
            with tarfile.open(fileobj=io.BytesIO(tar)) as tf:
                tf.extractall(tmp, filter="data")
            built = accuracy.build_rust(Path(tmp), Path(tmp) / "target")
            fd, part = tempfile.mkstemp(dir=CACHE)
            os.close(fd)
            shutil.copy2(built, part)
            os.replace(part, exe)
    return str(exe)


def evaluate(exe: str, lines: list[tuple[str, Args]], budget: Budget) -> list[Value]:
    """The results of the input lines, or why the program fails on them."""
    if len(lines) > CHUNK:
        return [
            v
            for i in range(0, len(lines), CHUNK)
            for v in evaluate(exe, lines[i : i + CHUNK], budget)
        ]
    if budget.left <= 0:
        return ["skipped"] * len(lines)
    timeout = RUN_TIMEOUT + LINE_TIMEOUT * len(lines)
    try:
        return list[Value](accuracy.run([exe, "xsf"], lines, timeout))
    except subprocess.CalledProcessError as e:
        fail = "unsupported" if "unknown function" in e.stderr else "panic"
    except subprocess.TimeoutExpired:
        fail = "timeout"
        budget.left -= timeout
    if len(lines) == 1:
        return [fail]
    mid = len(lines) // 2
    return evaluate(exe, lines[:mid], budget) + evaluate(exe, lines[mid:], budget)


class Case(NamedTuple):
    fuzzer: Fuzzer
    args: Args
    # a parameter i set to sum(sign * args[j]) + offset, as (i, [(sign, j), ...], offset)
    rel: tuple[int, list[tuple[int, int]], float] | None = None

    def related(self, a: list[float]) -> Case:
        """The case of the arguments `a`, with the parameter of `rel` set again."""
        if self.rel is not None:
            i, terms, offset = self.rel
            a[i] = sum(s * a[j] for s, j in terms) + offset
        return self._replace(args=tuple(a))


class Fuzzer:
    """Mutated samples of a domain of a spec."""

    def __init__(self, spec: specs.Spec, dom: specs.Domain, rng: random.Random) -> None:
        self.dom, self.name = dom, dom.call or spec.func
        # the arguments that are mutated (of an array, its first element), but not the integers
        # (by their distribution, or for a sampler, by its samples)
        n = len(spec.args)
        samples = [dom.sample(rng) for _ in range(20)]
        dists = dom.dists or [None] * n
        self.free = [
            i
            for i in range(n)
            if not isinstance(dists[i], specs.Int | specs.ILU | specs.Each)
            and not all(s[i].is_integer() for s in samples)
        ]
        self.imag = [i for i in self.free if spec.args[i] in IMAG]
        # the parameters: all but the complex arguments, or else the last argument
        last = -1 if spec.complex else len(spec.args) - 1
        self.params = [
            i
            for i in self.free
            if i != last and i not in self.imag and i + 1 not in self.imag
        ]
        # what a relation links: the parameters, and a real argument
        self.linked = self.params + [i for i in self.free if i == last]
        self.bounds = {i: [0.0, 1.0, -1.0, *_bounds(dists[i])] for i in self.free}

    def sample(self, rng: random.Random) -> Case:
        a = list(self.dom.sample(rng))
        if rng.random() < 0.3:
            self.scale(rng, a)
        if (r := rng.random()) < 0.3:
            k = 1 if r < 0.2 else 2
            for i in self.params:
                a[i] = round(k * a[i], 0) / k
        # (the arguments more often, e.g. to a singular point)
        for i in self.free:
            if rng.random() < (0.25 if i in self.params else 0.5):
                a[i] = self.mutate(rng, i, a[i])
        self.cut(rng, a, 0.5)
        if len(self.linked) < 2 or rng.random() >= 0.4:
            return Case(self, tuple(a))
        # e.g. a balanced c = a + b + m for hyp2f1, a2 = a1 + m, or x = v + m for Bessel
        k = 3 if len(self.linked) > 2 and rng.random() < 2 / 3 else 2
        i, *js = rng.sample(self.linked, k)
        terms = [(rng.choice((-1, 1)), j) for j in js]
        return Case(self, (), (i, terms, rng.randint(-2, 2) + _small(rng))).related(a)

    def near(self, rng: random.Random, case: Case) -> Case:
        """A mutant of `case` with one or two changes."""
        a, rel = list(case.args), case.rel
        for _ in range(rng.randint(1, 2) if self.free else 0):
            i, j = rng.choice(self.free), rng.choice(self.free)
            match rng.randrange(4):
                case 0:
                    self.scale(rng, a)
                case 1:
                    t = rng.choice((-1, 1)) * 10 ** rng.uniform(-3, 2)
                    a[i] += t
                    a[j] += t if j != i else 0
                case 2:
                    a[i] *= 1 + _small(rng)
                case _:
                    a[i] = self.mutate(rng, i, a[i])
            if rel is not None and rel[0] in (i, j):
                rel = None
        self.cut(rng, a, 0.25)
        return Case(self, (), rel).related(a)

    def scale(self, rng: random.Random, a: list[float]) -> None:
        s = 10 ** rng.uniform(-2, 2)
        for i in self.params:
            a[i] *= s

    def cut(self, rng: random.Random, a: list[float], p: float) -> None:
        """Set the imaginary parts to 0 or tiny values of either sign, each with probability p."""
        for i in self.imag:
            if rng.random() < p:
                tiny = rng.choice((0.0, 10 ** rng.uniform(-320, -1)))
                re = a[i - 1] if math.isfinite(a[i - 1]) else 1.0
                a[i] = rng.choice((-1, 1)) * tiny * max(1.0, abs(re))

    def mutate(self, rng: random.Random, i: int, x: float) -> float:
        if not math.isfinite(x):
            return x
        # (an argument more often to a boundary, e.g. a singular point)
        weights = (2, 2, 2, 2, 1, 0.2) if i in self.params else (1, 1, 1, 4, 1, 0.2)
        match rng.choices(range(6), weights)[0]:
            case 0:
                top = math.log10(max(*map(abs, self.bounds[i]), 1.0)) + 2
                mag = 10 ** min(rng.uniform(-8, top), sys.float_info.max_10_exp)
                return rng.choice((-1, 1)) * mag
            case 1:
                return _nudge(rng, round(x, 0))
            case 2:
                return _nudge(rng, round(2 * x, 0) / 2)
            case 3:
                return _nudge(rng, rng.choice(self.bounds[i]))
            case 4:
                return rng.choice(SPECIAL)
            case _:
                return rng.choice(SPECIAL_NONFINITE)


def _bounds(d: specs.Dist | None) -> list[float]:
    match d:
        case specs.U(lo, hi):
            return [lo, hi]
        case specs.LU(lo, hi, sign):
            return [s * v for v in (lo, hi) for s in ((1, -1) if not sign else (sign,))]
        case specs.Offset(c, _, _):
            return [c]
    return []


def _small(rng: random.Random) -> float:
    """0, or a power of 2 in [2^-52, 1/2] of a random sign."""
    if rng.random() < 1 / 3:
        return 0.0
    return rng.choice((-1, 1)) * 2.0 ** -rng.randint(1, 52)


def _nudge(rng: random.Random, t: float) -> float:
    """`t`, possibly moved by `_small` relative to it (or to 1, if smaller)."""
    return t + _small(rng) * max(1.0, abs(t)) if math.isfinite(t) else t


def _parts(y: complex) -> tuple[float, ...]:
    return (y.real, y.imag) if isinstance(y, complex) else (y,)


def diffs(x: Value, y: Value) -> list[float]:
    """How much the parts of two results differ: 0 if by at most `RTOL` relative to the larger,
    `RTOL` for zeros of different signs, else their relative difference, or infinity if they
    differ in kind (NaN, infinite, zero, or a failure)."""
    if "skipped" in (x, y):
        return [0.0]
    if isinstance(x, str) or isinstance(y, str):
        return [0.0 if x == y else math.inf]
    return list(map(_diff, _parts(x), _parts(y)))


def _diff(u: float, v: float) -> float:
    if math.isnan(u) or math.isnan(v):
        return 0.0 if math.isnan(u) and math.isnan(v) else math.inf
    if u == v:
        return 0.0 if math.copysign(1, u) == math.copysign(1, v) else RTOL
    if not (math.isfinite(u) and math.isfinite(v)) or u == 0 or v == 0:
        return math.inf
    d = abs(u - v) / max(abs(u), abs(v))
    return d if d > RTOL else 0.0


def ref_parts(task: tuple[str, Args]) -> RefParts:
    """The reference (its parts), or why there is none."""
    name, args = task
    spec = specs.BY_NAME[name]
    if name in ARB:
        xs = [
            math.copysign(1, x) * ARB_ZERO if x == 0 and arg in IMAG else flint.arb(x)
            for x, arg in zip(args, spec.args, strict=True)
        ]
        y = specs.arb_ball(lambda: flint.acb(ARB[name](*xs)), ARB_PREC)
        if y is not None:
            if spec.complex:
                parts = [y.real, y.imag]
            elif y.imag.contains(0):
                parts = [y.real]
            else:  # (a complex continuation)
                return ["nan"]
            return [p.mid().str(40, radius=False) for p in parts]
    match accuracy.ref_task(task):
        case [list(v), _]:
            return v
        case [str(v), _]:
            return [v]
        case e:
            return str(e)


def _worker(conn: Connection) -> None:
    # (until the parent closes its end)
    with contextlib.suppress(EOFError):
        while True:
            conn.send(ref_parts(conn.recv()))


class Refs:
    """Worker processes for `ref_parts`, each killed and replaced when a reference takes longer
    than `REF_TIMEOUT` (or the timeout of its spec), or after `WORKER_TASKS` references."""

    def __init__(self, n: int) -> None:
        self.n = n
        self.count = Counter[multiprocessing.Process]()
        self.idle: list[tuple[multiprocessing.Process, Connection]] = []

    def map(self, tasks: list[tuple[str, Args]]) -> list[RefParts]:
        out: list[RefParts] = ["timeout"] * len(tasks)
        todo = list(enumerate(tasks))[::-1]
        busy: dict[Connection, tuple[multiprocessing.Process, int, float]] = {}
        while todo or busy:
            while todo and len(busy) < self.n:
                proc, conn = self.idle.pop() if self.idle else _start()
                i, task = todo.pop()
                conn.send(task)
                limit = min(specs.BY_NAME[task[0]].timeout, REF_TIMEOUT)
                busy[conn] = proc, i, time.monotonic() + limit
            first = min(t for _, _, t in busy.values())
            ready = wait(list(busy), max(first - time.monotonic(), 0))
            for conn in [c for c in busy if c in ready]:
                proc, i, _ = busy.pop(conn)
                try:
                    out[i] = conn.recv()
                    self.count[proc] += 1
                    if self.count[proc] < WORKER_TASKS:
                        self.idle.append((proc, conn))
                    else:
                        proc.kill()
                except EOFError:
                    out[i] = "worker died"
            now = time.monotonic()
            for conn, (proc, _, t) in list(busy.items()):
                if t <= now:
                    proc.kill()
                    del busy[conn]
        return out


def _start() -> tuple[multiprocessing.Process, Connection]:
    conn, child = multiprocessing.Pipe()
    proc = multiprocessing.Process(target=_worker, args=(child,), daemon=True)
    proc.start()
    return proc, conn


def error(y: Value, ref: list[str], parts: list[bool]) -> float:
    """The largest relative error of the given parts from those of `ref` rounded to binary64,
    where the parts below 2^-TINY_BITS of the modulus (of the finite rounded parts) count as 0
    and their errors are relative to that."""
    if isinstance(y, str):
        return math.inf
    mod = mp.sqrt(sum(mp.mpf(r) ** 2 for r in map(float, ref) if math.isfinite(r)))
    return max(
        (_error(u, s, mod) for u, s, p in zip(_parts(y), ref, parts, strict=True) if p),
        default=0.0,
    )


def _error(u: float, s: str, mod: mp.mpf) -> float:
    r, v = float(s), mp.mpf(s)  # (r correctly rounded)
    if not (math.isfinite(u) and math.isfinite(r)):
        return 0.0 if u == r or (math.isnan(u) and math.isnan(r)) else math.inf
    if abs(v) <= mp.ldexp(mod, -TINY_BITS):
        return float(abs(u) / mod) if mod else (0.0 if u == 0 else math.inf)
    if r == 0:  # (an underflow, whose sign is that of v)
        same = u == 0 and math.copysign(1, u) == math.copysign(1, r)
        return 0.0 if same else float(abs(u - v) / abs(v))
    return abs(u - r) / abs(r)


def classify(eb: float, eh: float) -> Kind:
    if eh > max(WORSE * eb, FLOOR):
        return Kind.REGRESSED
    if eb > max(WORSE * eh, FLOOR):
        return Kind.IMPROVED
    return Kind.BOTH_WRONG if min(eb, eh) > FLOOR else Kind.BOTH_RIGHT


def show(y: Value) -> str:
    return y if isinstance(y, str) else " ".join(map(repr, _parts(y)))


class Verdict(NamedTuple):
    kind: Kind
    ref: list[str]
    eb: float
    eh: float


def check(
    refs: Refs,
    spec: specs.Spec,
    cases: list[Case],
    ys: list[list[Value]],
    todo: list[int],
) -> dict[int, Verdict]:
    """Compare the results of the cases `todo` with their references."""
    out: dict[int, Verdict] = {}
    tasks = [(spec.name, cases[i].args) for i in todo]
    for i, ref in zip(todo, refs.map(tasks), strict=True):
        if isinstance(ref, str):
            out[i] = Verdict(Kind.NO_REF, [], math.nan, math.nan)
            continue
        # (only the parts that differ, e.g. not the real part if only the sign of a zero
        # imaginary part does)
        x, y = ys[0][i], ys[1][i]
        parts = [d > 0 for d in diffs(x, y)]
        eb, eh = error(x, ref, parts), error(y, ref, parts)
        out[i] = Verdict(classify(eb, eh), ref, eb, eh)
    return out


def fuzz(
    spec: specs.Spec,
    opts: argparse.Namespace,
    run: Callable[[list[Case]], list[list[Value]]],
    refs: Refs,
) -> tuple[list[str | int], list[str]] | None:
    """The row of the table and the TSV lines of the regressions, or None if a revision lacks
    the function."""
    rng = random.Random(f"{opts.seed}:{spec.name}")
    fuzzers = [Fuzzer(spec, dom, rng) for dom in spec.domains]
    cases = [Case(fuzzers[0], spec.domains[0].sample(rng))]
    ys = run(cases)
    if "unsupported" in (ys[0][0], ys[1][0]):
        return None

    def diff(i: int) -> float:
        return max(diffs(ys[0][i], ys[1][i]))

    def failed(i: int) -> bool:
        return isinstance(ys[0][i], str) != isinstance(ys[1][i], str)

    seen: set[tuple[str, bytes]] = set()
    flagged: list[int] = []
    checked: dict[int, Verdict] = {}
    for r in range(NEAR_ROUNDS + 1):
        # mutants of the regressions, and of the flagged inputs, but not of those that only
        # differ in the sign of a zero, which can be many
        worse = [cases[i] for i, v in checked.items() if v.kind == Kind.REGRESSED]
        other = [cases[i] for i in flagged if diff(i) > RTOL]
        seeds = rng.choices(worse, k=opts.cases // 2) if worse else []
        seeds += rng.choices(other, k=opts.cases - len(seeds)) if other else []
        if not seeds:  # (nothing to mutate yet)
            seeds = [f.sample(rng) for f in rng.choices(fuzzers, k=opts.cases)]
        more = [c.fuzzer.near(rng, c) for c in seeds] if r else seeds
        done = len(cases) if r else 0
        cases += more
        ys = [y + z for y, z in zip(ys, run(more), strict=True)]
        for i in range(done, len(cases)):
            c = cases[i]
            key = (c.fuzzer.name, struct.pack(f"{len(c.args)}d", *c.args))
            if key not in seen:
                seen.add(key)
                if diff(i):
                    flagged.append(i)
        # a panic or timeout of one revision needs no reference; then a share of the others,
        # those that differ in kind first
        for i in flagged:
            if failed(i) and i not in checked:
                eb, eh = (
                    (math.inf, 0.0) if isinstance(ys[0][i], str) else (0.0, math.inf)
                )
                checked[i] = Verdict(classify(eb, eh), [], eb, eh)
        todo = [i for i in flagged if i not in checked]
        rng.shuffle(todo)
        todo.sort(key=lambda i: diff(i) < math.inf)
        share = max(opts.max_refs - len(checked), 0) // (NEAR_ROUNDS + 1 - r)
        checked |= check(refs, spec, cases, ys, todo[:share])

    n = Counter(v.kind for v in checked.values())
    row = [spec.name, len(seen), len(flagged), *(n[k] for k in Kind)]
    row.append(len(flagged) - len(checked))
    worst = sorted(
        ((i, v) for i, v in checked.items() if v.kind == Kind.REGRESSED),
        key=lambda iv: -iv[1].eh,
    )
    lines = [
        "\t".join(
            [
                spec.name,
                " ".join(map(repr, cases[i].args)),
                show(ys[0][i]),
                show(ys[1][i]),
                " ".join(mp.nstr(mp.mpf(p), 17) for p in v.ref),
                f"{v.eb:.1e}",
                f"{v.eh:.1e}",
            ]
        )
        for i, v in worst
    ]
    return row, lines


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None
    )
    ap.add_argument(
        "--base",
        default="origin/master",
        help="the base revision (default: %(default)s)",
    )
    ap.add_argument("--head", help="the revision to check (default: the working tree)")
    ap.add_argument("--only", nargs="+", help="spec or function names (default: all)")
    ap.add_argument(
        "--cases",
        type=int,
        default=20000,
        help="inputs per spec and round (default: %(default)s)",
    )
    ap.add_argument(
        "--seed", default="0", help="the seed of the inputs (default: %(default)s)"
    )
    ap.add_argument(
        "--max-refs",
        type=int,
        default=2000,
        help="references per spec (default: %(default)s)",
    )
    ap.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="reference worker processes, at least 1 (default: %(default)s)",
    )
    opts = ap.parse_args()
    if opts.jobs < 1:
        ap.error("--jobs must be at least 1")
    mp.mp.dps = specs.DPS
    only = {n for arg in opts.only or [] for n in arg.split(",")}
    chosen = [s for s in specs.SPECS if not only or only & {s.name, s.func}]
    if unknown := only - {n for s in chosen for n in (s.name, s.func)}:
        ap.error(f"unknown spec or function: {', '.join(sorted(unknown))}")
    CACHE.mkdir(parents=True, exist_ok=True)
    try:
        exes = [program(opts.base), program(opts.head)]
    except subprocess.CalledProcessError as e:
        ap.exit(
            2, f"error: building the revisions failed (exit status {e.returncode})\n"
        )
    threads = concurrent.futures.ThreadPoolExecutor(len(exes))
    budgets = [Budget() for _ in exes]

    def run(cases: list[Case]) -> list[list[Value]]:
        lines = [(c.fuzzer.name, c.args) for c in cases]
        return list(threads.map(evaluate, exes, [lines] * len(exes), budgets))

    rows: list[list[str | int]] = []
    regressions: list[str] = []
    refs = Refs(opts.jobs)
    for spec in chosen:
        for b in budgets:
            b.left = TIMEOUT_BUDGET
        if (res := fuzz(spec, opts, run, refs)) is None:
            accuracy.log(f"{spec.name}: skipped, since a revision lacks it")
            continue
        rows.append(res[0])
        regressions += res[1]
        accuracy.log(f"{spec.name}: {res[0][COLUMNS.index(Kind.REGRESSED)]} regressed")

    widths = [
        max(len(str(r[j])) for r in [COLUMNS, *rows]) for j in range(len(COLUMNS))
    ]
    for r in [COLUMNS, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(r, widths, strict=True)))
    head = "function\tinput\tbase\thead\treference\tbase error\thead error"
    TSV.write_text("".join(f"{line}\n" for line in [head, *regressions]))
    print(f"\n{len(regressions)} regressions, in {TSV}")
    for line in [head, *regressions[:10]]:
        print(line.replace("\t", "  "))
    sys.exit(1 if regressions else 0)


if __name__ == "__main__":
    main()
