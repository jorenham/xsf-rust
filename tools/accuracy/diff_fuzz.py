"""Find the inputs where the working tree (or `--head`) is less accurate than a base revision.

Usage (from the repository root; or `just diff-fuzz`):

    uv run --project tools tools/accuracy/diff_fuzz.py [--base REV] [--head REV]
        [--only NAME ...] [--cases N] [--seed S] [--max-refs N] [--jobs N]

Steps:

1. Build the dump program of `accuracy.py` for both revisions; a revision other than the working
   tree is exported with `git archive`, and its program is cached in `.cache/diff-fuzz/`.
2. Evaluate both on the same `--cases` seeded inputs per spec of `specs.py` (or of `--only`, by
   spec or function name): samples of its domains, with the parameters scaled log-uniformly or
   rounded to (half-)integers, single arguments moved to log-uniform, (half-)integer, boundary
   or special values, a parameter or real argument set to (nearly) a signed sum of others plus
   an integer, and the imaginary parts of complex arguments set to 0 or tiny values of either
   sign.
3. Flag the inputs where the results differ by more than `RTOL` relative to the larger, in a part
   of a complex result, or where one is NaN, infinite, zero, or a panic and the other is not, or
   where the signs of zeros differ.
4. Compute references for a share of the flagged inputs (at most `--max-refs` in all; those that
   differ in kind first, then spread over the signs and magnitudes of the arguments), with Arb
   where `specs.ARB` or `ARB` has the function, else with the reference of the spec. Classify
   each by the largest relative error of the parts that differ: regressed (more than `WORSE`
   times the base error, and more than `FLOOR`), improved (vice versa), both wrong, or both
   right.
5. For `NEAR_ROUNDS` rounds, repeat 3 and 4 for `--cases` mutants of the regressions and the
   flagged inputs so far, with the parameters scaled (keeping that sum), two arguments shifted
   together, or one moved slightly or as in 2 (or, if none are flagged yet, for new inputs as
   in 2).
6. Print the counts, write the regressions to `.cache/diff-fuzz.tsv`, and exit with 1 if there
   are any.
"""

import argparse
import concurrent.futures
import io
import json
import math
import multiprocessing
import os
import random
import shutil
import signal
import subprocess
import sys
import tarfile
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Final, NamedTuple

import flint
import mpmath as mp
import specs

import accuracy

RTOL: Final = 1e-12
WORSE: Final = 2.0
FLOOR: Final = 4e-16
# a part of a reference below this share of its modulus counts as 0, so that its error is
# relative to the modulus
TINY_PART: Final = 2.0**-200
ARB_PRECS: Final = (128, 256, 512, 1024, 2048)
# seconds per reference, so that a few slow ones cannot stall a run
REF_TIMEOUT: Final = 10
CACHE: Final = accuracy.CACHE / "diff-fuzz"
TSV: Final = accuracy.CACHE / "diff-fuzz.tsv"
# lines per run of a dump program; a failed run is bisected down to the panicking inputs
CHUNK: Final = 1000
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

type Args = tuple[float, ...]
type Value = Args | None  # None for a panic or a timeout


def build(manifest: Path, target: Path | None = None) -> Path:
    """Build a dump program in release mode; returns its path."""
    cmd = ["cargo", "build", "--release", "--manifest-path", str(manifest)]
    cmd.append("--message-format=json-render-diagnostics")
    env = os.environ | ({} if target is None else {"CARGO_TARGET_DIR": str(target)})
    out = subprocess.run(
        cmd, capture_output=True, text=True, check=True, env=env
    ).stdout
    exes = [m.get("executable") for m in map(json.loads, out.splitlines())]
    return Path([e for e in exes if e][-1])


def program(rev: str | None) -> Path:
    """The dump program of the working tree (`rev` None), or of a commit (cached)."""
    if rev is None:
        return build(accuracy.HERE / "Cargo.toml")
    git = ["git", "-C", str(accuracy.ROOT)]
    sha = subprocess.check_output(
        [*git, "rev-parse", f"{rev}^{{commit}}"], text=True
    ).strip()[:12]
    exe = CACHE / f"xsf-accuracy-{sha}"
    if not exe.exists():
        accuracy.log(f"building {rev} ({sha})")
        # (a source directory per commit, since cargo would take the old mtimes of the
        # exported files for a fresh build of another commit, and a target directory of its
        # own, since the program of the working tree has the same name)
        src = CACHE / f"src-{sha}"
        tar = subprocess.check_output([*git, "archive", sha])
        with tarfile.open(fileobj=io.BytesIO(tar)) as tf:
            tf.extractall(src, filter="data")
        shutil.copy2(
            build(src / "tools" / "accuracy" / "Cargo.toml", CACHE / "target"), exe
        )
        shutil.rmtree(src)
    return exe


def evaluate(exe: Path, lines: list[str]) -> list[Value]:
    """The results of the input lines, or None where the program fails on them."""
    if len(lines) > CHUNK:
        return [
            v
            for i in range(0, len(lines), CHUNK)
            for v in evaluate(exe, lines[i : i + CHUNK])
        ]
    try:
        out = subprocess.run(
            [str(exe), "xsf"],
            input="".join(lines),
            capture_output=True,
            text=True,
            check=True,
            timeout=10 + 0.02 * len(lines),
        ).stdout
        return [tuple(map(float, line.split())) for line in out.splitlines()]
    except subprocess.CalledProcessError, subprocess.TimeoutExpired:
        if len(lines) == 1:
            return [None]
        mid = len(lines) // 2
        return evaluate(exe, lines[:mid]) + evaluate(exe, lines[mid:])


@dataclass(frozen=True)
class Case:
    fuzzer: Fuzzer
    args: Args
    # a parameter i set to sum(sign * args[j]) + offset, as (i, [(sign, j), ...], offset)
    rel: tuple[int, list[tuple[int, int]], float] | None = None

    @property
    def line(self) -> str:
        return f"{self.fuzzer.name} {' '.join(map(repr, self.args))}\n"

    def related(self, a: list[float]) -> Case:
        """The case of the arguments `a`, with the parameter of `rel` set again."""
        if self.rel is not None:
            i, terms, offset = self.rel
            a[i] = sum(s * a[j] for s, j in terms) + offset
        return Case(self.fuzzer, tuple(a), self.rel)


class Fuzzer:
    """Mutated samples of a domain of a spec."""

    def __init__(self, spec: specs.Spec, dom: specs.Domain, rng: random.Random) -> None:
        self.dom, self.name = dom, dom.call or spec.func
        # the integer arguments (by their distribution, or for a sampler, by its samples) stay
        samples = [dom.sample(rng) for _ in range(20)]
        n = len(samples[0])
        dists = dom.dists or [None] * n
        self.free = [
            i
            for i, d in enumerate(dists)
            if not isinstance(d, specs.Int | specs.ILU | specs.Each)
            and not all(s[i] == int(s[i]) for s in samples)
        ]
        names = spec.args
        self.imag = [i for i in self.free if names[i] in IMAG]
        # the parameters: all but the complex arguments, or else the last argument
        last = -1 if spec.complex else n - 1
        self.params = [
            i
            for i in self.free
            if i != last and names[i] not in IMAG and i + 1 not in self.imag
        ]
        # what a relation links: the parameters, and a real argument
        self.linked = self.params + [i for i in self.free if i == last]
        self.bounds = {i: [0.0, 1.0, -1.0, *_bounds(d)] for i, d in enumerate(dists)}
        self.top = {i: max(map(abs, self.bounds[i])) for i in self.free}

    def sample(self, rng: random.Random) -> Case:
        a = list(self.dom.sample(rng))
        if rng.random() < 0.3:
            self.scale(rng, a)
        if (r := rng.random()) < 0.3:
            k = 1 if r < 0.2 else 2
            for i in self.params:
                a[i] = round(k * a[i]) / k
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
        for _ in range(rng.randint(1, 2)):
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
                mag = 10 ** rng.uniform(-8, math.log10(max(self.top[i], 1.0)) + 2)
                return rng.choice((-1, 1)) * mag
            case 1:
                return _nudge(rng, round(x))
            case 2:
                return _nudge(rng, round(2 * x) / 2)
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


def diffs(x: Value, y: Value) -> list[float]:
    """How much the parts of two results differ: 0 if by at most `RTOL` relative to the larger,
    `RTOL` for zeros of different signs, else their relative difference, or infinity if they
    differ in kind (NaN, infinite, zero, or a panic)."""
    if x is None or y is None:
        return [0.0 if x == y else math.inf]
    return list(map(_diff, x, y))


def _diff(u: float, v: float) -> float:
    if math.isnan(u) or math.isnan(v):
        return 0.0 if math.isnan(u) and math.isnan(v) else math.inf
    if u == v:
        return 0.0 if math.copysign(1, u) == math.copysign(1, v) else RTOL
    if not (math.isfinite(u) and math.isfinite(v)) or u == 0 or v == 0:
        return math.inf
    d = abs(u - v) / max(abs(u), abs(v))
    return d if d > RTOL else 0.0


def _arb(name: str, args: Args) -> list[str] | None:
    """The Arb value at the lowest of `ARB_PRECS` where each part is accurate, or tiny."""
    spec = specs.BY_NAME[name]
    # (a signed zero imaginary part as a tiny one, for the side of a branch cut)
    args_arb = [
        math.copysign(1, x) * flint.arb(2) ** -(2 * ARB_PRECS[-1])
        if x == 0 and arg in IMAG
        else flint.arb(x)
        for x, arg in zip(args, spec.args, strict=True)
    ]
    for prec in ARB_PRECS:
        with flint.ctx.workprec(prec):
            y = flint.acb(ARB[name](*args_arb))
        parts = [y.real, y.imag] if spec.complex else [y.real]
        if not all(p.is_finite() for p in parts):
            continue
        mod = abs(y).abs_lower()
        if all(
            p.rel_accuracy_bits() >= 64
            or (p.contains(0) and p.abs_upper() <= TINY_PART * mod)
            or p.is_zero()
            for p in parts
        ) and (mod > 0 or y.is_zero()):
            return [p.mid().str(40, radius=False) for p in parts]
    return None


def _cap_timeouts() -> None:
    for spec in specs.BY_NAME.values():
        spec.timeout = min(spec.timeout, REF_TIMEOUT)


def ref_task(task: tuple[str, Args]) -> list[str] | str:
    """The reference (its parts), or why there is none."""
    name, args = task
    if name in ARB:
        signal.signal(signal.SIGALRM, accuracy._alarm)
        signal.alarm(specs.BY_NAME[name].timeout)
        try:
            if (y := _arb(name, args)) is not None:
                return y
        except accuracy.Timeout:
            return "timeout"
        finally:
            signal.alarm(0)
    match accuracy.ref_task(task):
        case [list(v), _]:
            return v
        case [str(v), _]:
            return [v]
        case e:
            return str(e)


def error(y: Value, ref: list[str], parts: list[bool]) -> float:
    """The largest relative error of the given parts, where the parts of `ref` below `TINY_PART`
    of its modulus are relative to that."""
    if y is None:
        return math.inf
    f = [mp.mpf(v) for v in ref]
    mod = mp.sqrt(sum(v * v for v in f))
    errs = []
    for u, v, part in zip(y, f, parts, strict=True):
        if not part:
            continue
        if abs(v) <= TINY_PART * mod:
            v = mp.mpf(0)
        if math.isnan(u):
            errs.append(math.inf)
        elif abs(v) > accuracy.DBL_MAX:
            errs.append(0.0 if math.isinf(u) and (u > 0) == (v > 0) else math.inf)
        elif math.isinf(u):
            errs.append(math.inf)
        elif abs(v) < accuracy.DBL_MIN and abs(u) < accuracy.DBL_MIN:
            errs.append(0.0)
        else:
            errs.append(float(abs(u - v) / max(abs(v) or mod, accuracy.DBL_MIN)))
    return max(errs, default=0.0)


def classify(eb: float, eh: float) -> str:
    if eh > max(WORSE * eb, FLOOR):
        return "regressed"
    if eb > max(WORSE * eh, FLOOR):
        return "improved"
    return "both wrong" if min(eb, eh) > FLOOR else "both right"


def spread(rng: random.Random, todo: list[int], cases: list[Case]) -> list[int]:
    """The cases `todo` in a random order that takes one of each cell in turn: of the same
    signs and magnitudes (in steps of 16, from 2^-12 to 2^12) of the arguments."""

    def cell(x: float) -> object:
        if not math.isfinite(x) or x == 0:
            return x
        return x > 0, min(max(math.floor(math.log2(abs(x)) / 4), -3), 3)

    ranks: dict[object, int] = {}
    order = {}
    for i in rng.sample(todo, len(todo)):
        key = (cases[i].fuzzer, *map(cell, cases[i].args))
        ranks[key] = ranks.get(key, -1) + 1
        order[i] = ranks[key]
    return sorted(todo, key=order.__getitem__)


def show(y: Value) -> str:
    return "panic" if y is None else " ".join(map(repr, y))


class Verdict(NamedTuple):
    kind: str  # one of `classify`, or "no ref"
    ref: list[str]
    eb: float
    eh: float


def check(
    pool: Pool,
    spec: specs.Spec,
    cases: list[Case],
    ys: list[list[Value]],
    todo: list[int],
) -> dict[int, Verdict]:
    """Compare the results of the cases `todo` with their references."""
    refs = pool.map(ref_task, [(spec.name, cases[i].args) for i in todo], chunksize=1)
    out: dict[int, Verdict] = {}
    for i, ref in zip(todo, refs, strict=True):
        if isinstance(ref, str):
            out[i] = Verdict("no ref", [], math.nan, math.nan)
            continue
        # (only the parts that differ, e.g. not the real part if only the sign of a zero
        # imaginary part does)
        x, y = ys[0][i], ys[1][i]
        both = x is not None and y is not None
        parts = [d > 0 for d in diffs(x, y)] if both else [True] * len(ref)
        eb, eh = error(x, ref, parts), error(y, ref, parts)
        out[i] = Verdict(classify(eb, eh), ref, eb, eh)
    return out


def fuzz(
    spec: specs.Spec,
    opts: argparse.Namespace,
    run: Callable[[list[Case]], list[list[Value]]],
    pool: Pool,
) -> tuple[list[str | int], list[str]] | None:
    """The row of the table and the TSV lines of the regressions, or None if a revision fails
    on an unmutated sample (e.g. lacks the function)."""
    rng = random.Random(f"{opts.seed}:{spec.name}")
    fuzzers = [Fuzzer(spec, dom, rng) for dom in spec.domains]
    cases = [Case(fuzzers[0], spec.domains[0].sample(rng))]
    cases += [f.sample(rng) for f in rng.choices(fuzzers, k=opts.cases)]
    ys = run(cases)
    if ys[0][0] is None or ys[1][0] is None:
        return None

    def diff(i: int) -> float:
        return max(diffs(ys[0][i], ys[1][i]))

    seen: set[Args] = set()
    flagged: list[int] = []
    checked: dict[int, Verdict] = {}
    done = 0
    for r in range(NEAR_ROUNDS + 1):
        for i in range(done, len(cases)):
            if cases[i].args not in seen:
                seen.add(cases[i].args)
                if diff(i):
                    flagged.append(i)
        done = len(cases)
        # a share of the references: for those that differ in kind first, then spread over the
        # cells of the inputs, since most can be in one region
        todo = spread(rng, [i for i in flagged if i not in checked], cases)
        todo.sort(key=lambda i: diff(i) < math.inf)
        share = (opts.max_refs - len(checked)) // (NEAR_ROUNDS + 1 - r)
        checked |= check(pool, spec, cases, ys, todo[:share])
        if r == NEAR_ROUNDS:
            break
        # mutants of the regressions, and of the flagged inputs, but not of those that only
        # differ in the sign of a zero, which can be many
        worse = [cases[i] for i, v in checked.items() if v.kind == "regressed"]
        other = [cases[i] for i in flagged if diff(i) > RTOL]
        seeds = rng.choices(worse, k=opts.cases // 2) if worse else []
        seeds += rng.choices(other, k=opts.cases - len(seeds)) if other else []
        if seeds:
            more = [c.fuzzer.near(rng, c) for c in seeds]
        else:  # (nothing to mutate yet)
            more = [f.sample(rng) for f in rng.choices(fuzzers, k=opts.cases)]
        cases += more
        ys = [y + z for y, z in zip(ys, run(more), strict=True)]

    n = Counter(v.kind for v in checked.values())
    row = [spec.name, len(seen), len(flagged), *(n[c] for c in COLUMNS[3:-1])]
    row.append(len(flagged) - len(checked))
    worst = sorted(
        ((i, v) for i, v in checked.items() if v.kind == "regressed"),
        key=lambda iv: -iv[1].eh,
    )
    lines = [
        "\t".join(
            [
                spec.name,
                show(cases[i].args),
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
    ap.add_argument("--base", default="origin/master", help="the base revision")
    ap.add_argument("--head", help="the revision to check (default: the working tree)")
    ap.add_argument("--only", nargs="+", help="spec or function names")
    ap.add_argument(
        "--cases", type=int, default=20000, help="inputs per spec and round"
    )
    ap.add_argument("--seed", default="0", help="the seed of the inputs")
    ap.add_argument("--max-refs", type=int, default=2000, help="references per spec")
    ap.add_argument("--jobs", type=int, default=8, help="reference worker processes")
    opts = ap.parse_args()
    mp.mp.dps = specs.DPS
    only = {n for arg in opts.only or [] for n in arg.split(",")}
    chosen = [s for s in specs.SPECS if not only or only & {s.name, s.func}]
    if unknown := only - {n for s in chosen for n in (s.name, s.func)}:
        sys.exit(f"unknown: {', '.join(sorted(unknown))}")
    CACHE.mkdir(parents=True, exist_ok=True)
    exes = [program(opts.base), program(opts.head)]
    threads = concurrent.futures.ThreadPoolExecutor(len(exes))

    def run(cases: list[Case]) -> list[list[Value]]:
        lines = [c.line for c in cases]
        return list(threads.map(evaluate, exes, [lines] * len(exes)))

    rows: list[list[str | int]] = []
    regressions: list[str] = []
    with multiprocessing.Pool(opts.jobs, _cap_timeouts) as pool:
        for spec in chosen:
            if (res := fuzz(spec, opts, run, pool)) is None:
                accuracy.log(f"{spec.name}: skipped, since a revision fails on it")
                continue
            rows.append(res[0])
            regressions += res[1]
            accuracy.log("  ".join(map(str, res[0])))

    widths = [
        max(len(str(r[j])) for r in [COLUMNS, *rows]) for j in range(len(COLUMNS))
    ]
    for r in [COLUMNS, *rows]:
        print("  ".join(str(v).rjust(w) for v, w in zip(r, widths, strict=True)))
    head = "function\tinput\tbase\thead\treference\tbase error\thead error"
    TSV.write_text("".join(f"{line}\n" for line in [head, *regressions]))
    print(f"\n{len(regressions)} regressions, in {TSV}")
    for line in regressions[:10]:
        print(line.replace("\t", "  "))
    sys.exit(1 if regressions else 0)


COLUMNS: Final = [
    "function",
    "cases",
    "flagged",
    "regressed",
    "improved",
    "both wrong",
    "both right",
    "no ref",
    "unchecked",
]


if __name__ == "__main__":
    main()
