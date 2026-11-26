"""Find inputs where the working tree (or `--head`) is less accurate than a base revision.

Usage (from the repository root; or `just diff-fuzz`):

    uv run --project tools tools/accuracy/diff_fuzz.py [--base REV] [--head REV]
        [--only NAME ...] [--cases N] [--seed S] [--max-refs N] [--jobs N]

Steps:

1. Build the dump program of `accuracy.py` for both revisions; that of a commit is built from
   `git archive` and cached in `.cache/diff-fuzz/`.
2. Evaluate both on the same `--cases` seeded inputs per spec of `specs.py` (or of `--only`, by
   spec or function name): samples of its domains, with arguments moved to (near) integers,
   half-integers, domain bounds or special values, one set to (near) an integer plus a signed
   sum of others, and imaginary parts set to signed zeros or tiny values. Skip a spec if a
   revision lacks its function or times out on it.
3. Flag the inputs whose results differ by more than `RTOL`, or where one revision panics.
   Repeat 2 and 3 for `NEAR_ROUNDS` rounds of `--cases` mutants of the flagged inputs (else of
   those whose results differ at all, else of new samples).
4. Compute the references of at most `--max-refs` flagged inputs per spec, with Arb where `ARB`
   has the function, else as `accuracy.py` does, and compare the errors of both revisions.
5. Print the counts, write the regressions to `.cache/diff-fuzz.tsv`, and exit with 1 if there
   are any (including a timeout of the head alone).
"""

import argparse
import math
import multiprocessing
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Final

import flint
import mpmath as mp
import specs

import accuracy
from accuracy import DBL_MAX, DBL_MIN

RTOL: Final = 1e-12
# a regression has more than `WORSE` times the error of the base, and more than `FLOOR` (4 units
# of roundoff)
WORSE: Final = 2
FLOOR: Final = 2.0**-51
ARB_PRECS: Final = (128, 512, 2048)
# a signed zero imaginary argument for Arb, for the side of a branch cut
ARB_ZERO: Final = flint.arb(2) ** (-2 * ARB_PRECS[-1])
# seconds per run of a dump program, and per line in it
RUN_TIMEOUT: Final = 10
LINE_TIMEOUT: Final = 0.02
# rounds of mutants of the flagged inputs
NEAR_ROUNDS: Final = 3
CACHE: Final = accuracy.CACHE / "diff-fuzz"
TSV: Final = accuracy.CACHE / "diff-fuzz.tsv"
IMAG: Final = {"im", "xi", "yi", "zi", "pi"}
# (of a random sign)
SPECIAL: Final = (0.0, math.ulp(0), DBL_MIN, DBL_MAX, math.inf, math.nan)
KINDS: Final = (
    "regressed",
    "improved",
    "both wrong",
    "both right",
    "no ref",
    "unchecked",
)

type Line = tuple[str, accuracy.Args]
type Value = complex | None  # None for a panic


def _hyp1f1_c(a: flint.arb, b: flint.arb, re: flint.arb, im: flint.arb) -> flint.acb:
    return flint.acb(re, im).hypgeom_1f1(a, b)


def _hyp2f1_c(
    a: flint.arb, b: flint.arb, c: flint.arb, re: flint.arb, im: flint.arb
) -> flint.acb:
    return flint.acb(re, im).hypgeom_2f1(a, b, c)


# with the complex functions that `specs.ARB` lacks, since their specs use mpmath
ARB: Final = specs.ARB | {"hyp1f1#c": _hyp1f1_c, "hyp2f1#c": _hyp2f1_c}


def build(root: Path, target: Path) -> Path:
    """Build the dump program of the source tree `root` into `target`; returns its path."""
    manifest = root / accuracy.HERE.relative_to(accuracy.ROOT) / "Cargo.toml"
    cmd = ["cargo", "build", "--release", "--quiet", "--manifest-path", str(manifest)]
    subprocess.run([*cmd, "--target-dir", str(target)], check=True)
    return target / "release" / "xsf-accuracy"


def program(rev: str | None) -> Path:
    """The dump program of the working tree (`rev` None), or of a commit."""
    if rev is None:
        return build(accuracy.ROOT, CACHE / "target")
    git = ["git", "-C", str(accuracy.ROOT)]
    sha = subprocess.check_output(
        [*git, "rev-parse", "--verify", f"{rev}^{{commit}}"], text=True
    ).strip()
    exe = CACHE / f"xsf-accuracy-{sha}"
    if not exe.exists():
        accuracy.log(f"building {rev} ({sha[:12]})")
        with tempfile.TemporaryDirectory() as tmp:
            tar = subprocess.check_output([*git, "archive", sha])
            subprocess.run(["tar", "-x", "-C", tmp], input=tar, check=True)
            # (renamed into place, for concurrent runs)
            part = exe.with_suffix(f".{os.getpid()}")
            shutil.copy2(build(Path(tmp), Path(tmp) / "target"), part)
            part.replace(exe)
    return exe


def run(exe: Path, lines: list[Line]) -> list[Value]:
    """The results of the lines, bisecting to find those that panic."""
    timeout = RUN_TIMEOUT + LINE_TIMEOUT * len(lines)
    try:
        return list[Value](accuracy.run([str(exe), "xsf"], lines, timeout))
    except subprocess.CalledProcessError as e:
        if "unknown function" in e.stderr:
            raise
        if len(lines) == 1:
            return [None]
    mid = len(lines) // 2
    return run(exe, lines[:mid]) + run(exe, lines[mid:])


def evaluate(exe: Path, lines: list[Line]) -> list[Value] | str:
    """The results of the lines, or "unsupported" or "timeout"."""
    try:
        return run(exe, lines)
    except subprocess.CalledProcessError:
        return "unsupported"
    except subprocess.TimeoutExpired:
        return "timeout"


def _marks(doms: list[specs.Domain], i: int) -> tuple[float, ...]:
    """0, 1, -1, and the bounds of the distributions of argument `i`."""
    marks = [0.0, 1.0, -1.0]
    for d in doms:
        match d.dists[i] if d.dists else None:
            case specs.U(lo, hi):
                marks += [lo, hi]
            case specs.Offset(c):
                marks.append(c)
    return tuple(marks)


def _small(rng: random.Random) -> float:
    """0, or a power of 2 in [2^-52, 1/2] of a random sign."""
    return rng.choice((0.0, rng.choice((-1, 1)) * 2.0 ** -rng.randint(1, 52)))


class Mutator:
    """Mutations of the arguments of a function, from samples of its domains."""

    def __init__(
        self, spec: specs.Spec, doms: list[specs.Domain], rng: random.Random
    ) -> None:
        samples = [d.sample(rng) for d in doms for _ in range(20)]
        # the arguments that are not integers (by their samples), their imaginary parts, the
        # others, and of those the parameters (not real parts)
        free = [i for i in range(len(spec.args)) if any(s[i] % 1 for s in samples)]
        self.imag = [i for i in free if spec.args[i] in IMAG]
        self.real = [i for i in free if i not in self.imag]
        self.params = [i for i in self.real if i + 1 not in self.imag]
        self.marks = {i: _marks(doms, i) for i in self.real}

    def __call__(self, rng: random.Random, args: accuracy.Args) -> accuracy.Args:
        a = list(args)
        if rng.random() < 0.2:
            k = rng.choice((1, 2))
            for i in self.params:
                a[i] = round(k * a[i], 0) / k
        # (the arguments more often, e.g. to a singular point)
        for i in self.real:
            if rng.random() < (0.25 if i in self.params else 0.5):
                t = rng.choice((a[i], round(a[i], 0), *self.marks[i]))
                a[i] = t + _small(rng) * max(1.0, abs(t))
        if self.real and rng.random() < 0.1:
            a[rng.choice(self.real)] = rng.choice((-1, 1)) * rng.choice(SPECIAL)
        if len(self.params) > 1 and rng.random() < 0.5:
            # e.g. a balanced c = a + b + m for hyp2f1, a2 = a1 + m, or x = v + m for Bessel
            k = min(len(self.params), rng.choice((2, 3, 3)))
            i, *js = rng.sample(self.params, k)
            a[i] = sum(rng.choice((-1, 1)) * a[j] for j in js)
            a[i] += rng.randint(-2, 2) + _small(rng)
        for i in self.imag:
            if rng.random() < 0.5:
                tiny = rng.choice((0.0, 2.0 ** -rng.randint(1, 1074)))
                a[i] = rng.choice((-1, 1)) * tiny
        return tuple(a)


def differ(x: Value, y: Value, rtol: float) -> bool:
    """Whether a part of two results differs by more than `rtol` (relative to the larger), or
    only one is a panic."""
    if x is None or y is None:
        return (x is None) != (y is None)
    x, y = complex(x), complex(y)
    return not all(
        math.isclose(u, v, rel_tol=rtol) or (math.isnan(u) and math.isnan(v))
        for u, v in ((x.real, y.real), (x.imag, y.imag))
    )


def ref(line: Line) -> accuracy.Ref | str:
    """The reference, with Arb where `ARB` has the function, else as `accuracy.ref_task`."""
    name, args = line
    spec = specs.BY_NAME[name]
    if (f := ARB.get(name)) is not None:
        xs = [
            math.copysign(1, x) * ARB_ZERO if x == 0 and a in IMAG else flint.arb(x)
            for x, a in zip(args, spec.args, strict=True)
        ]
        for prec in ARB_PRECS:
            with flint.ctx.workprec(prec):
                y = flint.acb(f(*xs))
            if y.rel_accuracy_bits() >= 64:
                re, im = (p.str(30, radius=False) for p in (y.real, y.imag))
                if spec.complex:
                    return [[re, im], None]
                return [re, None] if y.imag.contains(0) else f"complex result {y}"
    return accuracy.ref_task(line)


def verdict(eb: float, eh: float) -> str:
    if eh > max(WORSE * eb, FLOOR):
        return "regressed"
    if eb > max(WORSE * eh, FLOOR):
        return "improved"
    return "both wrong" if min(eb, eh) > FLOOR else "both right"


def fuzz(
    spec: specs.Spec, exes: list[Path], pool: Pool, opts: argparse.Namespace
) -> tuple[list[int], list[tuple[float, str]]]:
    """The counts of `spec` (none if skipped) and its regressions with their head errors."""
    rng = random.Random(f"{opts.seed}:{spec.name}")
    calls = [d.call or spec.func for d in spec.domains]
    mutate = {
        c: Mutator(spec, [d for d, e in zip(spec.domains, calls) if e == c], rng)
        for c in dict.fromkeys(calls)
    }
    n = -(-opts.cases // len(spec.domains))

    def samples() -> list[Line]:
        return [
            (c, d.sample(rng)) for c, d in zip(calls, spec.domains) for _ in range(n)
        ]

    new = samples()
    seen: set[str] = set()
    lines: list[Line] = []
    base: list[Value] = []
    head: list[Value] = []
    flagged: list[int] = []
    for _ in range(NEAR_ROUNDS + 1):
        fresh = {repr(c): c for c in ((c, mutate[c](rng, a)) for c, a in new)}
        fresh = [c for k, c in fresh.items() if k not in seen]
        seen.update(map(repr, fresh))
        with ThreadPoolExecutor() as threads:
            b, h = threads.map(evaluate, exes, [fresh] * len(exes))
        if isinstance(b, str) or isinstance(h, str):
            status = [y if isinstance(y, str) else "ok" for y in (b, h)]
            accuracy.log(f"{spec.name}: skipped, base {status[0]}, head {status[1]}")
            hung = [(math.inf, f"{spec.name}\t\t\ttimeout")]
            return [], hung if status == ["ok", "timeout"] else []
        lines += fresh
        base += b
        head += h
        pairs = list(enumerate(zip(base, head)))
        flagged = [i for i, (x, y) in pairs if differ(x, y, RTOL)]
        # then mutants of the flagged inputs, else of those with any difference, if any
        seeds = flagged or [i for i, (x, y) in pairs if differ(x, y, 0.0)]
        new = (
            rng.choices([lines[i] for i in seeds], k=opts.cases) if seeds else samples()
        )
    todo = [i for i in flagged if None not in (base[i], head[i])]
    todo = rng.sample(todo, min(len(todo), opts.max_refs))
    tasks = [(spec.name, lines[i][1]) for i in todo]
    refs = dict(zip(todo, pool.map(ref, tasks, chunksize=1), strict=True))
    count = Counter[str]()
    worse: list[tuple[float, str]] = []
    for i in flagged:
        x, y, r = base[i], head[i], refs.get(i)
        if x is None or y is None:
            eb, eh = (math.inf, 0.0) if x is None else (0.0, math.inf)
        elif isinstance(r, list):
            eb, eh = accuracy.error(x, r), accuracy.error(y, r)
        else:
            count["unchecked" if r is None else "no ref"] += 1
            continue
        count[kind := verdict(eb, eh)] += 1
        if kind == "regressed":
            v = r[0] if isinstance(r, list) else ""
            v = " ".join(v) if isinstance(v, list) else v or ""
            show = ["panic" if u is None else repr(u) for u in (x, y)]
            cols = [spec.name, " ".join(map(repr, lines[i][1])), *show, v]
            worse.append((eh, "\t".join([*cols, f"{eb:.1e}", f"{eh:.1e}"])))
    return [len(lines), len(flagged), *(count[k] for k in KINDS)], worse


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None
    )
    ap.add_argument("--base", default="origin/master", help="(default: %(default)s)")
    ap.add_argument("--head", help="the revision to check (default: the working tree)")
    ap.add_argument("--only", nargs="+", default=[], help="spec or function names")
    ap.add_argument("--cases", type=int, default=20000, help="inputs per spec")
    ap.add_argument("--seed", default="0", help="the seed of the inputs")
    ap.add_argument("--max-refs", type=int, default=2000, help="references per spec")
    ap.add_argument("--jobs", type=int, default=8, help="reference worker processes")
    opts = ap.parse_args()
    only = set(opts.only)
    chosen = [s for s in specs.SPECS if not only or only & {s.name, s.func}]
    if unknown := only - {n for s in chosen for n in (s.name, s.func)}:
        ap.error(f"unknown spec or function: {', '.join(sorted(unknown))}")
    CACHE.mkdir(parents=True, exist_ok=True)
    exes = [program(opts.base), program(opts.head)]
    mp.mp.dps = specs.DPS

    w = max(len(s.name) for s in chosen)
    print(f"{'function':<{w}}", *(f"{c:>10}" for c in ["cases", "flagged", *KINDS]))
    regressions: list[str] = []
    with multiprocessing.Pool(opts.jobs, maxtasksperchild=500) as pool:
        for spec in chosen:
            row, worse = fuzz(spec, exes, pool, opts)
            if row:
                print(f"{spec.name:<{w}}", *(f"{c:>10}" for c in row), flush=True)
            regressions += [line for _, line in sorted(worse, reverse=True)]

    head = "function\tinput\tbase\thead\treference\tbase error\thead error"
    TSV.write_text("".join(f"{line}\n" for line in [head, *regressions]))
    print(f"\n{len(regressions)} regressions, in {TSV}")
    for line in [head, *regressions[:10]] if regressions else []:
        print(line.replace("\t", "  "))
    sys.exit(1 if regressions else 0)


if __name__ == "__main__":
    main()
