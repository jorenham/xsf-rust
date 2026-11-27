"""Compare the accuracy of the working tree (or `--head`) with a base revision over whole domains.

Usage (from the repository root; or `just domain-diff`):

    uv run --project tools tools/accuracy/domain_diff.py [--base REV] [--head REV]
        [--only NAME ...] [--points N] [--max-refs N] [--seed S] [--jobs N]

Steps:

1. Sample each domain of the specs of `specs.py` (or of `--only`, by spec or function name) at
   `--points` scrambled Sobol points (rounded up to a power of 2), which its distributions map
   into the domain, so log-uniformly where they are log-uniform; repeated inputs count once.
   A sampler that rejects draws or takes a variable number of them (e.g. of an array) runs out
   of coordinates of its point and draws pseudo-randomly instead; the share of such points is
   in the column "pseudo-random".
2. Evaluate both revisions with the dump programs of `fuzz_diff.py`, and keep the inputs whose
   results differ in their bits (with all NaNs equal).
3. Compute the references (as `fuzz_diff.ref`) of all differing inputs where a finite result
   became NaN or the reverse, and of a random sample of at most `--max-refs` of the others;
   they are cached by input in `.cache/domain-diff/`.
4. Print a markdown table with a row per domain, and write the flagged inputs (those that got
   worse, fail, panic, or whose reference failed) to `.cache/domain-diff.tsv`.

Of the sample of the others, an input got better or worse if one error is more than twice the
other and above 4 units of roundoff (`fuzz_diff.verdict`), which ignores last-bit changes, and
one of its errors is above `ACCURATE`; one that got worse with both below it is a flag instead.
Inputs of the sample whose reference failed are not counted; their number is in the column "refs".
The log10 ratios (of the errors clamped to [2^-53, 1]), the share of errors above `ACCURATE`,
and the gmean and peak (as in `docs/accuracy.md`) are also those of this sample.

Exit with 1 if a domain has

- (a) a finite result that is now NaN (or a panic), with a base error of at most `NAN_BASE`,
- (b) a NaN (not a panic) that is now finite, with an error above `NAN_HEAD`,
- (c) a net loss: `worse - better > SIGMAS * sqrt(worse + better)`; a heuristic of a sign test
  with more than `SIGMAS` standard deviations, so that a loss in one domain fails whatever the
  gains in the others, but a few losses among many gains do not (the inputs of one Sobol net
  are not independent, and failed references drop votes),
- (d) a new panic,

or if the head times out where the base does not.

Exit with 2 (and 1 takes precedence) if nothing failed, but a comparison could not be made: a
spec was skipped (an unsupported or timed out revision, other than a head-only timeout), or the
reference of an input of (a) or (b) failed (marked "?" in "fails").
"""

import argparse
import cmath
import fcntl
import hashlib
import inspect
import json
import math
import multiprocessing
import os
import random
import sys
import zlib
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Final, Self, override

import flint
import fuzz_diff
import mpmath as mp
import specs
from fuzz_diff import Line, Value
from scipy.stats import qmc

import accuracy

NAN_BASE: Final = 1e-3
NAN_HEAD: Final = 1e-2
SIGMAS: Final = 3
ACCURATE: Final = 1e-13
CACHE: Final = accuracy.CACHE / "domain-diff"
TSV: Final = accuracy.CACHE / "domain-diff.tsv"
COLUMNS: Final = (
    "function",
    "domain",
    "differ / inputs",
    "pseudo-random",
    "refs",
    "better",
    "worse",
    "flags",
    "log10 ratio p5 / p50 / p95",
    f"> {ACCURATE:.0e}",
    "gmean",
    "peak",
    "finite -> NaN",
    "NaN -> finite",
    "panics",
    "fails",
)


class Quasi(random.Random):
    """A `Random` whose `random()` returns the coordinates of the point of `at` in turn, then
    pseudo-random numbers, so that `Domain.sample` maps the point into the domain."""

    # (`Domain.sample` must draw only through `random()` and `_randbelow`, as `uniform`,
    # `choice` and `randint` do)

    def __init__(self, seed: int) -> None:
        super().__init__(seed)
        self.coords: Iterator[float] = iter(())
        self.fell = False  # whether `random()` ran out of coordinates since `at`

    def at(self, point: Iterable[float]) -> Self:
        self.coords = iter(point)
        self.fell = False
        return self

    @override
    def random(self) -> float:
        u = next(self.coords, None)
        self.fell |= u is None
        return super().random() if u is None else u

    def _randbelow(self, n: int) -> int:
        # (monotone in `random()`, unlike the default of a subclass without `getrandbits`)
        return min(int(self.random() * n), n - 1)


def key(args: accuracy.Args) -> str:
    """The inputs as the dump programs read them, which keeps the signs of zeros."""
    return " ".join(map(repr, args))


def sample(
    spec: specs.Spec, dom: specs.Domain, n: int, seed: int
) -> tuple[list[accuracy.Args], float]:
    """The distinct inputs at `n` Sobol points (rounded up to a power of 2), and the share of
    the points that needed pseudo-random draws."""
    # (a value and a sign per argument, as for `LU(..., sign=0)`)
    sobol = qmc.Sobol(2 * len(spec.args), rng=seed)
    rng = Quasi(seed)
    points = sobol.random_base2((n - 1).bit_length())
    xs: list[accuracy.Args] = []
    fell = 0
    for p in points:
        xs.append(dom.sample(rng.at(p.tolist())))
        fell += rng.fell
    return list({key(a): a for a in xs}.values()), fell / len(points)


def bits(y: Value) -> tuple[bytes | str, ...] | None:
    return None if y is None else accuracy.bits(y)


def finite(y: Value) -> bool:
    return y is not None and cmath.isfinite(y)


def nan(y: Value) -> bool:
    return y is not None and cmath.isnan(y)


def lost(b: Value, h: Value) -> bool:
    """Rule (a): a finite result that is now NaN or a panic."""
    return finite(b) and (h is None or nan(h))


def gained(b: Value, h: Value) -> bool:
    """Rule (b): a NaN (not a panic) that is now finite."""
    return nan(b) and finite(h)


def error(y: Value, r: accuracy.Ref) -> float:
    return math.inf if y is None else accuracy.error(y, r)


def log10_ratio(eb: float, eh: float) -> float:
    def clamp(e: float) -> float:
        return min(max(e, accuracy.GMEAN_FLOOR), 1.0)

    return math.log10(clamp(eh) / clamp(eb))


def load(path: Path) -> dict[str, accuracy.Ref]:
    try:
        return json.loads(path.read_text())
    except OSError, json.JSONDecodeError:
        return {}


def references(
    spec: specs.Spec, args: list[accuracy.Args], pool: Pool
) -> list[accuracy.Ref | str]:
    """The references of `args`, cached by input, in a file of the spec keyed by the sources
    of the references and the versions."""
    arb = fuzz_diff.ARB.get(spec.name)
    # (the name of a method of flint, whose repr has its address)
    arb_src = getattr(arb, "__qualname__", None)
    if inspect.isfunction(arb):
        arb_src = inspect.getsource(arb)
    src = [
        accuracy.ref_key(spec, []),
        inspect.getsource(accuracy.ref_task),
        inspect.getsource(fuzz_diff),
        arb_src,
        flint.__version__,
        spec.complex,
        # (also the globals that the key of the reference misses)
        hashlib.sha1(Path(specs.__file__).read_bytes()).hexdigest(),
    ]
    digest = hashlib.sha1(json.dumps(src).encode()).hexdigest()
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{spec.name}-{digest}.json"
    cache = load(path)
    keys = [key(a) for a in args]
    todo = {k: a for k, a in zip(keys, args, strict=True) if k not in cache}
    new: dict[str, accuracy.Ref | str] = {}
    if todo:
        accuracy.log(f"{spec.name}: computing {len(todo)} references")
        tasks = [(spec.name, a) for a in todo.values()]
        new = dict(zip(todo, pool.map(fuzz_diff.ref, tasks, chunksize=1), strict=True))
        # (merged with what concurrent runs wrote, and renamed into place, under a lock)
        with path.with_suffix(".lock").open("w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            merged = load(path) | {k: r for k, r in new.items() if isinstance(r, list)}
            part = path.with_suffix(f".{os.getpid()}")
            part.write_text(json.dumps(merged))
            part.replace(path)
    return [cache[k] if k in cache else new[k] for k in keys]


def fmt(x: float | None) -> str:
    return "-" if x is None else f"{x:.1e}"


def percentiles(xs: list[float]) -> str:
    """The 5th, 50th and 95th percentiles (nearest rank), or "-"."""
    xs = sorted(xs)
    ps = (xs[min(int(q * len(xs)), len(xs) - 1)] for q in (0.05, 0.5, 0.95))
    return " / ".join(f"{p:+.1f}" for p in ps) if xs else "-"


def side(es: list[float]) -> tuple[str, str, str]:
    """The share of the errors above `ACCURATE`, their gmean, and their peak."""
    st = accuracy.stats(es)
    if st is None:
        return "-", "-", "-"
    peak, gmean, _ = st
    return f"{sum(e > ACCURATE for e in es) / len(es):.0%}", fmt(gmean), fmt(peak)


def moved(ids: list[int], fails: list[int], noref: set[int]) -> str:
    """The cell of the inputs whose results became NaN or finite."""
    counts = ((len(fails), "fail"), (sum(i in noref for i in ids), "no ref"))
    extra = [f"{n} {s}" for n, s in counts if n]
    return f"{len(ids)} ({', '.join(extra)})" if extra else str(len(ids))


def summarize(
    n: int,
    fallback: float,
    diff: list[int],
    rest: list[int],
    errs: dict[int, tuple[float, float]],
    noref: set[int],
    b: list[Value],
    h: list[Value],
) -> tuple[list[str], str, dict[int, list[str]], bool]:
    """The cells of a domain, its failures, and its flagged inputs (with notes), from its
    number of inputs, the share of pseudo-random points, its differing inputs (as indices), the
    sample of those whose results are both NaN or neither, the base and head errors of the
    inputs with a reference, and the inputs whose reference failed. Also whether a reference
    of (a) or (b) is missing."""
    pairs = {i: errs[i] for i in rest if i in errs}
    kinds = {i: fuzz_diff.verdict(*e) for i, e in pairs.items()}
    minor = {i for i, e in pairs.items() if max(e) < ACCURATE}
    flags = [i for i in minor if kinds[i] == "regressed"]
    worse = [i for i, k in kinds.items() if k == "regressed" and i not in minor]
    better = sum(k == "improved" and i not in minor for i, k in kinds.items())
    loss = len(worse) - better > SIGMAS * math.sqrt(len(worse) + better)
    sides = [side([e[j] for e in pairs.values()]) for j in (0, 1)]

    to_nan = [i for i in diff if lost(b[i], h[i])]
    to_fin = [i for i in diff if gained(b[i], h[i])]
    fail_a = [i for i in to_nan if i in errs and errs[i][0] <= NAN_BASE]
    fail_b = [i for i in to_fin if i in errs and errs[i][1] > NAN_HEAD]
    panics = [i for i in diff if h[i] is None and b[i] is not None]
    unknown = any(i in noref for i in [*to_nan, *to_fin])
    bad = (fail_a, fail_b, loss, panics)
    fails = "".join(c for c, x in zip("abcd", bad, strict=True) if x)

    lost_refs = sum(i in noref for i in rest)
    refs = f"{len(pairs)} ({lost_refs} fail)" if lost_refs else str(len(pairs))
    cells = [
        f"{len(diff)} / {n}",
        f"{fallback:.0%}",
        refs,
        str(better),
        str(len(worse)),
    ]
    cells += [str(len(flags)), percentiles([log10_ratio(*e) for e in pairs.values()])]
    cells += [" -> ".join(s) for s in zip(*sides)] if pairs else ["-"] * 3
    cells += [
        moved(to_nan, fail_a, noref),
        moved(to_fin, fail_b, noref),
        str(len(panics)),
        fails + "?" * unknown,
    ]
    notes: dict[int, list[str]] = {}
    labels = [(worse, "worse"), (flags, "flag"), (fail_a, "a"), (fail_b, "b")]
    labels += [
        (panics, "panic"),
        ([i for i in [*to_nan, *to_fin, *rest] if i in noref], "no ref"),
    ]
    for ids, label in labels:
        for i in ids:
            notes.setdefault(i, []).append(label)
    return cells, fails, notes, unknown


def compare(
    spec: specs.Spec, exes: list[Path], pool: Pool, opts: argparse.Namespace
) -> tuple[list[list[str]], list[tuple[float, str]], bool, bool]:
    """The rows of the domains of `spec`, its flagged inputs (as TSV lines, with their head
    errors), whether it fails, and whether it is inconclusive."""
    descs = [d.desc(spec.args) for d in spec.domains]
    seeds = [zlib.crc32(f"{opts.seed}:{spec.name}: {d}".encode()) for d in descs]
    spans: list[range] = []
    fallbacks: list[float] = []
    lines: list[Line] = []
    for dom, seed in zip(spec.domains, seeds, strict=True):
        xs, fallback = sample(spec, dom, opts.points, seed)
        spans.append(range(len(lines), len(lines) + len(xs)))
        fallbacks.append(fallback)
        lines += [(dom.call or spec.func, a) for a in xs]
    with ThreadPoolExecutor() as threads:
        b, h = threads.map(fuzz_diff.evaluate, exes, [lines] * len(exes))
    if isinstance(b, str) or isinstance(h, str):
        status = [y if isinstance(y, str) else "ok" for y in (b, h)]
        accuracy.log(f"{spec.name}: skipped, base {status[0]}, head {status[1]}")
        head_only = status == ["ok", "timeout"]
        note = (
            "timeout" if head_only else f"skipped: base {status[0]}, head {status[1]}"
        )
        row = [spec.name, "-", *["-"] * (len(COLUMNS) - 3), note]
        return (
            [row],
            [(math.inf, f"{spec.name}{'\t' * 8}{note}")],
            head_only,
            not head_only,
        )

    # per domain, the differing inputs, the inputs that became NaN or finite, and a sample of
    # the others
    diffs = [[i for i in span if bits(b[i]) != bits(h[i])] for span in spans]
    samples: list[list[int]] = []
    todo: set[int] = set()
    for diff, seed in zip(diffs, seeds, strict=True):
        rest = [i for i in diff if not (lost(b[i], h[i]) or gained(b[i], h[i]))]
        todo |= set(diff) - set(rest)
        sampled = random.Random(seed).sample(rest, min(opts.max_refs, len(rest)))
        samples.append(sampled)
        todo |= set(sampled)
    todo_sorted = sorted(todo)
    rr = dict(
        zip(
            todo_sorted,
            references(spec, [lines[i][1] for i in todo_sorted], pool),
            strict=True,
        )
    )
    errs = {
        i: (error(b[i], r), error(h[i], r))
        for i, r in rr.items()
        if isinstance(r, list)
    }
    noref = {i for i, r in rr.items() if isinstance(r, str)}

    rows: list[list[str]] = []
    flagged: list[tuple[float, str]] = []
    failed = inconclusive = False
    for desc, span, diff, rest, fb in zip(
        descs, spans, diffs, samples, fallbacks, strict=True
    ):
        cells, fails, notes, unknown = summarize(
            len(span), fb, diff, rest, errs, noref, b, h
        )
        rows.append([spec.name, desc, *cells])
        failed |= bool(fails)
        inconclusive |= unknown
        for i, labels in notes.items():
            show = ["panic" if y is None else repr(y) for y in (b[i], h[i])]
            ref = json.dumps(rr[i][0]) if i in errs else "-"
            eb, eh = errs.get(i, (math.nan, math.inf))
            cols = [spec.name, desc, key(lines[i][1]), *show, ref, fmt(eb), fmt(eh)]
            flagged.append((eh, "\t".join([*cols, ", ".join(labels)])))
    return rows, flagged, failed, inconclusive


def positive(s: str) -> int:
    n = int(s)
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1: {s}")
    return n


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None
    )
    ap.add_argument("--base", default="origin/master", help="(default: %(default)s)")
    ap.add_argument("--head", help="the revision to check (default: the working tree)")
    ap.add_argument("--only", nargs="+", default=[], help="spec or function names")
    ap.add_argument(
        "--points",
        type=positive,
        default=2**14,
        help="inputs per domain (default: %(default)s)",
    )
    ap.add_argument(
        "--max-refs",
        type=positive,
        default=500,
        help="references per domain and kind (default: %(default)s)",
    )
    ap.add_argument("--seed", default="0", help="the seed of the inputs")
    ap.add_argument(
        "--jobs", type=positive, default=2, help="reference worker processes"
    )
    opts = ap.parse_args()
    only = set(opts.only)
    chosen = [s for s in specs.SPECS if not only or only & {s.name, s.func}]
    if unknown := only - {n for s in chosen for n in (s.name, s.func)}:
        ap.error(f"unknown spec or function: {', '.join(sorted(unknown))}")
    fuzz_diff.CACHE.mkdir(parents=True, exist_ok=True)
    exes = [fuzz_diff.program(opts.base), fuzz_diff.program(opts.head)]
    mp.mp.dps = specs.DPS

    print("|", " | ".join(COLUMNS), "|")
    print("|", " | ".join("---" for _ in COLUMNS), "|")
    flagged: list[tuple[float, str]] = []
    failed = inconclusive = False
    with multiprocessing.Pool(opts.jobs, maxtasksperchild=500) as pool:
        for spec in chosen:
            rows, f, bad, unknown = compare(spec, exes, pool, opts)
            for row in rows:
                print("|", " | ".join(row), "|", flush=True)
            flagged += f
            failed |= bad
            inconclusive |= unknown

    head = (
        "function\tdomain\tinput\tbase\thead\treference\tbase error\thead error\tnote"
    )
    lines = [line for _, line in sorted(flagged, reverse=True)]
    TSV.write_text("".join(f"{line}\n" for line in [head, *lines]))
    code = 1 if failed else 2 if inconclusive else 0
    verdict = ("ok", "FAIL", "INCONCLUSIVE")[code]
    print(f"\n{len(lines)} inputs flagged, in {TSV}; {verdict}")
    sys.exit(code)


if __name__ == "__main__":
    main()
