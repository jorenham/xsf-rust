"""Compare the accuracy of the working tree (or `--head`) with a base revision over whole domains.

Each domain of the specs of `specs.py` is sampled at `--points` scrambled Sobol points (rounded
up to a power of 2, and lowered so that all domains have at most `--max-inputs`), which its
distributions map into the domain; repeated inputs count once. A sampler that runs out of
coordinates of its point (e.g. when it rejects draws) draws pseudo-randomly instead. With
`--affected`, only the specs that depend on the files of `src/` changed since `--base` (see
`affected.py`) are compared.

Both revisions are evaluated with the dump programs of `fuzz_diff.py`. The references (as
`fuzz_diff.ref`) of all inputs whose results differ and are of (a) or (b), and of a random sample
of at most `--max-refs` of the others, are cached by input in `.cache/domain-diff/`. Of that
sample, an input got better or worse if `fuzz_diff.verdict` says so and one of its errors is at
least `ACCURATE`. The peak errors are those of that sample. The table has a row per domain, and
`.cache/domain-diff.tsv` has the inputs that got worse, fail, panic, or whose reference failed.

Exit with 1 if a domain has

- (a) a base result that is neither NaN nor a panic (so an infinity counts), with a base error
  of at most `NAN_BASE`, and a head result that is NaN or a panic,
- (b) a base NaN (not a panic) and a head result that is neither NaN nor a panic (so a wrong
  infinity counts), with a head error above `NAN_HEAD`,
- (c) a net loss: `worse - better > SIGMAS * sqrt(worse + better)`; a heuristic of a sign test
  with more than `SIGMAS` standard deviations, so that a loss in one domain fails whatever the
  gains in the others, but a few losses among many gains do not (the inputs of one Sobol net
  are not independent),
- (d) a new panic.

Otherwise exit with 2 if a comparison could not be made: a spec was skipped (an unsupported
revision, or a timeout, which is retried once if only one revision timed out), the reference
of an input of (a) or (b) failed or is NaN (marked "?" in "fails"), or the sampled inputs
without a reference would make (c) a net loss if they all got worse (marked with their count).
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
import subprocess
import sys
import zlib
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Final, override

import affected
import flint
import fuzz_diff
import mpmath as mp
import specs
from fuzz_diff import Value
from scipy.stats import qmc

import accuracy

NAN_BASE: Final = 1e-3
NAN_HEAD: Final = 1e-2
SIGMAS: Final = 3
ACCURATE: Final = 1e-13
MIN_POINTS: Final = 2**10
CACHE: Final = accuracy.CACHE / "domain-diff"
TSV: Final = accuracy.CACHE / "domain-diff.tsv"
COLUMNS: Final = (
    "function",
    "domain",
    "differ / inputs",
    "better",
    "worse",
    "peak",
    "value -> NaN",
    "NaN -> value",
    "panics",
    "fails",
)


class Quasi(random.Random):
    """A `Random` whose `random()` returns the coordinates in `coords` in turn, then
    pseudo-random numbers, so that `Domain.sample` maps a Sobol point into the domain."""

    # (`Domain.sample` must draw only through `random()` and `_randbelow`, as `uniform`,
    # `choice` and `randint` do)

    coords: Iterator[float] = iter(())

    @override
    def random(self) -> float:
        u = next(self.coords, None)
        return super().random() if u is None else u

    def _randbelow(self, n: int) -> int:
        # (monotone in `random()`, unlike the default of a subclass without `getrandbits`)
        return min(int(self.random() * n), n - 1)


def key(args: accuracy.Args) -> str:
    """The inputs as the dump programs read them, which keeps the signs of zeros."""
    return " ".join(map(repr, args))


def sample(
    spec: specs.Spec, dom: specs.Domain, n: int, seed: int
) -> list[accuracy.Args]:
    """The distinct inputs at `n` Sobol points (rounded up to a power of 2)."""
    # (a value and a sign per argument, as for `LU(..., sign=0)`)
    points = qmc.Sobol(2 * len(spec.args), rng=seed).random_base2((n - 1).bit_length())
    rng = Quasi(seed)
    xs: dict[str, accuracy.Args] = {}
    for p in points:
        rng.coords = iter(p.tolist())
        x = dom.sample(rng)
        xs[key(x)] = x
    return list(xs.values())


def nan(y: Value) -> bool:
    return y is not None and cmath.isnan(y)


def error(y: Value, r: accuracy.Ref) -> float:
    return math.inf if y is None else accuracy.error(y, r)


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
    arb_src = (
        inspect.getsource(arb)
        if inspect.isfunction(arb)
        else getattr(arb, "__qualname__", None)
    )
    src = [
        accuracy.ref_key(spec, []),
        inspect.getsource(accuracy.ref_task),
        inspect.getsource(fuzz_diff),
        arb_src,
        flint.__version__,
        spec.complex,
        # (the contents of specs.py, for globals that the reference key misses)
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


def net_loss(worse: int, better: int) -> bool:
    return worse - better > SIGMAS * math.sqrt(worse + better)


def fmt(x: float | None) -> str:
    return "-" if x is None else f"{x:.1e}"


def compare(
    spec: specs.Spec, exes: list[Path], pool: Pool, opts: argparse.Namespace
) -> tuple[list[list[str]], list[tuple[float, str]], bool, bool]:
    """The rows of the domains of `spec`, its flagged inputs (as TSV lines, with their head
    errors), whether it fails, and whether it is inconclusive."""
    lines: list[fuzz_diff.Line] = []
    spans: list[range] = []
    seeds: list[int] = []
    descs = [d.desc(spec.args) for d in spec.domains]
    for dom, desc in zip(spec.domains, descs, strict=True):
        seeds.append(zlib.crc32(f"{opts.seed}:{spec.name}: {desc}".encode()))
        xs = sample(spec, dom, opts.points, seeds[-1])
        spans.append(range(len(lines), len(lines) + len(xs)))
        lines += [(dom.call or spec.func, x) for x in xs]
    with ThreadPoolExecutor() as threads:
        res = list(threads.map(fuzz_diff.evaluate, exes, [lines] * len(exes)))
    # a timeout of one revision only is probably the load of the machine, so it gets a retry
    timed_out = [y == "timeout" for y in res]
    if sum(timed_out) == 1 and not any(
        isinstance(y, str) and y != "timeout" for y in res
    ):
        i = timed_out.index(True)
        res[i] = fuzz_diff.evaluate(exes[i], lines)
    b, h = res
    if isinstance(b, str) or isinstance(h, str):
        status = [y if isinstance(y, str) else "ok" for y in (b, h)]
        accuracy.log(f"{spec.name}: skipped, base {status[0]}, head {status[1]}")
        note = f"skipped: base {status[0]}, head {status[1]}"
        # The reason takes the place of the counts.
        row = [spec.name, "-", note, *["-"] * (len(COLUMNS) - 3)]
        return [row], [], False, True

    rows: list[list[str]] = []
    flagged: list[tuple[float, str]] = []
    failed = inconclusive = False
    for desc, span, seed in zip(descs, spans, seeds, strict=True):
        # the differing inputs, the candidates of (a) and (b), and a sample of the others
        diff = [i for i in span if fuzz_diff.differ(b[i], h[i], 0.0)]
        cand_a = [
            i
            for i in diff
            if b[i] is not None and not nan(b[i]) and (h[i] is None or nan(h[i]))
        ]
        cand_b = [i for i in diff if nan(b[i]) and h[i] is not None and not nan(h[i])]
        rest = [i for i in diff if i not in {*cand_a, *cand_b}]
        rest = random.Random(seed).sample(rest, min(opts.max_refs, len(rest)))
        ids = [*cand_a, *cand_b, *rest]
        accuracy.log(
            f"{spec.name}: {len(cand_a) + len(cand_b)} inputs of (a) or (b) need references"
        )
        refs = references(spec, [lines[i][1] for i in ids], pool)
        errs = {
            i: (error(b[i], r), error(h[i], r))
            for i, r in zip(ids, refs, strict=True)
            if isinstance(r, list)
        }
        noref = {i for i in ids if i not in errs}

        # (c), from the sample, ignoring inputs whose errors are both small
        kinds = {
            i: fuzz_diff.verdict(*errs[i])
            for i in rest
            if i in errs and max(errs[i]) >= ACCURATE
        }
        better = sum(k == "improved" for k in kinds.values())
        worse_ids = [i for i, k in kinds.items() if k == "regressed"]
        worse = len(worse_ids)
        loss = net_loss(worse, better)
        # (inputs without a reference vote for neither side, so they could hide a loss)
        noref_rest = sum(i in noref for i in rest)
        hidden = not loss and net_loss(worse + noref_rest, better)
        sample_errs = [errs[i] for i in rest if i in errs]
        peak = " -> ".join(
            fmt(max((e[j] for e in sample_errs), default=None)) for j in (0, 1)
        )
        fail_a = [i for i in cand_a if i in errs and errs[i][0] <= NAN_BASE]
        fail_b = [i for i in cand_b if i in errs and errs[i][1] > NAN_HEAD]
        panics = [i for i in diff if h[i] is None and b[i] is not None]
        # (a NaN error comes from a NaN reference)
        unknown = any(
            i in noref or math.isnan(sum(errs[i])) for i in [*cand_a, *cand_b]
        )
        fails = "".join(
            c
            for c, bad in zip("abcd", (fail_a, fail_b, loss, panics), strict=True)
            if bad
        )
        failed |= bool(fails)
        inconclusive |= unknown or hidden
        rows.append(
            [
                spec.name,
                desc,
                f"{len(diff)} / {len(span)}",
                str(better),
                str(worse),
                peak,
                str(len(cand_a)),
                str(len(cand_b)),
                str(len(panics)),
                fails + "?" * unknown + f" ({noref_rest} no ref)" * hidden,
            ]
        )

        labels = {
            "worse": worse_ids,
            "a": fail_a,
            "b": fail_b,
            "panic": panics,
            "no ref": noref,
        }
        for i in sorted({i for v in labels.values() for i in v}):
            show = ["panic" if y is None else repr(y) for y in (b[i], h[i])]
            eb, eh = errs.get(i, (math.nan, math.inf))
            ref = json.dumps(refs[ids.index(i)][0]) if i in errs else "-"
            note = ", ".join(k for k, v in labels.items() if i in v)
            cols = [
                spec.name,
                desc,
                key(lines[i][1]),
                *show,
                ref,
                fmt(eb),
                fmt(eh),
                note,
            ]
            flagged.append((eh, "\t".join(cols)))
    return rows, flagged, failed, inconclusive


def positive(s: str) -> int:
    if (n := int(s)) < 1:
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
        "--affected",
        action="store_true",
        help="only the specs that depend on the files of src/ changed since --base"
        " (see affected.py); exit at once if there are none",
    )
    ap.add_argument(
        "--points",
        type=positive,
        default=2**14,
        help="inputs per domain (default: %(default)s)",
    )
    ap.add_argument(
        "--max-inputs",
        type=positive,
        help="lower --points to a power of 2 so that all domains have at most this many"
        " inputs (at least MIN_POINTS per domain)",
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
    if unknown := only - {n for s in specs.SPECS for n in (s.name, s.func)}:
        ap.error(f"unknown spec or function: {', '.join(sorted(unknown))}")
    chosen = [s for s in specs.SPECS if not only or only & {s.name, s.func}]
    if opts.affected:
        # (both sides of a rename, and with no --head also the untracked files)
        git = ["git", "-C", str(accuracy.ROOT)]
        revs = [opts.base, *([opts.head] if opts.head else [])]
        diff = [*git, "diff", "--name-only", "--no-renames", "-z", *revs, "--", "src"]
        out = subprocess.check_output(diff, text=True)
        if not opts.head:
            new = [*git, "ls-files", "-oz", "--exclude-standard", "--", "src"]
            out += subprocess.check_output(new, text=True)
        files = [accuracy.ROOT / f for f in out.split("\0") if f]
        names = set(affected.affected(files, opts.base, opts.head))
        chosen = [s for s in chosen if names & {s.name, s.func}]
        if not chosen:
            print("no function is affected")
            sys.exit(0)
    if opts.max_inputs:
        each = opts.max_inputs // sum(len(s.domains) for s in chosen)
        if each < MIN_POINTS:
            ap.error(f"--max-inputs is below {MIN_POINTS} inputs per domain")
        opts.points = min(opts.points, 1 << (each.bit_length() - 1))
    fuzz_diff.CACHE.mkdir(parents=True, exist_ok=True)
    exes = [fuzz_diff.program(opts.base), fuzz_diff.program(opts.head)]
    mp.mp.dps = specs.DPS

    print("|", " | ".join(COLUMNS), "|")
    print("|", " | ".join("---" for _ in COLUMNS), "|")
    flagged: list[tuple[float, str]] = []
    failed = inconclusive = False
    with multiprocessing.Pool(opts.jobs, maxtasksperchild=500) as pool:
        for spec in chosen:
            rows, f, bad, unsure = compare(spec, exes, pool, opts)
            for row in rows:
                print("|", " | ".join(row), "|", flush=True)
            flagged += f
            failed |= bad
            inconclusive |= unsure

    head = (
        "function\tdomain\tinput\tbase\thead\treference\tbase error\thead error\tnote"
    )
    lines = [line for _, line in sorted(flagged, reverse=True)]
    TSV.write_text("".join(f"{line}\n" for line in [head, *lines]))
    code = 1 if failed else 2 if inconclusive else 0
    print(
        f"\n{len(lines)} inputs flagged, in {TSV}; {('ok', 'FAIL', 'INCONCLUSIVE')[code]}"
    )
    sys.exit(code)


if __name__ == "__main__":
    main()
