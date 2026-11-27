"""Find inputs where the working tree (or `--head`) is less accurate than a base revision.

Usage (from the repository root; or `just fuzz-diff`):

    uv run --project tools tools/accuracy/fuzz_diff.py [--base REV] [--head REV]
        [--only NAME ...] [--cases N] [--seed S] [--max-refs N] [--jobs N]

Steps:

1. Build the dump program of `accuracy.py` for both revisions; that of a commit is built from
   `git archive` and cached in `.cache/fuzz-diff/`.
2. Evaluate both on the same `--cases` seeded inputs per spec of `specs.py` (or of `--only`, by
   spec or function name): samples of its domains, with arguments moved to (near) integers,
   half-integers, domain bounds or special values, one set to (near) an integer plus a signed
   sum of others, and imaginary parts set to signed zeros or tiny values. Skip a spec if a
   revision lacks its function or times out on it.
3. Flag the inputs whose results differ by more than `RTOL` or in the sign of a zero, or where
   one revision panics. Repeat 2 and 3 for `NEAR_ROUNDS` rounds of `--cases` mutants of the
   flagged inputs (else of those whose results differ at all, else of new samples).
4. Compute the references of at most `--max-refs` flagged inputs per spec (those with a panic
   first), with Arb where `ARB` has the function, else as `accuracy.py` does, and compare the
   errors of both revisions (see `judge`).
5. Print the counts, write the regressions to a new `.cache/fuzz-diff-*.tsv`, and exit with 1 if
   there are any (including a timeout of the head alone).

`test_fuzz_diff.py` checks the references and verdicts of known hard cases.
"""

import argparse
import cmath
import hashlib
import math
import multiprocessing
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Final, TypeGuard

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
# the relative accuracy of an Arb reference
BITS: Final = 64
# the exponents of 2 of the numbers that round to a nonzero finite f64 (or to 0, from a midpoint)
F64_EXPS: Final = range(
    sys.float_info.min_exp - sys.float_info.mant_dig, sys.float_info.max_exp + 1
)
# replaces a zero argument, for the side of a branch cut or the sign of a zero result
TINY: Final = math.ulp(0.0)
# the offsets of a zero imaginary argument: each the square of the last, so that one is also
# small next to a tiny real argument
SIDES: Final = tuple(mp.mpf(TINY) ** 2**k for k in range(3))
# the reference where the function rejects the input (and is NaN)
REJECTED: Final = "rejected"
# the smallest part of a complex mpmath reference (relative to its absolute value) that is
# accurate, as for a real reference in `accuracy.ref_task`
MP_TOL: Final = mp.mpf(10) ** -40
# the prefix of a reference that fails with an exception
FAILED: Final = "failed: "
# the exit code of a panic
PANIC: Final = 101
# seconds per run of a dump program, and per line in it
RUN_TIMEOUT: Final = 10
LINE_TIMEOUT: Final = 0.02
# rounds of mutants of the flagged inputs
NEAR_ROUNDS: Final = 3
CACHE: Final = accuracy.CACHE / "fuzz-diff"
IMAG: Final = {"im", "xi", "yi", "zi", "pi"}
# (of a random sign)
SPECIAL: Final = (0.0, TINY, DBL_MIN, DBL_MAX, math.inf, math.nan)
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


def _elliprc(x: flint.arb, y: flint.arb) -> flint.arb | None:
    return flint.acb.elliptic_rc(x, y).real if x >= 0 else None


def _elliprj(
    x: flint.arb, y: flint.arb, z: flint.arb, p: flint.arb
) -> flint.arb | None:
    return flint.acb.elliptic_rj(x, y, z, p).real if min(x, y, z) >= 0 else None


def _hypu(a: flint.arb, b: flint.arb, x: flint.arb) -> flint.arb | None:
    u = x.hypgeom_u(a, b)
    return u if x >= 0 or u.is_finite() else None


def _lambertw_c(re: flint.arb, im: flint.arb) -> flint.acb:
    return flint.acb(re, im).lambertw()


# with the complex functions that `specs.ARB` lacks, since their specs use mpmath (whose
# `lambertw` ignores an imaginary part of e.g. 2^-1100, and so the side of its cut), and the
# real functions outside the domains of their specs: None where they reject the input (and are
# NaN, e.g. U where it is not real), and the Cauchy principal value (the real part, the mean of
# both sides of the cut) for Carlson's integrals with a negative last argument
ARB: Final[dict[str, Callable[..., flint.arb | flint.acb | None]]] = {
    **specs.ARB,
    "hyp1f1#c": _hyp1f1_c,
    "hyp2f1#c": _hyp2f1_c,
    "lambertw#c": _lambertw_c,
    "elliprc": _elliprc,
    "elliprj": _elliprj,
    "hypu": _hypu,
}


def build(root: Path, target: Path) -> Path:
    """Build the dump program of the source tree `root` into `target`; returns its path."""
    manifest = root / accuracy.HERE.relative_to(accuracy.ROOT) / "Cargo.toml"
    cmd = ["cargo", "build", "--release", "--quiet", "--manifest-path", str(manifest)]
    subprocess.run([*cmd, "--target-dir", str(target)], check=True)
    return target / "release" / "xsf-accuracy"


def _install(exe: Path, dst: Path) -> Path:
    """Copy `exe` to `dst`, unless it exists (renamed into place, for concurrent runs)."""
    if not dst.exists():
        part = dst.with_suffix(f".{os.getpid()}")
        shutil.copy2(exe, part)
        part.replace(dst)
    return dst


def program(rev: str | None) -> Path:
    """The dump program of the working tree (`rev` None), or of a commit."""
    if rev is None:
        # built in this checkout (not in `CACHE`, which worktrees can share), and run from a
        # copy named by its hash, which no later build replaces
        exe = build(accuracy.ROOT, accuracy.HERE / "target")
        with exe.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256").hexdigest()
        return _install(exe, CACHE / f"xsf-accuracy-{digest[:16]}")
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
            _install(build(Path(tmp), Path(tmp) / "target"), exe)
    return exe


def run(exe: Path, lines: list[Line]) -> list[Value]:
    """The results of the lines (None where it panics), within the time of one run."""
    deadline = time.monotonic() + RUN_TIMEOUT + LINE_TIMEOUT * len(lines)
    inp = [accuracy.dump_line(*line) for line in lines]
    out: list[Value] = []
    # after the program dies without unwinding, the next run is of one line, to find that line
    one = False
    while len(out) < len(lines):
        todo = inp[len(out) : len(out) + (1 if one else len(lines))]
        one = False
        try:
            res = subprocess.run(
                [str(exe), "xsf"],
                input="".join(todo),
                capture_output=True,
                text=True,
                check=True,
                timeout=deadline - time.monotonic(),
            )
            done = accuracy.parse(res.stdout)
            assert len(done) == len(todo), (exe, len(done), len(todo))
            out += done
        except subprocess.CalledProcessError as e:
            if "unknown function" in e.stderr:
                raise
            # (drops a partial last line)
            done = accuracy.parse(e.stdout[: e.stdout.rfind("\n") + 1])
            # a panic writes the results before it as it unwinds; else (e.g. an abort or a
            # signal) the unwritten results are lost
            panic = len(todo) == 1 or (e.returncode == PANIC and len(done) < len(todo))
            out += done[: len(todo) - 1]
            if panic:
                out.append(None)
            one = not panic
    return out


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


def _parts(z: complex) -> tuple[float, float]:
    z = complex(z)
    return z.real, z.imag


def differ(x: Value, y: Value, rtol: float) -> bool:
    """Whether a part of two results differs by more than `rtol` (relative to the larger) or in
    the sign of a zero, or only one is a panic."""
    if x is None or y is None:
        return (x is None) != (y is None)
    return not all(
        (
            math.isclose(u, v, rel_tol=rtol)
            and math.copysign(1, u) == math.copysign(1, v)
        )
        or (math.isnan(u) and math.isnan(v))
        for u, v in zip(_parts(x), _parts(y))
    )


def _arb_str(x: flint.arb) -> str:
    """A part of an Arb reference: exact if the ball is, within `F64_EXPS` (so that it rounds
    correctly to f64, e.g. at a midpoint), "0" (of an unknown sign) if it contains 0, else 30
    digits."""
    if x.is_exact():
        man, exp = (int(v) for v in x.man_exp())
        if exp + man.bit_length() in F64_EXPS:
            return str(man << exp) if exp >= 0 else f"{man * 5**-exp}e{exp}"
    return "0" if x.contains(0) else x.str(30, radius=False)


def _value(name: str, args: accuracy.Args) -> accuracy.Ref | str:
    """The reference at `args`, with Arb where `ARB` has the function, else as
    `accuracy.ref_task`, or `REJECTED` (where Arb rejects it at the highest precision), or an
    error string (e.g. for a NaN of mpmath where the function has a limit); the parts of a
    complex reference are 0 where their sign is unknown."""
    spec = specs.BY_NAME[name]
    if (f := ARB.get(name)) is not None:
        best: flint.acb | None = None
        for prec in ARB_PRECS:
            with flint.ctx.workprec(prec):
                v = f(*map(flint.arb, args))
            if v is None:
                continue
            y = flint.acb(v)
            if y.rel_accuracy_bits() >= BITS:
                best = y
                if all(p.is_exact() or not p.contains(0) for p in (y.real, y.imag)):
                    break
        if best is not None:
            re, im = map(_arb_str, (best.real, best.imag))
            if spec.complex:
                return [[re, im], None]
            return [re, None] if best.imag.contains(0) else f"complex result {best}"
        if v is None:
            return REJECTED
    r = accuracy.ref_task((name, args))
    if isinstance(r, list) and any(mp.isnan(p) for p in _mpfs(r)):
        return "NaN reference"
    if spec.complex and isinstance(r, list) and isinstance(r[0], list):
        tol = abs(mp.mpc(*r[0])) * MP_TOL
        r[0] = ["0" if abs(mp.mpf(p)) < tol else p for p in r[0]]
    return r


def _sided(args: accuracy.Args, at: set[int], tiny: float) -> accuracy.Args:
    """The arguments, with the zeros at `at` replaced by `tiny` times their sign; `tiny` can be
    an mpf below the f64 range (as in `SIDES`), which `_value` takes exactly, as for a float."""
    return tuple(
        math.copysign(1, x) * tiny if i in at else x for i, x in enumerate(args)
    )


def _values(r: accuracy.Ref) -> list[str]:
    """The value of a reference, as a list of its parts."""
    v = r[0]
    return [v] if isinstance(v, str) else v or []


def _mpfs(r: accuracy.Ref) -> list[mp.mpf]:
    return [mp.mpf(p) for p in _values(r)]


def ref(line: Line) -> accuracy.Ref | str:
    """The reference (see `_value`) or an error string; with a signed zero imaginary argument,
    that on its side of a branch cut (see `_side`); and a part that is 0 has the sign of the
    limit from the side of the real zero arguments, as `TINY`, where that changes with them."""
    try:
        return _ref(line)
    except Exception as e:  # noqa: BLE001
        return f"{FAILED}{type(e).__name__}: {e}"


def _finite(r: accuracy.Ref | str) -> TypeGuard[accuracy.Ref]:
    return isinstance(r, list) and all(map(mp.isfinite, _mpfs(r)))


def _side(
    name: str, args: accuracy.Args, imag: set[int], r: accuracy.Ref
) -> accuracy.Ref:
    """The limit from the side of the signed zeros at `imag`: the value at the first offset of
    `SIDES` that agrees with that at the previous one (to `BITS` bits of the larger of its
    absolute value and `DBL_MIN`), else `r`, that at the point (as near a branch point of a
    small power, e.g. (1 - z)^0.01). If the reference fails, or is the same non-real value on
    both sides (which f(conj(z)) = conj(f(z)) rules out for most functions), that at the
    previous offset (or `r`)."""
    last = r
    for i, d in enumerate(SIDES):
        s = _value(name, _sided(args, imag, d))
        if not _finite(s):
            return last
        if i:
            u, v = (mp.mpc(*_mpfs(x)) for x in (last, s))
            if abs(u - v) <= mp.ldexp(max(abs(u), DBL_MIN), -BITS):
                return s
            if v.imag and s == _value(name, _sided(args, imag, -d)):
                return last
        last = s
    return r


def _ref(line: Line) -> accuracy.Ref | str:
    name, args = line
    spec = specs.BY_NAME[name]
    zeros = {i for i, x in enumerate(args) if x == 0}
    imag = {i for i in zeros if spec.args[i] in IMAG}
    real = zeros - imag
    r = _value(name, args)
    if imag and _finite(r):
        r = _side(name, args, imag, r)
    if not real or not isinstance(r, list) or all(_mpfs(r)):
        return r
    near = _sided(args, imag, TINY)
    s, t = (_value(name, _sided(near, real, d)) for d in (TINY, -TINY))
    if not isinstance(s, list) or not isinstance(t, list):
        return r
    parts = [
        repr(math.copysign(TINY, sign))
        if x == 0 and (sign := mp.sign(q)) == -mp.sign(u) != 0
        else p
        for p, x, q, u in zip(_values(r), _mpfs(r), _mpfs(s), _mpfs(t), strict=True)
    ]
    return [parts if spec.complex else parts[0], r[1]]


def verdict(eb: float, eh: float) -> str:
    if eh > max(WORSE * eb, FLOOR):
        return "regressed"
    if eb > max(WORSE * eh, FLOOR):
        return "improved"
    return "both wrong" if min(eb, eh) > FLOOR else "both right"


def _zeros(z: complex, r: accuracy.Ref) -> bool:
    """Whether the parts of `z` are 0 (below `DBL_MIN`, as in an underflow) where those of a
    finite reference are, and then of the same sign; a reference part that is exactly 0 has an
    unknown sign."""
    ps = _mpfs(r)
    return any(abs(p) > DBL_MAX for p in ps) or all(
        p == 0
        or (abs(u) >= DBL_MIN and abs(p) >= DBL_MIN)
        or (abs(u) < DBL_MIN and abs(p) < DBL_MIN and math.copysign(1, u) == mp.sign(p))
        for u, p in zip(_parts(z), ps)
    )


def judge(x: Value, y: Value, r: accuracy.Ref | str | None) -> tuple[str, float, float]:
    """The kind (see `KINDS`) of the results of the base and the head, with their errors; a
    panic counts as a NaN, which is right only for `REJECTED`. If both are right, but only one
    has the zeros of the reference (see `_zeros`), the other is wrong."""
    if r == REJECTED:
        eb, eh = (0.0 if v is None or cmath.isnan(v) else math.inf for v in (x, y))
        return verdict(eb, eh), eb, eh
    if not isinstance(r, list):
        return "unchecked" if r is None else "no ref", math.nan, math.nan
    eb, eh = (accuracy.error(math.nan if v is None else v, r) for v in (x, y))
    kind = verdict(eb, eh)
    if (
        kind == "both right"
        and x is not None
        and y is not None
        and (zx := _zeros(x, r)) != _zeros(y, r)
    ):
        kind = "regressed" if zx else "improved"
    return kind, eb, eh


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
    # those with a panic first: it is a regression only where the function has a value (not
    # where it rejects the input)
    todo = sorted(flagged, key=lambda i: (None not in (base[i], head[i]), rng.random()))
    todo = todo[: opts.max_refs]
    tasks = [(spec.name, lines[i][1]) for i in todo]
    refs = dict(zip(todo, pool.map(ref, tasks, chunksize=1), strict=True))
    # (a bug, which would otherwise only show as "no ref")
    failed = (r for r in refs.values() if isinstance(r, str) and r.startswith(FAILED))
    if (e := next(failed, None)) is not None:
        accuracy.log(f"{spec.name}: reference {e}")
    count = Counter[str]()
    worse: list[tuple[float, str]] = []
    for i in flagged:
        x, y, r = base[i], head[i], refs.get(i)
        kind, eb, eh = judge(x, y, r)
        count[kind] += 1
        if kind == "regressed":
            show = ["panic" if u is None else repr(u) for u in (x, y)]
            cols = [
                spec.name,
                " ".join(map(repr, lines[i][1])),
                *show,
                " ".join(_values(r)) if isinstance(r, list) else str(r),
            ]
            worse.append((eh, "\t".join([*cols, f"{eb:.1e}", f"{eh:.1e}"])))
    return [len(lines), len(flagged), *(count[k] for k in KINDS)], worse


def report(regressions: list[str]) -> Path | None:
    """Print the regressions, and write them to a new TSV file (for concurrent runs) in
    `accuracy.CACHE`, if there are any; returns its path."""
    if not regressions:
        print("\nno regressions")
        return None
    head = "function\tinput\tbase\thead\treference\tbase error\thead error"
    with tempfile.NamedTemporaryFile(
        "w", prefix="fuzz-diff-", suffix=".tsv", dir=accuracy.CACHE, delete=False
    ) as f:
        f.writelines(f"{line}\n" for line in [head, *regressions])
    print(f"\n{len(regressions)} regressions, in {f.name}")
    for line in [head, *regressions[:10]]:
        print(line.replace("\t", "  "))
    return Path(f.name)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None
    )
    ap.add_argument("--base", default="origin/master", help="(default: %(default)s)")
    ap.add_argument("--head", help="the revision to check (default: the working tree)")
    ap.add_argument("--only", nargs="+", default=[], help="spec or function names")
    ap.add_argument(
        "--cases", type=int, default=20000, help="inputs per spec and round"
    )
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

    sys.exit(1 if report(regressions) else 0)


if __name__ == "__main__":
    main()
