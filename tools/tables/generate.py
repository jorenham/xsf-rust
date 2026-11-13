# /// script
# requires-python = ">=3.15"
# dependencies = ["mpmath==1.4.1", "gmpy2==2.3.2", "python-flint==0.9.0"]
# ///
"""Replace xsref's references and tolerances in tables of `tests/data/xsref` by our own.

Usage (from the repository root; BASELINE holds xsref's tables as `convert.py` writes them):

    uv run tools/tables/generate.py BASELINE [--only TABLE,...] [--no-cache] [--jobs N]

For each entry of `mapping.ENTRIES`, every logical output of every row gets a status:

- convention (2): the inputs match a `mapping.Rule`, which xsf-rust must follow, and which must
  pass xsref's test (so xsref may differ where its tolerance allows it).
- computed (1): the mpmath reference at both `DPS`, which must agree to `AGREE`, round to the
  same f64, and be confirmed by Arb where its Arb version applies (or by an identity, without
  one). xsf-rust must be NaN or an infinity exactly where the reference is. The tolerance is
  the error of xsf (C++) in the table metric (see `tolerance`).
- legacy (0): the reference timed out, so xsref's value and tolerance stay.

Anything else is an error, and then no table is written. The inputs, their order, and the row
flags stay as they are. The references are cached in `.cache/refs/` (see `stamp`), which does
not change the output, and a report per table goes to `.cache/reports/`.
"""

import argparse
import hashlib
import json
import math
import multiprocessing
import signal
import struct
import sys
from collections import Counter
from collections.abc import Callable
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Final, NamedTuple

HERE: Final = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "accuracy"))

import flint
import mapping
import mpmath as mp
import specs

import accuracy

ROOT: Final = HERE.parent.parent
TABLES: Final = ROOT / "tests/data/xsref"
CACHE: Final = HERE / ".cache"

MAGIC: Final = b"XSFT"
VERSION: Final = 1
HEADER: Final = struct.Struct("<4sIIHHHH")
STATUS: Final = 1  # header flag: a status byte per logical output
SCALES: Final = 2  # header flag: an error scale per logical output
LEGACY, COMPUTED, CONVENTION = range(3)

DPS: Final = (50, 100)  # the digits of the two evaluations of a reference
AGREE: Final = 1e-30  # their largest difference, relative to max(|f|, scale)
ARB_DPS: Final = (75, 150, 300)
TIMEOUT: Final = 60  # seconds per evaluation
IDENTITY_ROWS: Final = 32  # the inputs that the identity of a table is checked on
CHANGES: Final = 12  # the changed references of an output that a report lists
# xsf counts as wrong from this error (half of the digits lost)
XSF_WRONG: Final = 2.0**-26

EPS: Final = sys.float_info.epsilon
DBL_MIN: Final = sys.float_info.min
DBL_MAX: Final = sys.float_info.max
DENORM_MIN: Final = 5e-324
ULP_AT_MAX: Final = EPS * DBL_MAX / 2

type Args = tuple[float, ...]
type Ref = list[str | None]  # [value, scale or None]
# table, sig, output, part, inputs, digits
type Task = tuple[str, str, int, int, Args, int]


def log(*a: object) -> None:
    print(*a, file=sys.stderr, flush=True)


def parts(z: complex) -> tuple[float, ...]:
    return (z.real, z.imag) if isinstance(z, complex) else (z,)


def abs_error(y: float, f: float) -> float:
    """xsf's `extended_absolute_error`: an infinity is one ulp past the largest f64."""
    if y == f or (math.isnan(y) and math.isnan(f)):
        return 0.0
    if math.isnan(y) or math.isnan(f):
        return math.inf
    for a, b in ((y, f), (f, y)):
        if math.isinf(a):
            return abs((math.copysign(DBL_MAX, a) - b) + math.copysign(ULP_AT_MAX, a))
    return abs(y - f)


def magnitude(f: float) -> float:
    """What xsf's `extended_relative_error` is relative to, with the sign of `f`."""
    if f == 0:
        return math.copysign(DENORM_MIN, f)
    if math.isinf(f):
        return math.copysign(DBL_MAX, f)
    return 1.0 if math.isnan(f) else f


def table_error(y: complex, f: complex, scale: float = 0.0) -> float:
    """The table metric of `y` for the reference `f`, like in `testing/src/xsref.rs`."""
    d = math.hypot(*map(abs_error, parts(y), parts(f)))
    m = [magnitude(x) for x in parts(f)]
    e = d / math.hypot(*m)
    if math.isinf(math.hypot(*m)):  # the modulus of a complex reference overflows
        e = d / 2 / math.hypot(*(x / 2 for x in m))
    return min(e, d / scale) if scale else e


def same_class(y: complex, f: complex, zeros: bool = True) -> bool:
    """Whether each part of `y` is NaN, the same infinity, or finite, like that of `f`; with
    `zeros`, also the same zero or nonzero (like in `testing/src/xsref.rs`)."""

    def same(a: float, e: float) -> bool:
        if math.isnan(e):
            return math.isnan(a)
        if math.isinf(e) or zeros and e == 0:
            return a == e and math.copysign(1, a) == math.copysign(1, e)
        return math.isfinite(a) and not (zeros and a == 0)

    return all(map(same, parts(y), parts(f), strict=True))


def number(s: str) -> complex:
    """A value of `testing/metric_vectors.tsv`: an f64, or `re,im`."""
    return complex(*map(float, s.split(","))) if "," in s else float(s)


def check_metric() -> None:
    """Check `table_error` on the test vectors that the reader checks too."""
    for line in (ROOT / "testing/metric_vectors.tsv").read_text().splitlines():
        if line and not line.startswith("#"):
            y, f, scale, e = line.split("\t")
            assert table_error(number(y), number(f), float(scale)) == float(e), line


class Baseline(NamedTuple):
    """A table of xsref, as `convert.py` writes it."""

    n_in: int
    n_out: int
    n_tol: int
    rows: list[bytes]  # the f64s of each row: inputs, outputs, tolerances
    flags: bytes  # the row flags


def read_baseline(path: Path) -> Baseline:
    data = path.read_bytes()
    magic, version, n_rows, n_in, n_out, n_tol, flags = HEADER.unpack_from(data)
    assert (magic, version, flags) == (MAGIC, VERSION, 0), (
        f"{path} is not xsref's table"
    )
    width = 8 * (n_in + n_out + n_tol)
    body = data[HEADER.size :]
    rows = [body[r * width : (r + 1) * width] for r in range(n_rows)]
    return Baseline(n_in, n_out, n_tol, rows, body[n_rows * width :])


def f64s(data: bytes) -> tuple[float, ...]:
    return struct.unpack(f"<{len(data) // 8}d", data)


def value(x: tuple[float, ...]) -> complex:
    """The output with the f64s `x`."""
    return complex(*x) if len(x) == 2 else x[0]


def with_timeout[T](f: Callable[[], T]) -> T | None:
    """`f()`, or None if it takes longer than `TIMEOUT`."""
    signal.signal(signal.SIGALRM, accuracy._alarm)
    signal.alarm(TIMEOUT)
    try:
        return f()
    except accuracy.Timeout:
        return None
    finally:
        signal.alarm(0)


def digits(v: mp.mpf, dps: int) -> str:
    """`v` with `dps` digits, or exactly if it is a short binary fraction in the range of the
    f64s, such as the midpoint of two f64s, which `dps` digits could round to either side."""
    man, exp = v.man_exp
    if abs(man).bit_length() > 64 or not -1200 < exp < 1100:
        return mp.nstr(v, dps)
    return mp.nstr(v, max(dps, len(str(abs(man) * (5**-exp if exp < 0 else 2**exp)))))


def evaluate(task: Task) -> Ref | None:
    """A part of a reference as decimal strings with all digits, or None if it timed out."""
    table, sig, i, part, args, dps = task
    ref = mapping.ENTRIES[table, sig].outs[i].refs[part]

    def f() -> Ref:
        mp.mp.dps = dps
        v = ref(*map(mp.mpf, args))
        v, scale = v if isinstance(v, tuple) else (v, None)
        return [digits(mp.mpf(v), dps), None if scale is None else mp.nstr(scale, dps)]

    try:
        return with_timeout(f)
    except Exception as e:
        e.add_note(f"{table} {sig} [{i}] at {args}")
        raise


def residual(task: tuple[str, str, Args]) -> float | None:
    """The residual of the identity of a table, or None if it timed out."""
    table, sig, args = task

    def f() -> float:
        mp.mp.dps = DPS[1]
        identity = mapping.ENTRIES[table, sig].identity
        assert identity
        return float(abs(identity.residual(*map(mp.mpf, args))))

    return with_timeout(f)


def stamp() -> dict[str, str]:
    """What a cached reference depends on, other than its key."""
    return {
        "generate.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "specs.py": hashlib.sha256(Path(specs.__file__).read_bytes()).hexdigest(),
        "mapping.py": hashlib.sha256(Path(mapping.__file__).read_bytes()).hexdigest(),
        "mpmath": mp.__version__,
        "python-flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
        "dtype": "f64",
    }


def references(tasks: list[Task], use_cache: bool, pool: Pool) -> dict[Task, Ref]:
    """The references of the tasks of one table, without those that timed out."""
    table, sig = tasks[0][:2]
    path = CACHE / "refs" / f"{table}-{sig}.json"

    def key(task: Task) -> str:
        _, _, i, part, args, dps = task
        return f"{i}.{part}/{dps}/{struct.pack(f'<{len(args)}d', *args).hex()}"

    cached: dict[str, Ref] = {}
    if use_cache and path.exists():
        data = json.loads(path.read_text())
        cached = data["refs"] if data["stamp"] == stamp() else {}

    todo = [t for t in tasks if key(t) not in cached]
    if todo:
        log(f"{table} {sig}: computing {len(todo)} references")
    new = zip(todo, pool.imap(evaluate, todo, chunksize=4), strict=True)
    cached |= {key(t): ref for t, ref in new if ref is not None}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"stamp": stamp(), "refs": cached}))
    return {t: cached[key(t)] for t in tasks if key(t) in cached}


class Problem(Exception):
    """Something that needs a decision in the mapping."""


def require(ok: bool, why: str) -> None:
    if not ok:
        raise Problem(why)


def arb_confirms(arb: mapping.Arb, args: Args, ref: Ref) -> bool | None:
    """Whether the Arb ball of a reference is within `AGREE` of it; None without a ball."""
    v, scale = ref
    for dps in ARB_DPS:
        flint.ctx.dps = dps
        y = arb(*map(flint.arb, args))
        if y is None:
            return None
        x = flint.arb(v)
        diff = abs(y - x) - AGREE * abs(x).max(flint.arb(scale or "0"))
        if diff <= 0:
            return True
        if diff > 0:
            return False
    return False


def validated(
    lo: Ref, hi: Ref, arb: mapping.Arb | None, args: Args
) -> tuple[float, float, bool]:
    """The f64 value and scale (0 for none) of a part of a reference, and whether Arb
    confirms it. Its evaluations `lo` and `hi` must agree to `AGREE` and round to the same
    f64s."""
    mp.mp.dps = DPS[1]
    (v_lo, s_lo), (v_hi, s_hi) = ([mp.mpf(x or "0") for x in r] for r in (lo, hi))
    require(mp.isfinite(v_hi), f"{v_hi} is not finite, which needs a rule")
    bound = AGREE * max(abs(v_hi), s_hi)
    close = abs(v_lo - v_hi) <= bound and abs(s_lo - s_hi) <= bound
    # the nearest f64s of the decimal strings
    (f_lo, scale_lo), (f, scale) = ([float(x or "0") for x in r] for r in (lo, hi))
    require(close and (f_lo, scale_lo) == (f, scale), f"{lo} and {hi} differ")

    confirmed = None if arb is None else arb_confirms(arb, args, hi)
    require(confirmed is not False, f"Arb does not confirm {hi}")
    # a scale must be a normal f64, and a subnormal f has the absolute precision of DBL_MIN
    if not DBL_MIN <= scale <= DBL_MAX:
        scale = DBL_MIN if 0 < abs(f) < DBL_MIN else 0.0
    return f, scale, bool(confirmed)


class Computed(NamedTuple):
    """A computed output of a row."""

    ref: complex  # the serialized reference
    scale: float
    e_rust: float  # the errors of xsf-rust and xsf in the table metric
    e_xsf: float
    by_arb: bool  # confirmed by Arb
    xsref: complex  # what it replaces


class Cell(NamedTuple):
    """A logical output of a row, as it is written."""

    status: int
    ref: bytes  # its f64s
    tol: bytes
    scale: float = 0.0


class Result(NamedTuple):
    table: str
    sig: str
    data: bytes  # the table, unless there are problems
    report: list[str]
    problems: list[str]


def tolerance(
    c: Computed, override: mapping.Override | None, budget: float, floor: float
) -> float:
    """The tolerance of a computed output: the error of xsf, so that the test (which accepts
    4 max(eps, tol)) means "at most 4 times less accurate than xsf". Under an override it is
    at most `budget`, the largest error of xsf-rust there, and an error of xsf above
    `XSF_WRONG` needs one. It is at least `floor`, the floor of the mapping there."""
    tol = min(c.e_xsf, budget) if override else c.e_xsf
    errors = f"xsf-rust: {c.e_rust:.2e}, xsf: {c.e_xsf:.2e}"
    require(tol <= XSF_WRONG, f"xsf is wrong, which needs an override ({errors})")
    tol = max(tol, floor)
    require(c.e_rust <= 4 * max(EPS, tol), f"xsf-rust fails the test ({errors})")
    return tol


def change(old: complex, new: complex) -> str:
    """How a reference changed: "same", "rounding" (at most an ulp), "value", or "class"."""

    def kind(a: float, b: float) -> int:
        if math.isnan(a) and math.isnan(b) or a.hex() == b.hex():
            return 0
        if not same_class(a, b):
            return 3
        return 1 if abs(a - b) <= math.ulp(max(abs(a), abs(b))) else 2

    kinds = map(kind, parts(old), parts(new), strict=True)
    return ("same", "rounding", "value", "class")[max(kinds)]


def describe(cells: dict[Args, tuple[Computed, float]]) -> list[str]:
    """The report of the computed outputs of a logical output."""
    kinds = {args: change(c.xsref, c.ref) for args, (c, _) in cells.items()}
    changes = sorted(
        (
            (table_error(c.xsref, c.ref), args, c)
            for args, (c, _) in cells.items()
            if kinds[args] in ("value", "class")
        ),
        key=lambda x: x[0],
        reverse=True,
    )
    columns = zip(*((c.e_rust, c.e_xsf, tol) for c, tol in cells.values()), strict=True)
    return [
        "references: "
        + ", ".join(
            f"{n} {kind}" for kind, n in sorted(Counter(kinds.values()).items())
        ),
        f"confirmed by Arb: {sum(c.by_arb for c, _ in cells.values())}",
        "largest e_rust, e_xsf, tol: " + ", ".join(f"{max(x):.1e}" for x in columns),
        *(
            f"  {args}: {c.xsref!r} -> {c.ref!r} ({size:.1e}), xsf: {c.e_xsf:.1e}"
            for size, args, c in changes[:CHANGES]
        ),
    ]


def generate(
    table: str,
    sig: str,
    old: Baseline,
    cmds: tuple[list[str], list[str]],
    use_cache: bool,
    pool: Pool,
) -> Result:
    entry = mapping.ENTRIES[table, sig]
    # an f32 output ("f") needs an f32 evaluation path first
    widths = [{"d": 1, "cd": 2}[t] for t in sig.split("-")[1].split("_")]
    assert len(entry.outs) == len(widths) == old.n_tol and sum(widths) == old.n_out
    n = len(widths)
    starts = [old.n_in + sum(widths[:i]) for i in range(n)]
    rows = [f64s(row) for row in old.rows]
    inputs = [row[: old.n_in] for row in rows]
    cells = [(r, i) for r in range(len(rows)) for i in range(n)]

    # a rule sees the f64 inputs, so signed zeros and NaN
    rules = {
        (r, i): rule
        for r, i in cells
        if (rule := next((x for x in entry.outs[i].rules if x.when(inputs[r])), None))
    }
    # each part of each other output at both `DPS`
    tasks = {
        (r, i): [
            [(table, sig, i, part, inputs[r], dps) for dps in DPS]
            for part in range(widths[i])
        ]
        for r, i in cells
        if (r, i) not in rules
    }
    todo = [t for ts in tasks.values() for lo_hi in ts for t in lo_hi]
    refs = references(todo, use_cache, pool)
    lines = [(entry.outs[i].func, inputs[r]) for r, i in cells]
    rust, xsf = (
        dict(zip(cells, accuracy.run(cmd, lines), strict=True)) for cmd in cmds
    )

    new: dict[tuple[int, int], Cell] = {}
    computed: dict[tuple[int, int], Computed] = {}
    problems: list[str] = []
    notes = [Counter[str]() for _ in range(n)]
    for r, i in cells:
        out, a, tol = entry.outs[i], starts[i], old.n_in + old.n_out + i
        old_ref = value(rows[r][a : a + widths[i]])
        try:
            if rule := rules.get((r, i)):
                # a convention must pass xsref's test, unless an override marks xsf (and so
                # xsref) wrong there, and xsf-rust must follow it already
                over = next((o for o in out.overrides if o.when(inputs[r])), None)
                ok = table_error(rule.value, old_ref) <= 4 * max(EPS, rows[r][tol])
                require(
                    ok or over is not None,
                    f"xsref has {old_ref}, the rule {rule.value} ({rule.why})",
                )
                if over:
                    notes[i][f"override: {over.why}"] += 1
                ok = same_class(rust[r, i], rule.value)
                require(
                    ok, f"xsf-rust has {rust[r, i]}, the rule {rule.value} ({rule.why})"
                )
                notes[i][f"convention: {rule.why}"] += 1
                for name, x in (("xsref", old_ref), ("xsf", xsf[r, i])):
                    differs = not same_class(x, rule.value)
                    notes[i][f"conventions that {name} does not follow"] += differs
                ref = struct.pack(f"<{widths[i]}d", *parts(rule.value))
                new[r, i] = Cell(CONVENTION, ref, struct.pack("<d", 0.0))
            elif any(t not in refs for lo_hi in tasks[r, i] for t in lo_hi):
                notes[i]["legacy: the reference timed out"] += 1
                # the bytes of xsref's reference and tolerance
                ref = old.rows[r][8 * a :][: 8 * widths[i]]
                new[r, i] = Cell(LEGACY, ref, old.rows[r][8 * tol :][:8])
            else:
                arbs = out.arbs or (None,) * widths[i]
                checked = [
                    validated(refs[lo], refs[hi], arb, inputs[r])
                    for (lo, hi), arb in zip(tasks[r, i], arbs, strict=True)
                ]
                f = value(tuple(f for f, _, _ in checked))
                # no tolerance covers NaN or an infinity
                ok = same_class(rust[r, i], f, zeros=False)
                require(ok, f"xsf-rust has {rust[r, i]}, the reference {f}")
                # a complex output has no scale
                scale = checked[0][1] if widths[i] == 1 else 0.0
                e_rust, e_xsf = (table_error(y[r, i], f, scale) for y in (rust, xsf))
                if not same_class(xsf[r, i], f, zeros=False):
                    e_xsf = math.inf
                by_arb = all(confirmed for _, _, confirmed in checked)
                computed[r, i] = Computed(f, scale, e_rust, e_xsf, by_arb, old_ref)
        except Problem as e:
            problems.append(f"[{i}] at {inputs[r]}: {e}")

    # the largest error of xsf-rust under each override
    overrides = {
        (r, i): next((o for o in entry.outs[i].overrides if o.when(inputs[r])), None)
        for r, i in computed
    }
    floors = {
        (r, i): next((f for f in entry.outs[i].floors if f.when(inputs[r])), None)
        for r, i in computed
    }
    budgets: dict[tuple[int, mapping.Override | None], float] = {}
    for (r, i), c in computed.items():
        key = (i, overrides[r, i])
        budgets[key] = max(budgets.get(key, 0.0), c.e_rust)
    described: list[dict[Args, tuple[Computed, float]]] = [{} for _ in range(n)]
    for (r, i), c in computed.items():
        override, floor = overrides[r, i], floors[r, i]
        try:
            tol = tolerance(
                c, override, budgets[i, override], floor.tol if floor else 0
            )
        except Problem as e:
            problems.append(f"[{i}] at {inputs[r]}: {e}")
            continue
        ref = struct.pack(f"<{widths[i]}d", *parts(c.ref))
        new[r, i] = Cell(COMPUTED, ref, struct.pack("<d", tol), c.scale)
        if override:
            notes[i][f"override: {override.why}"] += tol != c.e_xsf
        if floor:
            notes[i][f"floor {floor.tol:g}: {floor.why}"] += tol == floor.tol
        described[i][inputs[r]] = (c, tol)

    # the identity, on an even spread of the computed inputs
    identity = entry.identity
    assert identity or all(out.arbs for out in entry.outs), "an identity is needed"
    picks = sorted(
        {inputs[r] for r, _ in computed if identity and identity.when(inputs[r])}
    )
    if len(picks) > IDENTITY_ROWS:
        step = (len(picks) - 1) / (IDENTITY_ROWS - 1)
        picks = [picks[round(j * step)] for j in range(IDENTITY_ROWS)]
    residuals = pool.map(residual, [(table, sig, args) for args in picks])
    for args, res in zip(picks, residuals, strict=True):
        if res is None or not res <= AGREE:
            problems.append(f"the identity fails at {args}: {res}")

    status = Counter(cell.status for cell in new.values())
    counts = [status[code] for code in (COMPUTED, CONVENTION, LEGACY)]
    largest = max(filter(None, residuals), default=0.0)
    report = [
        f"# {table} {sig}",
        "{} computed, {} convention, {} legacy".format(*counts),
        *(
            [f"identity: {len(picks)} inputs, residuals up to {largest:.1e}"]
            * bool(picks)
        ),
    ]
    for i, out in enumerate(entry.outs):
        report += ["", f"## [{i}] {out.func}"]
        report += [f"{k:6} {name}" for name, k in sorted(notes[i].items()) if k]
        report += describe(described[i])

    data = b""
    if not problems:
        scaled = any(cell.scale for cell in new.values())
        flags = STATUS | (SCALES if scaled else 0)
        chunks = [HEADER.pack(MAGIC, VERSION, len(rows), old.n_in, old.n_out, n, flags)]
        for r, row in enumerate(old.rows):
            outs = [new[r, i] for i in range(n)]
            chunks.append(row[: 8 * old.n_in])
            chunks += [c.ref for c in outs] + [c.tol for c in outs]
            if scaled:
                chunks.append(struct.pack(f"<{n}d", *(c.scale for c in outs)))
        data = b"".join([*chunks, old.flags, bytes(new[c].status for c in cells)])
    return Result(table, sig, data, report, problems)


def main() -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("baseline", type=Path, help="the tables that `convert.py` writes")
    ap.add_argument("--only", help="comma-separated tables")
    ap.add_argument("--no-cache", action="store_true", help="recompute all references")
    ap.add_argument("--jobs", type=int, default=2, help="reference worker processes")
    opts = ap.parse_args()
    baseline: Path = opts.baseline
    only = (
        opts.only.split(",") if opts.only else [table for table, _ in mapping.ENTRIES]
    )

    check_metric()
    lines = (TABLES / "MANIFEST.tsv").read_text().splitlines()
    manifest = [line.split("\t") for line in lines]
    index = {(line[0], line[1]): line for line in manifest}
    cmds = ([*accuracy.build_rust(), "xsf"], accuracy.build_cpp())

    results: list[Result] = []
    # recycle the workers, since mpmath caches grow at high precision
    with multiprocessing.Pool(opts.jobs, maxtasksperchild=500) as pool:
        for table, sig in mapping.ENTRIES:
            if table not in only:
                continue
            old = read_baseline(baseline / table / f"{sig}.bin")
            res = generate(table, sig, old, cmds, not opts.no_cache, pool)
            results.append(res)
            index[table, sig][-1] = hashlib.sha256(res.data).hexdigest()
            (CACHE / "reports").mkdir(parents=True, exist_ok=True)
            text = "\n".join([*res.report, "", *res.problems])
            (CACHE / "reports" / f"{table}-{sig}.txt").write_text(text + "\n")
            log(f"{table} {sig}: {res.report[1]}, {len(res.problems)} problems")

    problems = [f"{res.table} {res.sig} {p}" for res in results for p in res.problems]
    if problems:
        sys.exit("\n".join(["no table written:", *problems[:40]]))
    for res in results:
        (TABLES / res.table / f"{res.sig}.bin").write_bytes(res.data)
    text = "".join("\t".join(line) + "\n" for line in manifest)
    (TABLES / "MANIFEST.tsv").write_text(text)
    log(f"wrote {len(results)} tables; the reports are in {CACHE / 'reports'}")


if __name__ == "__main__":
    main()
