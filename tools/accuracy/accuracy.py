"""Generate docs/accuracy.md: compare xsf-rust, xsf (C++), and Rust std against mpmath.

Usage (from the repository root; or `just accuracy`):

    uv run --project tools tools/accuracy/accuracy.py [--only NAME,... | --section TITLE]
        [--no-cache] [--jobs N] [--n N] [--arb-check]

Steps:

1. Generate the inputs of every domain in `specs.py` (seeded, so reproducible).
2. Evaluate them with xsf-rust and Rust std (`cargo run` of the dump program in this
   directory), and with xsf (C++), using `harness.cpp` compiled with `g++ -O2` against
   the xsf commit `XSF_COMMIT`, which is downloaded into `.cache/`; where xsf lacks the
   function, with SciPy or NumPy instead (see `scipy_funcs.py`).
3. Compute the mpmath references in parallel (with `--jobs` worker processes, 8 by default);
   these are cached in `.cache/refs/`, keyed by the inputs, the source of the reference
   function (normalized, see `_source`), and the mpmath version.
4. Write the table to `docs/accuracy.md` (or to stdout with `--only`, `--section`, or `--n`), and
   the inputs with the largest errors to `.cache/worst.txt`.
5. With `--arb-check`, compare the references with Arb (see `specs.ARB`), and write those that
   Arb refutes or cannot confirm to `.cache/arb_check.txt`.
"""

import argparse
import ast
import functools
import hashlib
import inspect
import io
import itertools
import json
import math
import multiprocessing
import platform
import random
import signal
import struct
import subprocess
import sys
import tarfile
import time
import tomllib
import urllib.request
import zlib
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from operator import itemgetter
from pathlib import Path
from types import CodeType, FrameType
from typing import Final, NamedTuple, Never

import flint
import mpmath as mp
import scipy_funcs
import specs

XSF_COMMIT: Final = "7fb9276f46aaf2752c308cdbb77b089120933c5f"
XSF_URL: Final = f"https://codeload.github.com/scipy/xsf/tar.gz/{XSF_COMMIT}"

HERE: Final = Path(__file__).resolve().parent
ROOT: Final = HERE.parent.parent
CACHE: Final = HERE / ".cache"

DBL_MIN: Final = 2.2250738585072014e-308
DBL_MAX: Final = sys.float_info.max
# how far (relative to |f|, or to the scale) a reference may be from the Arb ball: well below
# the f64 errors in the table, but above the rounding to the 30 digits that the cache keeps
ARB_TOL: Final = 1e-20

type Args = tuple[float, ...]
# [value, scale or None], where the value is [real, imaginary] for a complex function
type Ref = list[str | list[str] | None]
type Stats = tuple[float, float | None, int]


class Job(NamedTuple):
    spec: specs.Spec
    di: int  # domain index
    args: list[Args]


@dataclass(frozen=True)
class Result:
    spec: specs.Spec
    domain: str
    n: int  # number of valid references
    stats: dict[str, Stats | None]  # by implementation
    same: float | None  # share of bit-identical xsf-rust and xsf (or SciPy) results

    @property
    def notes(self) -> list[str]:
        """Its footnotes: those of the spec, and SciPy's if it has SciPy instead of xsf."""
        scipy = [scipy_funcs.FUNCS[self.spec.name][1]] if "scipy" in self.stats else []
        return self.spec.note.split() + scipy

    @property
    def second(self) -> str | None:
        """The name of what is in the second column, if anything."""
        if "cpp" in self.stats:
            return IMPLS["cpp"]
        if "scipy" in self.stats:
            numpy = scipy_funcs.FUNCS[self.spec.name][1] == scipy_funcs.NUMPY
            return "NumPy" if numpy else "SciPy"
        return None


def log(*a: object) -> None:
    print(*a, file=sys.stderr, flush=True)


def fetch_xsf() -> Path:
    inc = CACHE / f"xsf-{XSF_COMMIT}" / "include"
    if not inc.is_dir():
        log(f"downloading {XSF_URL}")
        with urllib.request.urlopen(XSF_URL) as r:
            data = r.read()
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
            tf.extractall(CACHE, filter="data")
    return inc


def build_cpp() -> list[str]:
    inc = fetch_xsf()
    exe = CACHE / f"harness-{XSF_COMMIT[:12]}"
    src = HERE / "harness.cpp"
    log(f"compiling {src.name}")
    cmd = ["g++", "-O2", "-std=c++17", "-I", str(inc), str(src), "-o", str(exe)]
    subprocess.run(cmd, check=True)
    return [str(exe)]


def build_rust() -> list[str]:
    opts = ["--release", "--quiet", "--manifest-path", str(HERE / "Cargo.toml")]
    subprocess.run(["cargo", "build", *opts], check=True)
    return ["cargo", "run", *opts, "--"]


def run(cmd: list[str], lines: list[tuple[str, Args]]) -> list[complex]:
    """Feed `<name> <args>...` lines to a dump program; return the f64 (or complex) results."""
    inp = "".join(f"{name} {' '.join(map(repr, args))}\n" for name, args in lines)
    out = subprocess.run(
        cmd, input=inp, capture_output=True, text=True, check=True
    ).stdout
    res: list[complex] = []
    for line in out.splitlines():
        re_im = [float(s) for s in line.split()]
        res.append(re_im[0] if len(re_im) == 1 else complex(*re_im))
    assert len(res) == len(lines), (cmd, len(res), len(lines))
    return res


class Timeout(Exception):
    pass


def _alarm(signum: int, frame: FrameType | None) -> Never:
    raise Timeout


def ref_task(task: tuple[str, Args]) -> Ref | str:
    """Evaluate one reference; returns [value, scale or None], or an error string. The value
    is a string, or a [real, imaginary] pair of strings for a complex function."""
    name, args = task
    spec = specs.BY_NAME[name]
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(spec.timeout)
    try:
        mp.mp.dps = specs.DPS
        v = spec.ref(*map(mp.mpf, args))
        scale = None
        if isinstance(v, tuple):
            v, scale = v
        sc = None if scale is None else mp.nstr(scale, 10)
        if spec.complex:
            v = mp.mpc(v)
            return [[mp.nstr(v.real, 30), mp.nstr(v.imag, 30)], sc]
        if isinstance(v, mp.mpc):
            if abs(v.imag) > abs(v.real) * mp.mpf(10) ** -40:
                return f"complex result {v}"
            v = v.real
        return [mp.nstr(v, 30), sc]
    except Timeout:
        return "timeout"
    except Exception as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"
    finally:
        signal.alarm(0)


def timed_ref_task(task: tuple[str, Args]) -> tuple[Ref | str, float]:
    t = time.perf_counter()
    return ref_task(task), time.perf_counter() - t


@functools.cache
def _defs(path: str) -> dict[int, list[ast.FunctionDef | ast.Lambda]]:
    """The functions and lambdas defined in a file, by their first line."""
    defs: defaultdict[int, list[ast.FunctionDef | ast.Lambda]] = defaultdict(list)
    for node in ast.walk(ast.parse(Path(path).read_text())):
        if isinstance(node, ast.Lambda):
            defs[node.lineno].append(node)
        elif isinstance(node, ast.FunctionDef):
            # the code of a decorated function starts at its first decorator
            defs[min(d.lineno for d in [node, *node.decorator_list])].append(node)
    return defs


def _source(code: CodeType) -> str:
    """The source of a function or lambda, normalized by `ast.unparse`, so that it does not
    depend on the formatting."""
    defs = _defs(code.co_filename)[code.co_firstlineno]
    nodes = [d for d in defs if getattr(d, "name", "<lambda>") == code.co_name]
    if len(nodes) > 1:
        # lambdas on the same line: the one whose body is where `code` returns from
        pos = set(code.co_positions())
        nodes = [
            d
            for d in nodes
            if isinstance(d, ast.Lambda)
            and (
                d.body.lineno,
                d.body.end_lineno,
                d.body.col_offset,
                d.body.end_col_offset,
            )
            in pos
        ]
    (node,) = nodes
    return ast.unparse(node)


def _sources(code: CodeType, seen: dict[str, str]) -> dict[str, str]:
    """Sources (see `_source`) of the `specs` functions that `code` uses, recursively."""
    for name in code.co_names:
        f = getattr(specs, name, None)
        if inspect.isfunction(f) and name not in seen:
            seen[name] = _source(f.__code__)
            _sources(f.__code__, seen)
    for c in code.co_consts:
        if inspect.iscode(c):
            _sources(c, seen)
    return seen


def ref_key(spec: specs.Spec, args: list[Args]) -> str:
    """Cache key: the inputs, and the sources of the reference and of what it calls. For a
    lambda or closure, which can share its source with other references, the key also has its
    default arguments and captured variables (by name, for functions)."""
    ref = spec.ref
    src = _sources(ref.__code__, {"": _source(ref.__code__)})
    key: list[object] = [args, sorted(src.items()), specs.DPS, mp.__version__]
    bound = [
        *(ref.__defaults__ or ()),
        *(c.cell_contents for c in ref.__closure__ or ()),
    ]
    if bound:
        for f in bound:
            if inspect.isfunction(f) and getattr(specs, f.__name__, None) is f:
                _sources(f.__code__, src)
                src[f.__name__] = _source(f.__code__)
        key[1] = sorted(src.items())
        key.append([getattr(f, "__name__", repr(f)) for f in bound])
    return hashlib.sha256(json.dumps(key).encode()).hexdigest()


def compute_refs(jobs: list[Job], use_cache: bool, nproc: int) -> list[list[Ref | str]]:
    """The references of each job: cached, or computed with `nproc` worker processes."""
    (CACHE / "refs").mkdir(parents=True, exist_ok=True)
    keys = [ref_key(job.spec, job.args) for job in jobs]
    paths = [
        CACHE / "refs" / f"{job.spec.name}-{job.di}-{key[:16]}.json"
        for job, key in zip(jobs, keys, strict=True)
    ]
    out: dict[int, list[Ref | str]] = {}
    if use_cache:
        for j, path in enumerate(paths):
            if (
                path.exists()
                and (data := json.loads(path.read_text()))["key"] == keys[j]
            ):
                out[j] = data["refs"]
    todo = [j for j in range(len(jobs)) if j not in out]
    tasks = [(jobs[j].spec.name, a) for j in todo for a in jobs[j].args]
    if tasks:
        log(f"computing {len(tasks)} references ({len(todo)} domains)")
        t0 = time.time()
        res: list[Ref | str] = []
        cost: defaultdict[str, float] = defaultdict(float)
        # recycle the workers, since mpmath caches grow at high precision
        with multiprocessing.Pool(nproc, maxtasksperchild=500) as pool:
            timed = pool.imap(timed_ref_task, tasks, chunksize=4)
            for i, ((name, _), (r, dt)) in enumerate(zip(tasks, timed, strict=True), 1):
                res.append(r)
                cost[name] += dt
                if i % 20000 == 0:
                    log(f"  {i}/{len(tasks)} ({time.time() - t0:.0f} s)")
        top = sorted(cost.items(), key=itemgetter(1), reverse=True)[:10]
        log("most expensive: " + ", ".join(f"{name} {dt:.0f} s" for name, dt in top))
        it = iter(res)
        for j in todo:
            out[j] = list(itertools.islice(it, len(jobs[j].args)))
            paths[j].write_text(json.dumps({"key": keys[j], "refs": out[j]}))
    return [out[j] for j in range(len(jobs))]


def arb_task(task: tuple[str, Args, Ref]) -> str:
    """Compare a reference with the Arb ball at 2, 4, and 8 times its digits; returns "", or a
    line that says whether Arb refutes it or cannot confirm it."""
    name, args, (value, scale) = task
    assert isinstance(value, str) and not isinstance(scale, list)
    for k in (2, 4, 8):
        flint.ctx.dps = k * specs.DPS
        ref = flint.arb(value)
        tol = abs(ref) if scale is None else abs(ref).max(flint.arb(scale))
        y = specs.ARB[name](*map(flint.arb, args))
        diff = abs(y - ref) - ARB_TOL * tol
        if diff <= 0:
            return ""
        if diff > 0:
            return f"refuted {name}{args}: mpmath {value}, Arb {y.str(30)}"
    return f"unconfirmed {name}{args}: mpmath {value}, Arb {y.str(30)}"


def arb_check(jobs: list[Job], refs: list[list[Ref | str]], nproc: int) -> None:
    """Check the references that have an Arb version (`specs.ARB`)."""
    lines: list[str] = []
    with multiprocessing.Pool(nproc) as pool:
        for job, rr in zip(jobs, refs, strict=True):
            if job.spec.name not in specs.ARB:
                continue
            tasks = [
                (job.spec.name, a, r)
                for a, r in zip(job.args, rr, strict=True)
                if isinstance(r, list)
            ]
            res = pool.map(arb_task, tasks, chunksize=4)
            n = Counter(line.partition(" ")[0] for line in res)
            log(
                f"{job.spec.name} [{job.spec.domains[job.di].desc(job.spec.args)}]: Arb"
                f" refutes {n['refuted']} and cannot confirm {n['unconfirmed']} of"
                f" {len(res)} references"
            )
            lines += filter(None, res)
    (CACHE / "arb_check.txt").write_text("".join(f"{line}\n" for line in lines))


def error(y: complex, ref: Ref) -> float:
    """Relative error of `y` w.r.t. the reference (see the header of docs/accuracy.md)."""
    value, scale = ref
    assert not isinstance(scale, list)
    if isinstance(value, list):
        y, r = complex(y), mp.mpc(*value)
        if math.isnan(y.real) or math.isnan(y.imag):
            return math.inf
        if abs(r) > DBL_MAX:
            return 0.0 if math.isinf(abs(y)) else math.inf
        if math.isinf(abs(y)):
            return math.inf
    else:
        assert isinstance(y, float)
        r = mp.mpf(value)
        if math.isnan(y):
            return math.inf
        if mp.isinf(r) or abs(r) > DBL_MAX:
            return 0.0 if math.isinf(y) and (y > 0) == (r > 0) else math.inf
    if abs(r) < DBL_MIN:  # underflow
        return 0.0 if abs(y) < DBL_MIN else float(abs(y - r) / DBL_MIN)
    d = abs(y - r)
    e = d / abs(r)
    if scale is not None:
        e = min(e, d / mp.mpf(scale))
    return float(e)


# the unit roundoff, the smallest error the geometric mean counts, so that a few exact results do
# not make it 0 (it is only 0 if all are exact)
GMEAN_FLOOR: Final = 2.0**-53


def stats(errs: list[float]) -> Stats | None:
    """(peak, geometric mean, number of infinite errors); the mean is over the finite errors."""
    if not errs:
        return None
    fin = [e for e in errs if math.isfinite(e)]
    peak = max(fin, default=0.0)
    gmean = (
        math.exp(math.fsum(math.log(max(e, GMEAN_FLOOR)) for e in fin) / len(fin))
        if peak
        else 0.0
    )
    nbad = len(errs) - len(fin)
    return (math.inf if nbad else peak), (gmean if fin else None), nbad


def bits(x: complex) -> tuple[bytes | str, ...]:
    """The f64 bits (of both parts, if complex), except that all NaNs are equal (the dumps do
    not keep NaN payloads)."""
    parts = (x.real, x.imag) if isinstance(x, complex) else (x,)
    return tuple("nan" if math.isnan(p) else struct.pack("<d", p) for p in parts)


IMPLS: Final = {"rust": "xsf-rust", "cpp": "xsf", "std": "std"}
NA: Final = "\u2014"  # for functions that xsf does not have
# footnote marks: *, dagger, double dagger, section sign, double bar, pilcrow, **, two daggers
FOOTNOTES: Final = (
    "*",
    "\u2020",
    "\u2021",
    "\u00a7",
    "\u2016",
    "\u00b6",
    "**",
    "\u2020\u2020",
)
# the smallest error that is marked with `BAD_MARK`
BAD: Final = 1e-6
BAD_MARK: Final = "\u2757"  # heavy exclamation mark, a red emoji


def fmt(x: float | str | None, ninf: int = 0) -> str:
    """Error in a fixed-width format, e.g. 1.1e-15, in monospace so that the digits align."""
    if x is None:
        return ""
    if isinstance(x, str):  # NA
        return x
    if x == 0:
        s = "0"
    elif not math.isfinite(x):
        # the number of infinite errors, after a no-break space to stay on one line
        s = f"\u221e\u00a0({ninf})"
    else:
        s = f"{x:.1e}"
    return f"<samp>{s}</samp>"


def _err_key(cell: str) -> tuple[float, int] | None:
    """Sort key of a formatted error, or None for "\u2014" and empty cells."""
    if not cell.startswith("<samp>"):
        return None
    s = cell.removeprefix("<samp>").removesuffix("</samp>")
    if s.startswith("\u221e"):
        return (math.inf, int(s.partition("(")[2].rstrip(")")))
    return (float(s), 0)


def mark_bad(cell: str) -> str:
    """Mark an error of at least `BAD`, or infinite errors, with a red `BAD_MARK` (GitHub drops
    CSS colors, and its math can be neither bold nor monospace)."""
    k = _err_key(cell.removeprefix("<b>").removesuffix("</b>"))
    return cell if k is None or k[0] < BAD else BAD_MARK + cell


def mark_best(errs: Iterable[str]) -> list[str]:
    """Make the smallest peaks and the smallest gmean bold, unless all are equal at the displayed
    precision; `errs` alternates peak and gmean per implementation."""
    errs = list(errs)
    for j in range(2):
        keys = {
            c: k for c in range(j, len(errs), 2) if (k := _err_key(errs[c])) is not None
        }
        best = min(keys.values(), default=None)
        if any(k != best for k in keys.values()):
            for c, k in keys.items():
                if k == best:
                    errs[c] = f"<b>{errs[c]}</b>"
    return errs


def table(names: list[str], rows: list[list[str]]) -> str:
    """HTML table with a two-level head, with the implementation `names`; markdown tables
    cannot span columns."""
    head = [
        '<th rowspan="2">Function</th>',
        '<th rowspan="2">Domain</th>',
        '<th rowspan="2">N</th>',
    ]
    head += [f'<th colspan="2">{name}</th>' for name in names]
    head.append('<th rowspan="2">same</th>')
    lines = ["<table>", "<thead>", f"<tr>{''.join(head)}</tr>"]
    lines.append("<tr>" + "<th>peak</th><th>gmean</th>" * len(names) + "</tr>")
    lines += ["</thead>", "<tbody>"]
    for name, domain, n, *errs, same in rows:
        cells = [f"<td>{name}</td>", f"<td>{domain}</td>"]
        cells += [
            f'<td align="right" nowrap>{x}</td>'
            for x in [n, *map(mark_bad, mark_best(errs)), same]
        ]
        lines.append(f"<tr>{''.join(cells)}</tr>")
    lines += ["</tbody>", "</table>"]
    return "\n".join(lines)


def versions() -> dict[str, str]:
    def first_line(cmd: list[str]) -> str:
        return subprocess.check_output(cmd, text=True).splitlines()[0]

    return {
        "crate": tomllib.loads((ROOT / "Cargo.toml").read_text())["package"]["version"],
        "rustc": first_line(["rustc", "--version"]),
        "gxx": first_line(["g++", "--version"]),
        "mpmath": mp.__version__,
        "platform": f"{platform.machine()} {platform.system()}",
    }


HEADER: Final = """\
# Accuracy

Relative error versus mpmath of xsf-rust {crate}, the C++ [xsf](https://github.com/scipy/xsf)
it is translated from (commit [`{short}`](https://github.com/scipy/xsf/tree/{commit})), and
the Rust standard library where it has the function. Where xsf lacks a function, SciPy {scipy}
or NumPy {numpy} takes its place.

- **peak**: the largest error; **gmean**: the geometric mean of the finite errors, each at least
  $`2^{{-53}}`$ (0 if all are 0)
- **\u221e (k)**: $`k`$ infinite errors, e.g. NaN, or a spurious overflow
- **bold**: the smallest peak or gmean of a row (unless all are equal at the displayed precision)
- **\u2757**: an error of {bad} or more, or infinite errors
- **\u2014**: in neither xsf nor SciPy (with the same definition); footnotes mark the rows with
  SciPy or NumPy instead of xsf
- **same**: the share of inputs where xsf-rust and xsf (or SciPy) return the same f64
  (bit-identical)
- **(log)**: sampled log-uniformly; **$`\\pm`$**: with a random sign

For reference, $`\\varepsilon`$ = 2.2e-16. See [Method](#method) for the details.
"""

METHOD: Final = """\
## Method

- Each domain has $`N`$ seeded pseudo-random f64 inputs, uniform unless marked "(log)".
- The reference $`f`$ is computed with mpmath {mpmath} at 50 or more digits, for the exact f64
  input, or with an alternative formula where mpmath is inaccurate (see
  `tools/accuracy/specs.py`).
- The error of $`y`$ is $`\\lvert y - f\\rvert / \\lvert f\\rvert`$ (of the complex values, for complex functions). Where
  a function oscillates, it is the smaller of that and $`\\lvert y - f\\rvert / A`$, with the local
  amplitude $`A`$ given below the table. A result is exact if both overflow to the same infinity,
  or both underflow ($`\\lvert f\\rvert`$ < 2.2e-308).
- **same** is rounded down, and counts any two NaNs as equal.
- Built on {platform}: xsf-rust with `{rustc}` (release profile), xsf with `{gxx}`
  (`-O2 -std=c++17`).
- To regenerate, run `just accuracy` (needs `cargo`, `g++`, `uv`, and internet access to
  download xsf). The references are cached in `tools/accuracy/.cache/`.
"""


def render(results: list[Result], ver: dict[str, str]) -> str:
    out = [
        HEADER.format(
            short=XSF_COMMIT[:7],
            commit=XSF_COMMIT,
            bad=f"$`10^{{{math.log10(BAD):.0f}}}`$",
            **scipy_funcs.VERSIONS,
            **ver,
        ).rstrip("\n")
    ]
    group = None  # the heading of a group without a table, until one of its subgroups has one
    for title, level in specs.SECTIONS.items():
        heading = f"{'#' * level} {title}"
        if level == 2:
            group = heading
        sec = [r for r in results if r.spec.section == title]
        if not sec:
            continue
        # the rows of a function together, in the order of their first row
        rank = {f: i for i, f in enumerate(dict.fromkeys(r.spec.func for r in sec))}
        sec.sort(key=lambda r: rank[r.spec.func])
        if level == 3 and group:
            out.append(group)
        group = None
        impls = ["rust", "cpp"] + (["std"] if any(r.spec.std for r in sec) else [])
        # the notes of this section, as footnote marks in order of first use, SciPy's last
        notes = dict.fromkeys(note for r in sec for note in r.spec.note.split())
        notes |= dict.fromkeys(note for r in sec for note in r.notes)
        marks = dict(zip(notes, FOOTNOTES[: len(notes)], strict=True))
        rows: list[list[str]] = []
        for r in sec:
            sup = "".join(marks[note] for note in r.notes)
            name = f"<code>{r.spec.func}</code>" + (f"<sup>{sup}</sup>" if sup else "")
            row = [name, r.domain, f"{r.n:,}"]
            for impl in impls:
                # SciPy in the xsf column, where xsf lacks the function
                s = r.stats.get(
                    "scipy" if impl == "cpp" and "scipy" in r.stats else impl
                )
                if s is None:
                    s = (NA, NA, 0) if impl == "cpp" else (None, None, 0)
                row += [fmt(s[0], s[2]), fmt(s[1])]
            if r.same is None:
                row.append(NA)
            else:
                row.append(
                    "100%" if r.same == 1 else f"{math.floor(1000 * r.same) / 10:.1f}%"
                )
            rows.append(row)
        # the second column: xsf, or SciPy or NumPy where xsf lacks the function
        used = {r.second for r in sec}
        second = " / ".join(s for s in (IMPLS["cpp"], "SciPy", "NumPy") if s in used)
        names = [(second or IMPLS[i]) if i == "cpp" else IMPLS[i] for i in impls]
        out += [heading, table(names, rows)]
        texts = specs.NOTES | scipy_funcs.NOTES
        out += [f"<sup>{m}</sup> {texts[note]}" for note, m in marks.items()]
    out.append(specs.NOT_COVERED.rstrip("\n"))
    out.append(METHOD.format(**ver).rstrip("\n"))
    return "\n\n".join(out) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None
    )
    ap.add_argument("--only", help="comma-separated function names; print to stdout")
    ap.add_argument(
        "--section", help="only the functions of this section; print to stdout"
    )
    ap.add_argument("--no-cache", action="store_true", help="recompute all references")
    ap.add_argument("--jobs", type=int, default=8, help="reference worker processes")
    ap.add_argument(
        "--n", type=int, help="at most N inputs per domain, for quick checks"
    )
    ap.add_argument(
        "--arb-check", action="store_true", help="check the references with Arb"
    )
    opts = ap.parse_args()
    t0 = time.time()
    mp.mp.dps = 30

    chosen = specs.SPECS
    if opts.only:
        chosen = [specs.BY_NAME[n] for n in opts.only.split(",")]
    elif opts.section:
        chosen = [s for s in chosen if s.section == opts.section]

    jobs: list[Job] = []
    for spec in chosen:
        for di, dom in enumerate(spec.domains):
            desc = dom.desc(spec.args)
            rng = random.Random(zlib.crc32(f"{spec.name}: {desc}".encode()))
            n = min(dom.n, opts.n) if opts.n else dom.n
            jobs.append(Job(spec, di, [dom.sample(rng) for _ in range(n)]))

    def has(job: Job, impl: str) -> bool:
        return getattr(job.spec, impl) and (
            impl != "cpp" or job.spec.domains[job.di].cpp
        )

    def lines(impl: str) -> list[tuple[str, Args]]:
        return [
            (job.spec.domains[job.di].call or job.spec.func, a)
            for job in jobs
            if has(job, impl)
            for a in job.args
        ]

    log(
        f"{sum(len(job.args) for job in jobs)} inputs, {len(jobs)} domains, {len(chosen)} functions"
    )
    rust = build_rust()
    outs: dict[str, Iterator[complex]] = {
        "rust": iter(run([*rust, "xsf"], lines("rust"))),
        "cpp": iter(run(build_cpp(), lines("cpp"))),
        "std": iter(run([*rust, "std"], lines("std"))),
    }

    refs = compute_refs(jobs, not opts.no_cache, opts.jobs)
    if opts.arb_check:
        arb_check(jobs, refs, opts.jobs)

    results: list[Result] = []
    worst: list[str] = []
    for job, rr in zip(jobs, refs, strict=True):
        spec, args = job.spec, job.args
        dom = spec.domains[job.di]
        ys = {
            impl: [next(it) for _ in args]
            for impl, it in outs.items()
            if has(job, impl)
        }
        if "cpp" not in ys and spec.name in scipy_funcs.FUNCS:
            ys["scipy"] = scipy_funcs.evaluate(spec.name, args)
        ok = {i: r for i, r in enumerate(rr) if isinstance(r, list)}
        bad = [r for r in rr if isinstance(r, str)]
        if bad:
            log(
                f"{spec.name} [{dom.desc(spec.args)}]: {len(bad)} references failed, e.g. {bad[0]}"
            )
        st: dict[str, Stats | None] = {}
        for impl, y in ys.items():
            errs = [(error(y[i], r), i) for i, r in ok.items()]
            st[impl] = stats([e for e, _ in errs])
            for e, i in sorted(errs, reverse=True)[:3]:
                worst.append(
                    f"{spec.name}{args[i]} {impl}: {y[i]!r} ref={ok[i][0]} err={e:.2e}"
                )
        same = None
        if (other := ys.get("cpp", ys.get("scipy"))) is not None:
            pairs = zip(ys["rust"], other, strict=True)
            same = sum(bits(y) == bits(z) for y, z in pairs) / len(args)
        results.append(Result(spec, dom.pretty(spec.args), len(ok), st, same))

    (CACHE / "worst.txt").write_text("\n".join(worst) + "\n")
    md = render(results, versions())
    if opts.only or opts.section or opts.n:
        print(md)
    else:
        (ROOT / "docs" / "accuracy.md").write_text(md)
        log(f"wrote {ROOT / 'docs' / 'accuracy.md'}")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
