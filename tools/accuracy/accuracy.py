"""Generate ACCURACY.md: compare xsf-rust, xsf (C++), and Rust std against mpmath.

Usage (from the repository root; or `just accuracy`):

    uv run --no-project --with mpmath==1.4.1 --with gmpy2==2.3.2 \
        python tools/accuracy/accuracy.py [--only NAME,...] [--no-cache] [--jobs N]

Steps:

1. Generate the inputs of every domain in `specs.py` (seeded, so reproducible).
2. Evaluate them with xsf-rust and Rust std (`cargo run` of the dump program in this
   directory), and with xsf (C++), using `harness.cpp` compiled with `g++ -O2` against
   the xsf commit `XSF_COMMIT`, which is downloaded into `.cache/`.
3. Compute the mpmath references in parallel (with `--jobs` worker processes, 8 by default);
   these are cached in `.cache/refs/`, keyed by the inputs, the source of the reference
   function, and the mpmath version.
4. Write the table to `ACCURACY.md` (or to stdout with `--only`), and the inputs with the
   largest errors to `.cache/worst.txt`.
"""

import argparse
import hashlib
import inspect
import io
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
import urllib.request
import zlib
from pathlib import Path

import mpmath as mp
import specs

XSF_COMMIT = "b6cb42595c449c93989b9ad5eecfe797a67ddacb"
XSF_URL = f"https://codeload.github.com/scipy/xsf/tar.gz/{XSF_COMMIT}"

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CACHE = HERE / ".cache"

DBL_MIN = 2.2250738585072014e-308
DBL_MAX = sys.float_info.max


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- evaluation


def fetch_xsf():
    inc = CACHE / f"xsf-{XSF_COMMIT}" / "include"
    if not inc.is_dir():
        log(f"downloading {XSF_URL}")
        with urllib.request.urlopen(XSF_URL) as r:
            data = r.read()
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
            tf.extractall(CACHE, filter="data")
    return inc


def build_cpp():
    inc = fetch_xsf()
    exe = CACHE / f"harness-{XSF_COMMIT[:12]}"
    src = HERE / "harness.cpp"
    log(f"compiling {src.name}")
    cmd = ["g++", "-O2", "-std=c++17", "-I", str(inc), str(src), "-o", str(exe)]
    subprocess.run(cmd, check=True)
    return [str(exe)]


def build_rust():
    manifest = HERE / "Cargo.toml"
    cmd = ["cargo", "build", "--release", "--quiet", "--manifest-path", str(manifest)]
    subprocess.run(cmd, check=True)
    cmd = ["cargo", "run", "--release", "--quiet", "--manifest-path", str(manifest), "--"]
    return cmd


def run(cmd, lines):
    """Feed `<name> <args>...` lines to a dump program; return the f64 (or complex) results."""
    inp = "".join(f"{name} {' '.join(map(repr, args))}\n" for name, args in lines)
    out = subprocess.run(cmd, input=inp, capture_output=True, text=True, check=True).stdout
    res = []
    for line in out.splitlines():
        re_im = [float(s) for s in line.split()]
        res.append(re_im[0] if len(re_im) == 1 else complex(*re_im))
    assert len(res) == len(lines), (cmd, len(res), len(lines))
    return res


# ---------------------------------------------------------------- references


class Timeout(Exception):
    pass


def _alarm(*_):
    raise Timeout


def ref_task(task):
    """Evaluate one reference; returns [value, scale or None], or an error string. The value
    is a string, or a [real, imaginary] pair of strings for a complex function."""
    name, args = task
    spec = specs.BY_NAME[name]
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(spec.timeout)
    try:
        mp.mp.dps = spec.dps
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


def timed_ref_task(task):
    t = time.perf_counter()
    return ref_task(task), time.perf_counter() - t


def _sources(code, seen):
    """Sources of the `specs` functions that `code` uses, recursively."""
    for name in code.co_names:
        f = getattr(specs, name, None)
        if inspect.isfunction(f) and name not in seen:
            seen[name] = inspect.getsource(f)
            _sources(f.__code__, seen)
    for c in code.co_consts:
        if inspect.iscode(c):
            _sources(c, seen)
    return seen


def ref_key(spec, args):
    """Cache key: the inputs, and the sources of the reference and of what it calls. For a
    lambda or closure, which can share its source with other references, the key also has its
    default arguments and captured variables (by name, for functions)."""
    src = _sources(spec.ref.__code__, {"": inspect.getsource(spec.ref)})
    key = [args, sorted(src.items()), spec.dps, mp.__version__]
    bound = list(spec.ref.__defaults__ or ()) + [c.cell_contents for c in spec.ref.__closure__ or ()]
    if bound:
        for f in bound:
            if inspect.isfunction(f) and getattr(specs, f.__name__, None) is f:
                _sources(f.__code__, src)
                src[f.__name__] = inspect.getsource(f)
        key[1] = sorted(src.items())
        key.append([getattr(f, "__name__", repr(f)) for f in bound])
    blob = json.dumps(key)
    return hashlib.sha256(blob.encode()).hexdigest()


def compute_refs(jobs, use_cache, nproc, tag=""):
    """jobs: list of (spec, domain index, args list); returns a list of refs per job."""
    (CACHE / "refs").mkdir(parents=True, exist_ok=True)
    out, todo = [], []
    for j, (spec, di, args) in enumerate(jobs):
        path = CACHE / "refs" / f"{spec.name}-{di}{tag}.json"
        key = ref_key(spec, args)
        if use_cache and path.exists():
            data = json.loads(path.read_text())
            if data["key"] == key:
                out.append(data["refs"])
                continue
        out.append(None)
        todo.append(j)
    tasks = [(jobs[j][0].name, a) for j in todo for a in jobs[j][2]]
    if tasks:
        log(f"computing {len(tasks)} references ({len(todo)} domains)")
        t0 = time.time()
        # recycle the workers, since mpmath caches grow at high precision
        with multiprocessing.Pool(nproc, maxtasksperchild=500) as pool:
            res, cost = [], {}
            for i, (r, dt) in enumerate(pool.imap(timed_ref_task, tasks, chunksize=4)):
                res.append(r)
                cost[tasks[i][0]] = cost.get(tasks[i][0], 0) + dt
                if (i + 1) % 20000 == 0:
                    log(f"  {i + 1}/{len(tasks)} ({time.time() - t0:.0f} s)")
        top = sorted(cost.items(), key=lambda kv: -kv[1])[:10]
        log("most expensive: " + ", ".join(f"{name} {dt:.0f} s" for name, dt in top))
        k = 0
        for j in todo:
            spec, di, args = jobs[j]
            refs = res[k : k + len(args)]
            k += len(args)
            out[j] = refs
            path = CACHE / "refs" / f"{spec.name}-{di}{tag}.json"
            path.write_text(json.dumps({"key": ref_key(spec, args), "refs": refs}))
    return out


# ---------------------------------------------------------------- errors


def error(y, ref):
    """Relative error of `y` w.r.t. the reference (see the header of ACCURACY.md)."""
    if isinstance(ref[0], list):
        return error_complex(y, mp.mpc(*ref[0]), ref[1])
    r = mp.mpf(ref[0])
    if math.isnan(y):
        return math.inf
    if mp.isinf(r) or abs(r) > DBL_MAX:
        return 0.0 if math.isinf(y) and (y > 0) == (r > 0) else math.inf
    if abs(r) < DBL_MIN:  # underflow
        return 0.0 if abs(y) < DBL_MIN else float(abs(y - r) / DBL_MIN)
    d = abs(mp.mpf(y) - r)
    e = d / abs(r)
    if ref[1] is not None:
        e = min(e, d / mp.mpf(ref[1]))
    return float(e)


def error_complex(y, r, scale=None):
    """|y - r| / |r|, with the same overflow, underflow, and scale rules as `error`."""
    if math.isnan(y.real) or math.isnan(y.imag):
        return math.inf
    if abs(r) > DBL_MAX:
        return 0.0 if math.isinf(abs(y)) else math.inf
    if math.isinf(abs(y)):
        return math.inf
    if abs(r) < DBL_MIN:
        return 0.0 if abs(y) < DBL_MIN else float(abs(y - r) / DBL_MIN)
    d = abs(mp.mpc(y) - r)
    e = d / abs(r)
    if scale is not None:
        e = min(e, d / mp.mpf(scale))
    return float(e)


def stats(errs):
    """(peak, rms, number of infinite errors); the rms is over the finite errors."""
    if not errs:
        return None
    fin = [e for e in errs if math.isfinite(e)]
    peak = max(fin, default=0.0)
    rms = peak * math.sqrt(math.fsum((e / peak) ** 2 for e in fin) / len(fin)) if peak else 0.0
    nbad = len(errs) - len(fin)
    return (math.inf if nbad else peak), (rms if fin else None), nbad


def bits(x):
    """The f64 bits, except that all NaNs are equal (the dumps do not keep NaN payloads)."""
    if isinstance(x, complex):
        return bits(x.real), bits(x.imag)
    return "nan" if math.isnan(x) else struct.pack("<d", x)


# ---------------------------------------------------------------- output


IMPLS = {"rust": "xsf-rust", "cpp": "xsf", "std": "std"}
NA = "\u2014"  # for functions that xsf does not have


def fmt(x, ninf=0):
    """Error in a fixed-width format, e.g. 1.1e-15, in monospace so that the digits align."""
    if x is None:
        return ""
    if x == NA:
        return NA
    if x == 0:
        s = "0"
    elif not math.isfinite(x):
        s = f"{x} ({ninf})" if ninf else str(x)
    else:
        s = f"{x:.1e}"
    return f"<samp>{s}</samp>"


def _err_key(cell):
    """Sort key of a formatted error, or None for "\u2014" and empty cells."""
    if not cell.startswith("<samp>"):
        return None
    s = cell.removeprefix("<samp>").removesuffix("</samp>")
    if s.startswith("inf"):
        return (math.inf, int(s.partition("(")[2].rstrip(")") or 0))
    return (float(s), 0)


def mark_best(errs):
    """Make the smallest peak and the smallest rms bold, if unique at the displayed precision;
    `errs` alternates peak and rms per implementation."""
    errs = list(errs)
    for j in range(2):
        keys = {c: _err_key(errs[c]) for c in range(j, len(errs), 2)}
        keys = {c: k for c, k in keys.items() if k is not None}
        if len(keys) > 1:
            best = min(keys.values())
            if list(keys.values()).count(best) == 1:
                c = next(c for c, k in keys.items() if k == best)
                errs[c] = f"<b>{errs[c]}</b>"
    return errs


def table(impls, rows):
    """HTML table with a two-level head; markdown tables cannot span columns."""
    head = [
        '<th rowspan="2">Function</th>',
        '<th rowspan="2">Domain</th>',
        '<th rowspan="2">N</th>',
    ]
    head += [f'<th colspan="2">{IMPLS[i]}</th>' for i in impls]
    head.append('<th rowspan="2">bit-identical</th>')
    lines = ["<table>", "<thead>", f"<tr>{''.join(head)}</tr>"]
    lines.append("<tr>" + "<th>peak</th><th>rms</th>" * len(impls) + "</tr>")
    lines += ["</thead>", "<tbody>"]
    for name, domain, *nums in rows:
        cells = [f"<td><code>{name}</code></td>", f"<td>{domain}</td>"]
        n, *errs, same = nums
        cells += [f'<td align="right" nowrap>{x}</td>' for x in [n, *mark_best(errs), same]]
        lines.append(f"<tr>{''.join(cells)}</tr>")
    lines += ["</tbody>", "</table>"]
    return "\n".join(lines)


def versions():
    def first_line(cmd):
        return subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()[0]

    crate = next(
        line.split('"')[1]
        for line in (ROOT / "Cargo.toml").read_text().splitlines()
        if line.startswith("version")
    )
    return {
        "crate": crate,
        "rustc": first_line(["rustc", "--version"]),
        "gxx": first_line(["g++", "--version"]),
        "mpmath": mp.__version__,
        "platform": f"{platform.machine()} {platform.system()}",
    }


HEADER = """\
# Accuracy

Relative error versus mpmath of xsf-rust {crate}, the C++ [xsf](https://github.com/scipy/xsf)
it is translated from (commit [`{short}`](https://github.com/scipy/xsf/tree/{commit})), and
the Rust standard library where it has the function.

- **peak**: the largest error; **rms**: the root mean square of the finite errors
- **inf (k)**: k infinite errors, e.g. NaN, or a spurious overflow
- **bold**: the smallest peak or rms of a row (if unique at the displayed precision)
- **\u2014**: not in xsf (mostly SciPy's Cython or Python code)
- **bit-identical**: the share of inputs where xsf-rust and xsf return the same f64
- **(log)**: sampled log-uniformly; **\u00b1**: with a random sign

For reference, \u03b5 = 2.2e-16. See [Method](#method) for the details.
"""

METHOD = """\
## Method

- Each domain has N seeded pseudo-random f64 inputs, uniform unless marked "(log)".
- The reference f is computed with mpmath {mpmath} at 50 or more digits, for the exact f64
  input, or with an alternative formula where mpmath is inaccurate (see
  `tools/accuracy/specs.py`).
- The error of y is |y \u2212 f| / |f| (of the complex values, for complex functions). Where
  a function oscillates, it is the smaller of that and |y \u2212 f| / A, with the local
  amplitude A given below the table. A result is exact if both overflow to the same infinity,
  or both underflow (|f| < 2.2e-308).
- bit-identical is rounded down, and counts any two NaNs as equal.
- Built on {platform}: xsf-rust with `{rustc}` (release profile), xsf with `{gxx}`
  (`-O2 -std=c++17`).
- To regenerate, run `just accuracy` (needs `cargo`, `g++`, `uv`, and internet access to
  download xsf). The references are cached in `tools/accuracy/.cache/`.
"""


def render(results, ver):
    out = [
        HEADER.format(short=XSF_COMMIT[:7], commit=XSF_COMMIT, **ver).rstrip("\n"),
    ]
    for title in specs.SECTIONS:
        sec = [r for r in results if r["spec"].section == title]
        if not sec:
            continue
        impls = ["rust", "cpp"] + (["std"] if any(r["spec"].std for r in sec) else [])
        rows = []
        for r in sec:
            spec = r["spec"]
            row = [spec.func, r["domain"], str(r["n"])]
            for impl in impls:
                s = r["stats"].get(impl)
                if s is None:
                    s = (NA, NA, 0) if impl == "cpp" else (None, None, 0)
                row += [fmt(s[0], s[2]), fmt(s[1])]
            same = r["same"]
            if same is None:
                row.append(NA)
            else:
                row.append("100%" if same == 1 else f"{math.floor(1000 * same) / 10:.1f}%")
            rows.append(row)
        out += [f"## {title}", table(impls, rows)]
        notes = {}
        for r in sec:
            for note in r["spec"].note.split():
                notes.setdefault(note, []).append(r["spec"].func)
        for note, names in notes.items():
            names = ", ".join(f"`{n}`" for n in dict.fromkeys(names))
            out.append(f"{names}: {specs.NOTES[note]}")
    out.append(specs.NOT_COVERED)
    out.append(METHOD.format(**ver).rstrip("\n"))
    return "\n\n".join(out) + "\n"


# ---------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", help="comma-separated function names; print to stdout")
    ap.add_argument("--section", help="only the functions of this section; print to stdout")
    ap.add_argument("--no-cache", action="store_true", help="recompute all references")
    ap.add_argument("--jobs", type=int, default=8, help="reference worker processes")
    ap.add_argument("--n", type=int, help="at most N inputs per domain, for quick checks")
    opts = ap.parse_args()
    t0 = time.time()
    mp.mp.dps = 30

    chosen = specs.SPECS
    if opts.only:
        names = opts.only.split(",")
        chosen = [specs.BY_NAME[n] for n in names]
    elif opts.section:
        chosen = [s for s in chosen if s.section == opts.section]

    jobs = []
    for spec in chosen:
        for di, dom in enumerate(spec.domains):
            desc = dom.desc(spec.args)
            rng = random.Random(zlib.crc32(f"{spec.name}: {desc}".encode()))
            n = min(dom.n, opts.n) if opts.n else dom.n
            jobs.append((spec, di, [dom.sample(rng) for _ in range(n)]))

    def has(spec, di, impl):
        return getattr(spec, impl) and (impl != "cpp" or spec.domains[di].cpp)

    def lines(impl):
        return [
            (spec.domains[di].call or spec.func, a)
            for spec, di, args in jobs
            if has(spec, di, impl)
            for a in args
        ]

    log(f"{sum(len(a) for *_, a in jobs)} inputs, {len(jobs)} domains, {len(chosen)} functions")
    outs = {
        "rust": iter(run(build_rust() + ["xsf"], lines("rust"))),
        "cpp": iter(run(build_cpp(), lines("cpp"))),
        "std": iter(run(build_rust() + ["std"], lines("std"))),
    }

    refs = compute_refs(jobs, not opts.no_cache, opts.jobs, f"-n{opts.n}" if opts.n else "")

    results, worst = [], []
    for (spec, di, args), rr in zip(jobs, refs):
        dom = spec.domains[di]
        ys = {impl: [next(it) for _ in args] for impl, it in outs.items() if has(spec, di, impl)}
        ok = [i for i, r in enumerate(rr) if isinstance(r, list)]
        bad = [r for r in rr if not isinstance(r, list)]
        if bad:
            log(f"{spec.name} [{dom.desc(spec.args)}]: {len(bad)} references failed, e.g. {bad[0]}")
        st = {}
        for impl, y in ys.items():
            errs = [(error(y[i], rr[i]), i) for i in ok]
            st[impl] = stats([e for e, _ in errs])
            for e, i in sorted(errs, reverse=True)[:3]:
                worst.append(
                    f"{spec.name}{tuple(args[i])} {impl}: {y[i]!r} ref={rr[i][0]} err={e:.2e}"
                )
        same = None
        if "cpp" in ys:
            same = sum(bits(y) == bits(z) for y, z in zip(ys["rust"], ys["cpp"])) / len(args)
        domain = dom.pretty(spec.args)
        results.append(dict(spec=spec, domain=domain, n=len(ok), stats=st, same=same))

    (CACHE / "worst.txt").write_text("\n".join(worst) + "\n")
    md = render(results, versions())
    if opts.only or opts.section or opts.n:
        print(md)
    else:
        (ROOT / "ACCURACY.md").write_text(md)
        log(f"wrote {ROOT / 'ACCURACY.md'}")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
