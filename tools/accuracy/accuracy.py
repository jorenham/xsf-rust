"""Generate ACCURACY.md: compare xsf-rust, xsf (C++), and Rust std against mpmath.

Usage (from the repository root; or `just accuracy`):

    uv run --no-project --with mpmath==1.4.1 --with gmpy2==2.3.2 \
        python tools/accuracy/accuracy.py [--only NAME,...] [--no-cache]

Steps:

1. Generate the inputs of every domain in `specs.py` (seeded, so reproducible).
2. Evaluate them with xsf-rust and Rust std (`cargo run` of the dump program in this
   directory), and with xsf (C++), using `harness.cpp` compiled with `g++ -O2` against
   the xsf commit `XSF_COMMIT`, which is downloaded into `.cache/`.
3. Compute the mpmath references in parallel; these are cached in `.cache/refs/`, keyed
   by the inputs, the source of the reference function, and the mpmath version.
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
import os
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
    """Feed `<name> <args>...` lines to a dump program; return the f64 results."""
    inp = "".join(f"{name} {' '.join(map(repr, args))}\n" for name, args in lines)
    out = subprocess.run(cmd, input=inp, capture_output=True, text=True, check=True).stdout
    res = [float(s) for s in out.split()]
    assert len(res) == len(lines), (cmd, len(res), len(lines))
    return res


# ---------------------------------------------------------------- references


class Timeout(Exception):
    pass


def _alarm(*_):
    raise Timeout


def ref_task(task):
    """Evaluate one reference; returns [value, scale or None], or an error string."""
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
        if isinstance(v, mp.mpc):
            if abs(v.imag) > abs(v.real) * mp.mpf(10) ** -40:
                return f"complex result {v}"
            v = v.real
        return [mp.nstr(v, 30), None if scale is None else mp.nstr(scale, 10)]
    except Timeout:
        return "timeout"
    except Exception as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"
    finally:
        signal.alarm(0)


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
    """Cache key: the inputs, and the sources of the reference and of what it calls."""
    src = _sources(spec.ref.__code__, {"": inspect.getsource(spec.ref)})
    blob = json.dumps([args, sorted(src.items()), spec.dps, mp.__version__])
    return hashlib.sha256(blob.encode()).hexdigest()


def compute_refs(jobs, use_cache):
    """jobs: list of (spec, domain index, args list); returns a list of refs per job."""
    (CACHE / "refs").mkdir(parents=True, exist_ok=True)
    out, todo = [], []
    for j, (spec, di, args) in enumerate(jobs):
        path = CACHE / "refs" / f"{spec.name}-{di}.json"
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
        with multiprocessing.Pool(os.cpu_count()) as pool:
            res = []
            for i, r in enumerate(pool.imap(ref_task, tasks, chunksize=4)):
                res.append(r)
                if (i + 1) % 20000 == 0:
                    log(f"  {i + 1}/{len(tasks)} ({time.time() - t0:.0f} s)")
        k = 0
        for j in todo:
            spec, di, args = jobs[j]
            refs = res[k : k + len(args)]
            k += len(args)
            out[j] = refs
            path = CACHE / "refs" / f"{spec.name}-{di}.json"
            path.write_text(json.dumps({"key": ref_key(spec, args), "refs": refs}))
    return out


# ---------------------------------------------------------------- errors


def error(y, ref):
    """Relative error of `y` w.r.t. the reference (see the header of ACCURACY.md)."""
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
    return "nan" if math.isnan(x) else struct.pack("<d", x)


# ---------------------------------------------------------------- output


IMPLS = {"rust": "xsf-rust", "cpp": "xsf", "std": "std"}


def fmt(x, ninf=0):
    """Error in a fixed-width format, e.g. 1.1e-15, in monospace so that the digits align."""
    if x is None:
        return ""
    if x == 0:
        s = "0"
    elif not math.isfinite(x):
        s = f"{x} ({ninf})" if ninf else str(x)
    else:
        s = f"{x:.1e}"
    return f"<samp>{s}</samp>"


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
        cells += [f'<td align="right" nowrap>{x}</td>' for x in nums]
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

This table compares the accuracy of the f64 functions of xsf-rust {crate}, of the C++
[xsf](https://github.com/scipy/xsf) library that it is a translation of, and, where
available, of the Rust standard library. The xsf version is commit
[`{short}`](https://github.com/scipy/xsf/tree/{commit}).

**Method.** For each function and domain, N pseudo-random f64 inputs are drawn (with a fixed
seed) uniformly, or log-uniformly where marked with "(log)"; "\u00b1" means a random sign. Each
output y is compared to a reference value f, computed with [mpmath](https://mpmath.org/)
{mpmath} at 50 or more significant digits for the exact f64 input; where mpmath is known to
be inaccurate, the reference uses an alternative formula instead (see
`tools/accuracy/specs.py`). The error is the relative error |y \u2212 f| / |f|, except where
noted below a table. A result counts as exact if f overflows and y is an infinity of the same
sign, or if f underflows (|f| < 2.2e-308) and so does y. "peak" is the largest error, and
"rms" the root mean square of the finite errors; "inf (k)" means that k results have an
infinite error, e.g. NaN, or a spurious overflow. "bit-identical" is the fraction of the
inputs for which xsf-rust and xsf return the same f64 (rounded down, and counting any two NaNs
as equal). For reference, the machine epsilon is 2.2e-16.

Built on {platform} with:

- xsf-rust: `{rustc}`, release profile
- xsf: `{gxx}`, `-O2 -std=c++17`

To regenerate this file, run `just accuracy`. This needs `cargo`, `g++`, `uv`, and an
internet connection to download xsf. The references are cached in `tools/accuracy/.cache/`.
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
            row = [spec.name, r["domain"], str(r["n"])]
            for impl in impls:
                s = r["stats"].get(impl)
                if s is None:
                    row += ["", ""]
                else:
                    row += [fmt(s[0], s[2]), fmt(s[1])]
            same = r["same"]
            row.append("100%" if same == 1 else f"{math.floor(1000 * same) / 10:.1f}%")
            rows.append(row)
        out += [f"## {title}", table(impls, rows)]
        notes = {}
        for r in sec:
            for note in r["spec"].note.split():
                notes.setdefault(note, []).append(r["spec"].name)
        for note, names in notes.items():
            names = ", ".join(f"`{n}`" for n in dict.fromkeys(names))
            out.append(f"{names}: {specs.NOTES[note]}")
    return "\n\n".join(out) + "\n"


# ---------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", help="comma-separated function names; print to stdout")
    ap.add_argument("--no-cache", action="store_true", help="recompute all references")
    opts = ap.parse_args()
    t0 = time.time()
    mp.mp.dps = 30

    chosen = specs.SPECS
    if opts.only:
        names = opts.only.split(",")
        chosen = [specs.BY_NAME[n] for n in names]

    jobs = []
    for spec in chosen:
        for di, dom in enumerate(spec.domains):
            desc = dom.desc(spec.args)
            rng = random.Random(zlib.crc32(f"{spec.name}: {desc}".encode()))
            jobs.append((spec, di, [dom.sample(rng) for _ in range(dom.n)]))

    lines = [(spec.name, a) for spec, _, args in jobs for a in args]
    std_lines = [(spec.name, a) for spec, _, args in jobs if spec.std for a in args]
    log(f"{len(lines)} inputs, {len(jobs)} domains, {len(chosen)} functions")
    outs = {
        "rust": run(build_rust() + ["xsf"], lines),
        "cpp": run(build_cpp(), lines),
    }
    std_out = iter(run(build_rust() + ["std"], std_lines)) if std_lines else iter(())

    refs = compute_refs(jobs, not opts.no_cache)

    results, worst, k = [], [], 0
    for (spec, di, args), rr in zip(jobs, refs):
        dom = spec.domains[di]
        ys = {impl: outs[impl][k : k + len(args)] for impl in outs}
        if spec.std:
            ys["std"] = [next(std_out) for _ in args]
        k += len(args)
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
        same = sum(bits(ys["rust"][i]) == bits(ys["cpp"][i]) for i in range(len(args))) / len(args)
        domain = dom.pretty(spec.args)
        results.append(dict(spec=spec, domain=domain, n=len(ok), stats=st, same=same))

    (CACHE / "worst.txt").write_text("\n".join(worst) + "\n")
    md = render(results, versions())
    if opts.only:
        print(md)
    else:
        (ROOT / "ACCURACY.md").write_text(md)
        log(f"wrote {ROOT / 'ACCURACY.md'}")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
