"""Checks of the references and verdicts of `fuzz_diff.py` in known hard cases, and of the
handling of panics, unknown functions, and slow programs.

Usage (from the repository root):

    uv run --project tools tools/accuracy/test_fuzz_diff.py
"""

import cmath
import contextlib
import io
import math
import tempfile
from pathlib import Path
from unittest import mock

import flint
import fuzz_diff as fd
import mpmath as mp
import specs
from fuzz_diff import Line, Value

# a dump program that echoes the argument of each line, writes them as it ends or panics on
# `panic` (as a buffer that is flushed as it unwinds), fails on `unknown`, dies with only a part
# of them written on `abort`, sleeps for `slow`, and counts its runs in `<path>.runs`
FAKE = """\
#!/bin/sh
echo >> "$0.runs"
out=
while read -r f x; do
  case $f in
    unknown) echo "unknown function: $f $x" >&2; exit 101 ;;
    panic) printf %s "$out"; echo "panicked" >&2; exit 101 ;;
    abort) printf %s "${out%??}"; kill -KILL $$ ;;
    slow) sleep 0.2 ;;
  esac
  out="$out$x
"
done
printf %s "$out"
"""


def kind(line: Line, base: Value, head: Value) -> str:
    return fd.judge(base, head, fd.ref(line))[0]


def test_principal_values() -> None:
    # the real functions return the Cauchy principal value
    assert kind(("elliprc", (1.0, -1.0)), 0.6232252401402306, math.nan) == "regressed"
    line = ("elliprj", (1.0, 2.0, 3.0, -1.0))
    assert kind(line, -0.09324045243867635, math.nan) == "regressed"


def test_rejected_input() -> None:
    # U is not real for x < 0, unless it is a polynomial
    assert kind(("hypu", (0.5, 1.0, -2.0)), math.nan, 1.0) == "regressed"
    assert kind(("hypu", (0.5, 1.0, -2.0)), math.nan, None) == "both right"
    assert kind(("hypu", (-1.0, 1.0, -2.0)), -3.0, math.nan) == "regressed"
    assert kind(("hypu", (-2.0, 0.5, -2.0)), math.nan, 10.75) == "improved"
    # a NaN of mpmath, where the function has a limit
    assert kind(("entr", (0.0,)), 0.0, math.nan) == "no ref"
    assert kind(("xlogy", (0.0, 0.0)), 0.0, math.nan) == "no ref"


def test_branch_cut() -> None:
    line = ("log1p#c", (-2.0, -0.0))
    assert kind(line, complex(0, -math.pi), complex(0, math.pi)) == "regressed"
    assert kind(line, complex(0, math.pi), complex(0, -math.pi)) == "improved"
    # W_0 below -1/e, also where a reference is the same on both sides of 2^-2148 (as mpmath),
    # or fails there
    line = ("lambertw#c", (-2.5, -0.0))
    w = complex(0.33408142401229419, -1.7585360826226356)
    assert kind(line, w, w.conjugate()) == "regressed"
    with mock.patch.dict(fd.ARB, clear=True):
        assert kind(line, w, w.conjugate()) == "regressed"
    value = fd._value

    def fails(name: str, args: tuple[float, ...]) -> object:
        return "timeout" if 0 < abs(args[-1]) < fd.TINY else value(name, args)

    with mock.patch.object(fd, "_value", fails):
        assert kind(line, w, w.conjugate()) == "regressed"
    # near a tiny real part, which `TINY` is not small next to
    z = complex(-fd.DBL_MIN, -0.0)
    line = ("xlogy#c", (1.0, 0.0, z.real, z.imag))
    assert kind(line, cmath.log(z), cmath.log(z).conjugate()) == "regressed"
    line = ("expi#c", (-math.ulp(0), -0.0))
    ei = complex(-743.8628562564797, -math.pi)
    assert kind(line, ei, ei.conjugate()) == "regressed"


def test_branch_point() -> None:
    # the value at z = 1, which the sides approach slowly, as (1 - z)^0.001
    line = ("hyp2f1#c", (0.5, 0.5, 1.001, 1.0, 0.0))
    assert kind(line, 319.19260460853, complex(167.995, 0.2375)) == "regressed"
    # J_v(z) ~ (z / 2)^v
    line = ("bessel_j#c", (0.01, 0.0, 0.0))
    assert kind(line, 0j, complex(5.84e-4, 9.2e-6)) == "regressed"


def test_pole() -> None:
    # no side of the pole at z = 1
    line = ("hyp2f1#c", (1.0, 1.0, 1.0, 1.0, 0.0))
    assert kind(line, complex(math.inf, 0), complex(0, math.inf)) not in {
        "improved",
        "regressed",
    }


def test_zero_signs() -> None:
    assert fd.differ(-0.0, 0.0, fd.RTOL)
    assert kind(("sinpi", (-0.0,)), -0.0, 0.0) == "regressed"
    assert kind(("sinpi", (0.0,)), -0.0, 0.0) == "improved"
    # even in x, so of an unknown sign
    assert kind(("cosm1", (0.0,)), -0.0, 0.0) == "both right"
    line = ("hyp1f1#c", (1.0, 1.0, 0.0, 1e-100))
    assert kind(line, complex(1, 1e-100), complex(1, 0)) == "regressed"
    # e^z on the real axis, with the sign of the zero imaginary part
    line = ("hyp1f1#c", (1.0, 1.0, 2.0, -0.0))
    e2 = math.exp(2)
    assert kind(line, complex(e2, -0.0), complex(e2, 0.0)) == "regressed"
    # a zero of the function on the real axis, with that of its imaginary part
    for line in [("cospi#c", (-2.5, -0.0)), ("sinpi#c", (3.0, 0.0))]:
        assert kind(line, complex(0, -0.0), complex(0, 0.0)) == "regressed"


def test_midpoint() -> None:
    # 1 - 2^-54, halfway between 1 and its predecessor, rounds to 1 (ties to even)
    r = fd.ref(("hyp1f1", (-1.0, 1.0, 2.0**-54)))
    assert isinstance(r, list) and isinstance(r[0], str)
    assert float(r[0]) == 1.0


def test_large_exact() -> None:
    # U(a, a + 1, x) = x^-a, exactly 2^-6444, of more digits than Python writes an int with
    # (as 30 digits)
    line = ("hypu", (-6.0, -5.0, math.ulp(0)))
    r = fd.ref(line)
    assert isinstance(r, list) and isinstance(r[0], str)
    assert mp.almosteq(mp.mpf(r[0]), mp.ldexp(1, -6444), 1e-29, 0)
    with mock.patch.object(fd, "_value", side_effect=ValueError):
        assert fd.ref(line) == f"{fd.FAILED}ValueError: "


def test_precision() -> None:
    # the last accurate one, if a higher precision is less accurate
    def f(*_: flint.arb) -> flint.acb:
        return flint.acb(1, flint.arb(0, 1e-30 if flint.ctx.prec < 512 else 1))

    with mock.patch.dict(fd.ARB, {"hyp1f1#c": f}):
        assert fd.ref(("hyp1f1#c", (1.0, 1.0, 1.0, 1.0))) == [["1", "0"], None]


def test_documented_panic() -> None:
    # an argument precondition, where the reference fails
    line = ("roots_gengegenbauer.1", (41.0, 33.0, -1.0, 1.0))
    assert kind(line, math.nan, None) == "no ref"


def test_panics_and_unknown() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / "fake"
        exe.write_text(FAKE)
        exe.chmod(0o755)
        lines: list[Line] = [
            (f, (float(i),)) for i, f in enumerate(["a", "panic", "b", "panic"])
        ]
        runs = Path(f"{exe}.runs")
        assert fd.evaluate(exe, lines) == [0.0, None, 2.0, None]
        # one run per panic
        assert len(runs.read_text().splitlines()) == 2
        # an abort, after which the lines that it lost run again
        died: list[Line] = [("a", (0.25,)), ("abort", (1.0,)), ("b", (2.0,))]
        assert fd.evaluate(exe, died) == [0.25, None, 2.0]
        # a variant of the function that the program lacks, after others
        assert fd.evaluate(exe, [*lines, ("unknown", (4.0,))]) == "unsupported"
        # the time of one run for all lines, also after panics
        slow: list[Line] = [(f, (0.0,)) for f in ["slow", "panic"] * 3]
        with (
            mock.patch.object(fd, "RUN_TIMEOUT", 0),
            mock.patch.object(fd, "LINE_TIMEOUT", 0.05),
        ):
            assert fd.evaluate(exe, slow) == "timeout"


def test_report() -> None:
    with contextlib.redirect_stdout(io.StringIO()):
        assert fd.report([]) is None
        paths = [fd.report(["x"]) for _ in range(2)]
    assert None not in paths and paths[0] != paths[1]
    for p in paths:
        if p is not None:
            p.unlink()


if __name__ == "__main__":
    mp.mp.dps = specs.DPS
    fd.CACHE.mkdir(parents=True, exist_ok=True)
    for name, f in list(globals().items()):
        if name.startswith("test_") and callable(f):
            f()
            print(f"{name}: ok")
