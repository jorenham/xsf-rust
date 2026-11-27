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
from typing import Final
from unittest import mock

import flint
import fuzz_diff as fd
import mpmath as mp
import specs
from fuzz_diff import Line, Value

# a dump program that echoes the argument of each line, writes them as it ends or panics on
# `panic` (as a buffer that is flushed as it unwinds), fails on `unknown`, dies by a signal on
# `abort`, sleeps for `slow`, and counts its runs in `<path>.runs`
FAKE = """\
#!/bin/sh
echo >> "$0.runs"
out=
while read -r f x; do
  case $f in
    unknown) echo "unknown function: $f $x" >&2; exit 101 ;;
    panic) printf %s "$out"; echo "panicked" >&2; exit 101 ;;
    abort) kill -KILL $$ ;;
    slow) sleep 0.2 ;;
  esac
  out="$out$x
"
done
printf %s "$out"
"""


def kind(line: Line, base: Value, head: Value) -> str:
    return fd.judge(base, head, fd.ref(line))[0]


def conj(z: complex) -> complex:
    return z.conjugate()


c = complex
NAN, INF, PI = math.nan, math.inf, math.pi
EI = c(-743.8628562564797, -PI)
W = c(0.33408142401229419, -1.7585360826226356)
W1 = 0.5671432904097838
EXP = cmath.exp(c(709.9, PI / 4))
M = 0.5043435602314388
J = c(0.008577626689470264, 0.0)
DBL_MAX, DBL_MIN = fd.DBL_MAX, fd.DBL_MIN
SUB = math.nextafter(DBL_MIN, 0)

# the line, the results of the base and the head, and the kind
CASES: Final[list[tuple[Line, Value, Value, str]]] = [
    # the real functions return the Cauchy principal value
    (("elliprc", (1.0, -1.0)), 0.6232252401402306, NAN, "regressed"),
    (("elliprj", (1.0, 2.0, 3.0, -1.0)), -0.09324045243867635, NAN, "regressed"),
    # the real U rejects x < 0, also where it is a polynomial (and a panic is a NaN)
    (("hypu", (-1.0, 1.0, -2.0)), NAN, None, "both right"),
    # and any NaN argument
    (("elliprj", (1.0, NAN, 2.0, 3.0)), NAN, 1.0, "regressed"),
    # 0 log 0 = 0, where mpmath is NaN
    (("entr", (0.0,)), 0.0, NAN, "regressed"),
    (("xlogy", (0.0, 0.0)), 0.0, NAN, "regressed"),
    # and no reference where it is NaN otherwise
    (("xlog1py", (0.0, -1.0)), 0.0, NAN, "no ref"),
    # a NaN is better than a large error
    (("sinpi", (0.5,)), NAN, 0.5, "regressed"),
    (("sinpi", (0.5,)), NAN, INF, "both wrong"),
    # over- and underflow as for the reference rounded to f64
    (("hyp0f0#c", (709.9, PI / 4)), EXP, c(INF, INF), "regressed"),
    (("hypu", (-1.0, -(2.0**969), DBL_MAX)), DBL_MAX, INF, "regressed"),
    (("xlogy", (3.2101030212800104e-308, 2.0)), DBL_MIN, SUB, "regressed"),
    (("xlogy", (3.2101030212800104e-308, 2.0)), 2 * DBL_MIN, 0.0, "both wrong"),
    (("erfcinv", (0.0,)), INF, -INF, "regressed"),
    # an argument precondition, where the reference fails
    (("roots_gengegenbauer.1", (41.0, 33.0, -1.0, 1.0)), NAN, None, "no ref"),
    (("log1p#c", (-2.0, -0.0)), c(0, -PI), c(0, PI), "regressed"),
    # near a tiny real part, which `TINY` is not small next to
    (("expi#c", (-fd.TINY, -0.0)), EI, conj(EI), "regressed"),
    # the value at the point, which the sides approach slowly, as J_v(z) ~ (z / 2)^v
    (("bessel_j#c", (0.01, 0.0, 0.0)), 0j, c(5.84e-4, 9.2e-6), "regressed"),
    # I_v(z) ~ (z / 2)^v, still 1 at the offsets of 0 for a tiny v
    (("bessel_i#c", (fd.TINY, 0.0, 0.0)), 0j, c(1, 0), "regressed"),
    # no side of the pole at z = 1
    (("hyp2f1#c", (1.0, 1.0, 1.0, 1.0, 0.0)), c(INF, 0), c(0, INF), "no ref"),
    (("sinpi", (-0.0,)), -0.0, 0.0, "regressed"),
    # even in x, so of an unknown sign
    (("cosm1", (0.0,)), -0.0, 0.0, "both right"),
    (("hyp1f1#c", (1.0, 1.0, 0.0, 1e-100)), c(1, 1e-100), c(1, 0), "regressed"),
    # the sign of a zero imaginary part, also where Arb cannot tell it at the smaller offsets
    (("hyp1f1#c", (0.5, 1.5, -3.0, -0.0)), c(M, -0.0), c(M, 0.0), "regressed"),
    # but not their magnitudes, which are larger there
    (("hyp1f1#c", (2.0, 2.0**-52, 0.0, 0.0)), c(1, 0), c(1, 2 * DBL_MIN), "regressed"),
    # and where mpmath cannot
    (("lambertw#c", (1.0, -0.0)), c(W1, -0.0), c(W1, 0.0), "regressed"),
    # a zero of the function on the real axis, with that of its imaginary part
    (("cospi#c", (-2.5, -0.0)), c(0, -0.0), c(0, 0.0), "regressed"),
    (("cospi#c", (-2.5, -0.0)), c(0, -0.0), c(DBL_MAX, DBL_MAX), "regressed"),
    # a part below the accuracy of mpmath, of an unknown sign
    (("sph_bessel_j#c", (38.0, 35.38766539036541, 0.0)), J, J + 1e-80j, "both right"),
]


def test_verdicts() -> None:
    assert fd.differ(-0.0, 0.0, fd.RTOL)
    wrong = [(ln, b, h) for ln, b, h, want in CASES if kind(ln, b, h) != want]
    assert not wrong, wrong


def test_lambertw_cut() -> None:
    # W_0 below -1/e, also where a reference is the same on both sides of 2^-2148 (as mpmath),
    # or fails there
    value = fd._value

    def fails(name: str, args: tuple[float, ...]) -> object:
        # (an error string of `accuracy.ref_task`)
        return "timeout" if 0 < abs(args[-1]) < fd.TINY else value(name, args)

    for ctx in [
        contextlib.nullcontext(),
        mock.patch.dict(fd.ARB, clear=True),
        mock.patch.object(fd, "_value", fails),
    ]:
        with ctx:
            assert kind(("lambertw#c", (-2.5, -0.0)), W, conj(W)) == "regressed"


def test_references() -> None:
    # 1 - 2^-54, halfway between 1 and its predecessor, rounds to 1 (ties to even)
    r = fd.ref(("hyp1f1", (-1.0, 1.0, 2.0**-54)))
    assert isinstance(r, list) and isinstance(r[0], str)
    assert float(r[0]) == 1.0
    # U(a, a + 1, x) = x^-a, exactly 2^-6444, of more digits than Python writes an int with
    # (as 30 digits)
    line = ("hypu", (-6.0, -5.0, fd.TINY))
    r = fd.ref(line)
    assert isinstance(r, list) and isinstance(r[0], str)
    assert mp.almosteq(mp.mpf(r[0]), mp.ldexp(1, -6444), 1e-29, 0)
    with (
        mock.patch.object(fd, "_value", side_effect=ValueError),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        assert fd.ref(line) == "failed"


def test_precisions() -> None:
    # the last accurate precision, if a higher one is less accurate
    def f(*_: flint.arb) -> flint.acb:
        return flint.acb(1, flint.arb(0, 1e-30 if flint.ctx.prec < 512 else 1))

    with mock.patch.dict(fd.ARB, {"hyp1f1#c": f}):
        assert fd.ref(("hyp1f1#c", (1.0, 1.0, 1.0, 1.0))) == [["1", "0"], None]


def test_worker_precision() -> None:
    # (as in a worker process, which does not inherit it)
    with mp.workdps(15):
        fd.ref(("hyp1f1#c", (1.0, 1.0, 1.0, 1.0)))
        assert mp.mp.dps == specs.DPS


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
        # a death without unwinding
        died: list[Line] = [("a", (0.25,)), ("abort", (1.0,)), ("b", (2.0,))]
        assert fd.evaluate(exe, died) == "exit -9"
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
