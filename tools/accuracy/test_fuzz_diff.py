"""Checks of the references and verdicts of `fuzz_diff.py` in known hard cases, and of the
handling of panics, unknown functions, and slow programs.

Usage (from the repository root):

    uv run --project tools tools/accuracy/test_fuzz_diff.py
"""

import contextlib
import io
import math
import tempfile
from pathlib import Path
from unittest import mock

import fuzz_diff as fd
import mpmath as mp
import specs
from fuzz_diff import Line, Value

# a dump program that echoes the argument of each line, and panics on `panic`, fails on
# `unknown`, and sleeps for `slow`
FAKE = """\
#!/bin/sh
while read -r f x; do
  case $f in
    unknown) echo "unknown function: $f $x" >&2; exit 101 ;;
    panic) echo "panicked" >&2; exit 101 ;;
    slow) sleep 0.2 ;;
  esac
  echo "$x"
done
"""


def kind(line: Line, base: Value, head: Value) -> str:
    return fd.judge(base, head, fd.ref(line))[0]


def test_principal_values() -> None:
    # the real functions return the Cauchy principal value
    assert kind(("elliprc", (1.0, -1.0)), 0.6232252401402306, math.nan) == "regressed"
    line = ("elliprj", (1.0, 2.0, 3.0, -1.0))
    assert kind(line, -0.09324045243867635, math.nan) == "regressed"


def test_rejected_input() -> None:
    # the real function is NaN for x < 0, even where U is a polynomial
    assert kind(("hypu", (-1.0, 1.0, -2.0)), math.nan, -3.0) == "regressed"
    assert kind(("hypu", (-1.0, 1.0, -2.0)), math.nan, None) == "both right"


def test_branch_cut() -> None:
    line = ("log1p#c", (-2.0, -0.0))
    assert kind(line, complex(0, -math.pi), complex(0, math.pi)) == "regressed"
    assert kind(line, complex(0, math.pi), complex(0, -math.pi)) == "improved"


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
    line = ("hyp1f1#c", (1.0, 1.0, 0.0, 1e-100))
    assert kind(line, complex(1, 1e-100), complex(1, 0)) == "regressed"
    # e^z on the real axis, with the sign of the zero imaginary part
    line = ("hyp1f1#c", (1.0, 1.0, 2.0, -0.0))
    e2 = math.exp(2)
    assert kind(line, complex(e2, -0.0), complex(e2, 0.0)) == "regressed"


def test_midpoint() -> None:
    # 1 - 2^-54, halfway between 1 and its predecessor, rounds to 1 (ties to even)
    r = fd.ref(("hyp1f1", (-1.0, 1.0, 2.0**-54)))
    assert isinstance(r, list) and isinstance(r[0], str)
    assert float(r[0]) == 1.0


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
        assert fd.evaluate(exe, lines) == [0.0, None, 2.0, None]
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
