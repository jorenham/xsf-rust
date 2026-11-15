# /// script
# requires-python = ">=3.15"
# dependencies = ["mpmath==1.4.1"]
# ///
"""Print the tables of the real `log1p` (src/xsf/log1p.rs).

Usage (from the repository root):

    uv run tools/log1p_tables.py

- `LN2_HI`, `LN2_LO`: ln(2) as a multiple of 2^-42 plus a remainder.
- `P`: log1p(r) = r - r^2 / 2 + r^3 P(r), interpolated at the Chebyshev nodes of [-R, R], where R
  is the largest |r| that `TABLE` leaves.
- `TABLE`: (invc, logc_hi, logc_lo) for each of the 64 intervals of m in [OFF, 2 OFF). invc is
  1 / m at the midpoint as a multiple of 2^-9 (exactly 1 in the interval of m = 1), and
  logc_hi + logc_lo = -ln(invc) with logc_hi a multiple of 2^-42.

The asserts are the conditions under which `rh` and `rh^2` in `log1p` are exact.
"""

import struct

import mpmath as mp

OFF = 0x3FE6A00000000000  # bits of the lower bound of m
MANTISSA_BITS = 52
TABLE_BITS = 6
SPLIT_BITS = 29  # low bits of m that are not in its head
INVC_BITS = 9
HI_BITS = 42
DEGREE = 6

mp.mp.prec = 300


def from_bits(bits: int) -> mp.mpf:
    return mp.mpf(struct.unpack("<d", struct.pack("<Q", bits))[0])


def hi_lo(v: mp.mpf) -> tuple[float, float]:
    hi = mp.nint(v * 2**HI_BITS) / 2**HI_BITS
    return float(hi), float(v - hi)


def table() -> tuple[list[tuple[float, float, float]], int, mp.mpf]:
    """The rows, the index of the interval around 1, and the largest |r| of the other rows."""
    shift = MANTISSA_BITS - TABLE_BITS
    i_one = (0x3FF0000000000000 - OFF) >> shift
    rows: list[tuple[float, float, float]] = []
    r_max = mp.mpf(0)
    for i in range(1 << TABLE_BITS):
        a, b = from_bits(OFF + (i << shift)), from_bits(OFF + ((i + 1) << shift))
        invc = mp.nint(2 ** (INVC_BITS + 1) / (a + b)) / 2**INVC_BITS
        if i == i_one:
            invc = mp.mpf(1)
        logc_hi, logc_lo = hi_lo(-mp.log(invc))
        r = max(abs(a * invc - 1), abs(b * invc - 1))

        # rh = head(m) invc - 1 is a multiple of `unit`, and must have at most 26 bits
        head_ulp = mp.mpf(2) ** (SPLIT_BITS - MANTISSA_BITS - (a < 1))
        invc_ulp = next(
            u for k in range(INVC_BITS + 1) if (invc / (u := mp.mpf(2) ** -k)) % 1 == 0
        )
        unit = head_ulp * invc_ulp
        assert r + 2 * head_ulp < unit * 2**26, i
        if i != i_one:
            # s = k LN2_HI + logc_hi + rh cannot cancel to less than rh^2 / 2
            assert abs(logc_hi) > r, i
            r_max = max(r_max, r)
        rows.append((float(invc), logc_hi, logc_lo))
    return rows, i_one, r_max


def poly(r_max: mp.mpf) -> tuple[list[float], mp.mpf]:
    """The coefficients of P, lowest degree first, and the relative error of log1p(r)."""

    def f(r: mp.mpf) -> mp.mpf:
        if abs(r) < mp.mpf(2) ** -40:
            # the middle node is 0 up to rounding
            return mp.mpf(1) / 3 - r / 4
        return (mp.log1p(r) - r + r * r / 2) / r**3

    n = DEGREE + 1
    nodes = [r_max * mp.cos(mp.pi * (2 * j + 1) / (2 * n)) for j in range(n)]
    vandermonde = mp.matrix([[x**p for p in range(n)] for x in nodes])
    coef = [
        float(c) for c in mp.lu_solve(vandermonde, mp.matrix([f(x) for x in nodes]))
    ]

    err = mp.mpf(0)
    m = 4000
    for j in range(m + 1):
        r = r_max * (2 * j - m) / m
        if r != 0:
            p = mp.polyval(coef[::-1], r)
            err = max(err, abs((r - r * r / 2 + r**3 * p) / mp.log1p(r) - 1))
    return coef, err


def main() -> None:
    rows, i_one, r_max = table()
    # the interval around 1 is not reduced: there r = x, up to its upper bound minus 1
    r_max = max(
        r_max, from_bits(OFF + ((i_one + 1) << (MANTISSA_BITS - TABLE_BITS))) - 1
    )
    r_max *= 1 + mp.mpf(2) ** -20
    coef, err = poly(r_max)
    ln2_hi, ln2_lo = hi_lo(mp.log(2))

    print(
        f"// |r| <= {float(r_max)!r}, relative error of log1p(r): 2^{float(mp.log(err, 2)):.2f}"
    )
    print(f"const OFF: u64 = {OFF:#018x};")
    print(f"const I_ONE: u64 = {i_one};")
    print(f"const LN2_HI: f64 = {ln2_hi!r};")
    print(f"const LN2_LO: f64 = {ln2_lo!r};")
    print(f"const P: [f64; {len(coef)}] = [")
    for c in coef:
        print(f"    {c!r},")
    print("];")
    print(f"const TABLE: [(f64, f64, f64); {len(rows)}] = [")
    for row in rows:
        print(f"    ({row[0]!r}, {row[1]!r}, {row[2]!r}),")
    print("];")


if __name__ == "__main__":
    main()
