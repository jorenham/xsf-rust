# /// script
# requires-python = ">=3.12"
# dependencies = ["mpmath==1.4.1"]
# ///
"""Print the coefficient tables of the uniform asymptotic expansions of `pbdv` (src/xsf/par_cyl.rs).

Usage (from the repository root):

    uv run tools/pbdv_tables.py

- `U_POLY`, `V_POLY`: the polynomials u_s(t), v_s(t) of DLMF 12.10.9-13 for s <= 10.
- `TURNING_POINT`: the Taylor series in eps = t - 1 of zeta, phi(zeta), and A_s, B_s, C_s, D_s
  (DLMF 12.10.39-46) for s <= 3.

Everything is computed with exact rationals, and only rounded to f64 at the end.
"""

from collections.abc import Callable, Iterable
from fractions import Fraction as Fr
from math import comb, factorial

import mpmath as mp

SMAX = 10  # u_s, v_s for s <= SMAX
S = 3  # A_s, B_s, C_s, D_s for s <= S
N = 12  # Taylor terms

mp.mp.dps = 50

# polynomials (and truncated power series) are lists of coefficients, lowest degree first
type Poly = list[Fr]


def padd(*ps: Poly) -> Poly:
    n = max(map(len, ps))
    return [sum((p[i] for p in ps if i < len(p)), Fr(0)) for i in range(n)]


def pmul(p: Poly, q: Poly, n: int | None = None) -> Poly:
    n = len(p) + len(q) - 1 if n is None else n
    r = [Fr(0)] * n
    for i, a in enumerate(p[:n]):
        for j, b in enumerate(q[: n - i]):
            r[i + j] += a * b
    return r


def pscale(p: Poly, c: Fr | int) -> Poly:
    return [c * a for a in p]


def pder(p: Poly) -> Poly:
    return [i * p[i] for i in range(1, len(p))] or [Fr(0)]


def trim(p: Poly) -> Poly:
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p


T = [Fr(0), Fr(1)]  # t
T2M1 = [Fr(-1), Fr(0), Fr(1)]  # t^2 - 1


def solve_u(s: int, r: Poly) -> Poly:
    """The u of degree 3 s with (t^2 - 1) u' - 3 s t u = r (DLMF 12.10.12)."""

    # The t^j coefficient gives c[j + 1] = ((j - 1 - 3 s) c[j - 1] - r[j]) / (j + 1). For odd s,
    # u is odd, and for even s, c[0] is fixed by c[3 s] = 0.
    def chain(c0: Fr) -> Poly:
        c = [c0] + [Fr(0)] * (3 * s)
        c[1] = -r[0]
        for j in range(1, 3 * s):
            c[j + 1] = ((j - 1 - 3 * s) * c[j - 1] - (r[j] if j < len(r) else 0)) / (j + 1)
        return c

    c = chain(Fr(0))
    if s % 2 == 0:
        c1 = chain(Fr(1))
        c = chain(-c[-1] / (c1[-1] - c[-1]))
    lhs = padd(pmul(T2M1, pder(c)), pscale(pmul(T, c), -3 * s))
    assert trim(padd(lhs, pscale(r, -1))) == [0], s
    return c


def uv_polys(smax: int) -> tuple[list[Poly], list[Poly]]:
    """u_s and v_s (DLMF 12.10.10, 12.10.11, 12.10.13) for s <= smax."""
    us: list[Poly] = [[Fr(1)]]
    rs: list[Poly] = []
    for s in range(smax + 1):
        if s > 0:
            us.append(solve_u(s, rs[s - 1]))
        rp = rs[s - 1] if s > 0 else [Fr(0)]
        # 8 r_s = (3 t^2 + 2) u_s - 12 (s + 1) t r_{s-1} + 4 (t^2 - 1) r_{s-1}'
        r = padd(
            pmul([Fr(2), Fr(0), Fr(3)], us[s]),
            pscale(pmul(T, rp), -12 * (s + 1)),
            pscale(pmul(T2M1, pder(rp)), 4),
        )
        rs.append(trim(pscale(r, Fr(1, 8))))
    vs = [[Fr(1)]]
    for s in range(1, smax + 1):
        rm2 = rs[s - 2] if s >= 2 else [Fr(0)]
        vs.append(trim(padd(us[s], pscale(pmul(T, us[s - 1]), Fr(1, 2)), pscale(rm2, -1))))
    return [trim(p) for p in us], vs


def spow(f: Poly, p: Fr, n: int) -> Poly:
    """f^p as a power series, for f[0] = 1."""
    assert f[0] == 1
    g = [Fr(1)] + [Fr(0)] * (n - 1)
    for k in range(1, n):
        terms = ((p * j - (k - j)) * f[j] * g[k - j] for j in range(1, min(k + 1, len(f))))
        g[k] = sum(terms, Fr(0)) / k
    return g


def shift(p: Poly, n: int) -> Poly:
    """p(1 + eps) as a power series in eps."""
    r = [Fr(0)] * n
    for k, c in enumerate(p):
        for j in range(min(k + 1, n)):
            r[j] += c * comb(k, j)
    return r


def alpha(m: int) -> Fr:
    # DLMF 12.10.43
    num = Fr(1)
    for k in range(2 * m + 1, 6 * m, 2):
        num *= k
    return num / (factorial(m) * Fr(144) ** m)


def beta(m: int) -> Fr:
    return -Fr(6 * m + 1, 6 * m - 1) * alpha(m)


def comb_half(k: int) -> Fr:
    """binom(1/2, k)"""
    r = Fr(1)
    for j in range(k):
        r *= (Fr(1, 2) - j) / (j + 1)
    return r


def turning_point(us: list[Poly], vs: list[Poly]) -> dict[str, tuple[mp.mpf, Poly]]:
    """The Taylor series in eps of zeta / (2^(1/3) eps), phi, and A_s, B_s, C_s, D_s for s <= S,
    each with its irrational factor."""
    n = 3 * S + 4 + N
    # zeta / (2^(1/3) eps), with (2/3) zeta^(3/2) = int_1^t sqrt(u^2 - 1) du (12.10.39)
    s1 = [Fr(3, 2) * comb_half(k) / 2**k / (k + Fr(3, 2)) for k in range(n)]
    w = spow(s1, Fr(2, 3), n)
    r = pmul(w, spow([Fr(1), Fr(1, 2)], Fr(-1), n), n)  # zeta / (t^2 - 1), up to 2^(-2/3)
    p6 = pscale(spow(r, Fr(3, 2), n), Fr(1, 2))  # phi^6 (12.10.40), up to 2^-1
    p6pow = [[Fr(1)] + [Fr(0)] * (n - 1)]
    for _ in range(2 * S + 2):
        p6pow.append(pmul(p6pow[-1], p6, n))
    u = [shift(p, n) for p in us]
    v = [shift(p, n) for p in vs]

    two, third = mp.mpf(2), mp.mpf(1) / 3
    out = {"zeta": (mp.mpf(1), w[:N]), "phi": (two ** (-third / 2), spow(r, Fr(1, 4), N))}

    def coef(
        name: str,
        k: int,
        weights: Callable[[int], Fr],
        uv: list[Poly],
        sign: int,
        e: int,
        fac: mp.mpf | int,
    ) -> None:
        # 12.10.42 or 12.10.46: sign sum_m weights(m) phi^(6 (k - m)) uv[k - m] / zeta^e, which is
        # divisible by eps^e, times its factor 2^-fac
        ser = [Fr(0)] * n
        for m in range(k + 1):
            ser = padd(ser, pscale(pmul(p6pow[k - m], uv[k - m], n), sign * weights(m)))
        assert all(c == 0 for c in ser[:e]), name
        out[name] = (two**-fac, pmul(ser[e:], spow(w, Fr(-e), n), N))

    for s in range(S + 1):
        if s > 0:
            coef(f"A{s}", 2 * s, beta, u, 1, 3 * s, s)
        coef(f"B{s}", 2 * s + 1, alpha, u, -1, 3 * s + 2, s + 2 * third)
        coef(f"C{s}", 2 * s + 1, beta, v, -1, 3 * s + 1, s + third)
        if s > 0:
            coef(f"D{s}", 2 * s, alpha, v, 1, 3 * s, s)
    return out


def f64(x: float | Fr | mp.mpf) -> str:
    r = repr(float(x))
    if "e" in r:
        m, e = r.split("e")
        return f"{m if '.' in m else m + '.0'}e{int(e)}"
    return r


def table(name: str, typ: str, rows: Iterable[str]) -> None:
    print("#[rustfmt::skip]")
    print(f"const {name}: {typ} = [")
    for row in rows:
        print(f"    {row},")
    print("];")


def main() -> None:
    us, vs = uv_polys(SMAX)
    for name, ps in (("U_POLY", us), ("V_POLY", vs)):
        rows: list[str] = []
        for s, p in enumerate(ps):
            # u_s, v_s are even or odd, so in t^2 (times t for odd s), highest degree first
            assert all(c == 0 for c in p[1 - s % 2 :: 2])
            rows.append("&[" + ", ".join(f64(c) for c in reversed(p[s % 2 :: 2])) + "]")
        table(name, f"[&[f64]; {SMAX + 1}]", rows)
        print()

    tp = turning_point(us, vs)
    keys = [
        "zeta",
        "phi",
        *(f"{c}{s}" for c in "ABCD" for s in range(S + 1) if not (c in "AD" and s == 0)),
    ]
    rows = []
    for key in keys:
        fac, cs = tp[key]
        vals = (mp.mpf(c.numerator) / c.denominator * fac for c in reversed(cs))
        rows.append("[" + ", ".join(map(f64, vals)) + "]")
    table("TURNING_POINT", f"[[f64; {N}]; {len(keys)}]", rows)


if __name__ == "__main__":
    main()
