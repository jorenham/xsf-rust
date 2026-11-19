"""The spheroidal wave functions of the test tables, by the references of `specs` (Flammer's
method and normalization, like SciPy), with kd = 1 for the prolate and -1 for the oblate ones,
and the identities that check them.
"""

from collections.abc import Callable
from typing import Final

import mpmath as mp
import specs

type Pair = tuple[mp.mpf, mp.mpf]
type Ref = Callable[..., mp.mpf | Pair]  # of the inputs as `mpf`s

# the largest relative distance of a cv input to the characteristic value
CV_RTOL: Final = 1e-12


def segv(kd: int) -> Callable[[mp.mpf, mp.mpf, mp.mpf], mp.mpf]:
    """The characteristic value of the mode (m, n), which only depends on c^2."""
    return lambda m, n, c: specs._sw_cv(int(m), int(n), kd * c * c)[0]


def _rad_at_0(kf: int, m: int, n: int, c: mp.mpf, kd: int) -> Pair:
    """The oblate radial function of the first or second kind (kf) and its derivative at x = 0,
    where the terms of the series of the first kind are 0/0: continued from x = 1/2 by Taylor
    series, like the second kind for x < 2."""
    assert kd == -1
    if kf == 2:
        return specs._sw_rad(2, m, n, c, mp.mpf(0), kd)[:2]
    with specs.wdps(10 + int(c) // 2):
        cv, e, top = specs._sw_cv(m, n, kd * c * c)
        x = mp.mpf(1) / 2
        y, dy, _ = specs._sw_series(1, m, n, c, x, kd, e, top)
        # R = (x^2 + 1)^(m/2) g
        p = x * x + 1
        g = y / p ** (mp.mpf(m) / 2)
        dg = dy / p ** (mp.mpf(m) / 2) - m * x * g / p
        g, dg = specs._sw_taylor(m, c, kd, cv, x, g, dg, mp.mpf(0))
    # R has the parity of n - m, like its series of j_(m + r)(c x) / x^m with r = n - m mod 2
    return (g, mp.mpf(0)) if (n - m) % 2 == 0 else (mp.mpf(0), dg)


def _rad(kf: int, m: int, n: int, c: mp.mpf, x: mp.mpf, kd: int) -> Pair:
    """The radial function of the first or second kind (kf) and its derivative."""
    if x == 0:
        return _rad_at_0(kf, m, n, c, kd)
    return specs._sw_rad(kf, m, n, c, x, kd)[:2]


def radial(kf: int, kd: int, i: int) -> Ref:
    """The radial function of the first or second kind (kf) (i = 0) or its derivative (i = 1)
    of (m, n, c, x); with the scale of `specs._sw_radial` where they oscillate, except at
    x = 0."""

    def ref(m: mp.mpf, n: mp.mpf, c: mp.mpf, x: mp.mpf) -> mp.mpf | Pair:
        if x == 0:
            return _rad_at_0(kf, int(m), int(n), c, kd)[i]
        return specs._sw_radial(kf, m, n, c, x, kd)[i]

    return ref


def angular(kd: int, i: int) -> Ref:
    """The angular function of the first kind (i = 0) or its derivative (i = 1) of
    (m, n, c, x), with the scale of `specs._sw_angular` where they oscillate."""

    def ref(m: mp.mpf, n: mp.mpf, c: mp.mpf, x: mp.mpf) -> mp.mpf | Pair:
        return specs._sw_angular(m, n, c, x, kd)[i]

    return ref


def with_cv(ref: Ref, kd: int, any_cv_at_0: bool = False) -> Ref:
    """`ref` of (m, n, c, cv, x), for a cv input that is the characteristic value up to
    rounding: the function is that of the mode (m, n), and only defined for one. With
    `any_cv_at_0`, any cv is fine at x = 0, where the normalization fixes the function and the
    parity its derivative (the angular functions)."""

    def f(m: mp.mpf, n: mp.mpf, c: mp.mpf, cv: mp.mpf, x: mp.mpf) -> mp.mpf | Pair:
        if not (any_cv_at_0 and x == 0):
            exact = segv(kd)(m, n, c)
            if abs(cv - exact) > CV_RTOL * max(1, abs(exact)):
                raise ValueError(f"{cv} is not the characteristic value {exact}")
        return ref(m, n, c, x)

    return f


def _ode_residual(kd: int, m: int, n: int, c: mp.mpf, x: mp.mpf) -> mp.mpf:
    """The residuals of S' = dS/dx and of (p S')' + (cv - kd c^2 x^2 - m^2 / p) S = 0 with
    p = 1 - x^2 (DLMF 30.2.1), relative to their largest terms, for the angular function S of
    `specs._sw_ang`, with the derivatives by numerical differentiation. This checks the
    characteristic value, the recurrence, and the derivative, but not the normalization."""

    def s(i: int) -> Callable[[mp.mpf], mp.mpf]:
        return lambda t: specs._sw_ang(m, n, c, t, kd)[i]

    y, dy, cv = specs._sw_ang(m, n, c, x, kd)
    p = 1 - x * x
    terms = (
        p * mp.diff(s(1), x),
        -2 * x * dy,
        (cv - kd * (c * x) ** 2 - m * m / p) * y,
    )
    ode = abs(mp.fsum(terms)) / max(map(abs, terms))
    return max(ode, abs(mp.diff(s(0), x) - dy) / (abs(y) + abs(dy)))


def angular_residual(kd: int) -> Callable[..., mp.mpf]:
    """`_ode_residual` of the inputs (m, n, c, x) or (m, n, c, cv, x)."""

    def residual(m: mp.mpf, n: mp.mpf, c: mp.mpf, *a: mp.mpf) -> mp.mpf:
        return _ode_residual(kd, int(m), int(n), c, a[-1])

    return residual


def wronskian_residual(kd: int) -> Callable[..., mp.mpf]:
    """The residual of the Wronskian c (x^2 - kd) (R1 R2' - R1' R2) = 1 (DLMF 30.11.7, with
    gamma = c and z = x for the prolate, and gamma = -i c and z = i x for the oblate ones) of
    the inputs (m, n, c, x) or (m, n, c, cv, x)."""

    def residual(m: mp.mpf, n: mp.mpf, c: mp.mpf, *a: mp.mpf) -> mp.mpf:
        x = a[-1]
        (r1, d1), (r2, d2) = (_rad(kf, int(m), int(n), c, x, kd) for kf in (1, 2))
        return abs(c * (x * x - kd) * (r1 * d2 - d1 * r2) - 1)

    return residual
