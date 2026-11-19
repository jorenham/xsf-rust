"""The references of the Mathieu tables, from those of `specs`, and the identities that check
them. mpmath and Arb have no Mathieu functions.
"""

import math
from collections.abc import Callable

import mpmath as mp
import specs


def characteristic(odd: int) -> Callable[[mp.mpf, mp.mpf], mp.mpf]:
    """a_m(q), or b_m(q) if `odd`."""
    return (specs.r_mathieu_a, specs.r_mathieu_b)[odd]


def angular(
    odd: int, derivative: bool
) -> Callable[[mp.mpf, mp.mpf, mp.mpf], tuple[mp.mpf, mp.mpf]]:
    """ce_m(z, q), or se_m(z, q) if `odd`, or its z-derivative, at z = x degrees, with its
    scale."""
    return [
        (specs.r_mathieu_cem, specs.r_mathieu_cemp),
        (specs.r_mathieu_sem, specs.r_mathieu_semp),
    ][odd][derivative]


def cv_residual(odd: int) -> Callable[[mp.mpf, mp.mpf], mp.mpf]:
    """The difference of the characteristic value and eigenvalue floor((m - odd) / 2) of the
    recurrence matrix of `specs._mathieu_matrix`, from `mp.eigsy`, relative to max(|a|, 1). For
    m <= 30 and |q| <= 100, where the Fourier coefficients past 50 more rows are negligible."""

    def residual(m: mp.mpf, q: mp.mpf) -> mp.mpf:
        i = int(m - odd) // 2
        d, e, _ = specs._mathieu_matrix(odd, int(m), q, i + 50)
        a = mp.diag(d)
        for k in range(len(d) - 1):
            a[k, k + 1] = a[k + 1, k] = e[k]
        lam = sorted(mp.eigsy(a, eigvals_only=True))[i]
        cv = characteristic(odd)(m, q)
        return abs(lam - cv) / max(abs(cv), 1)

    return residual


def angular_residual(odd: int) -> Callable[[mp.mpf, mp.mpf, mp.mpf], mp.mpf]:
    """The residuals of y' = dy/dz and of the Mathieu equation y'' + (a - 2 q cos 2z) y = 0
    (DLMF 28.2.1) for y = ce_m (or se_m, if `odd`) at z = x degrees, relative to the scale of
    `specs._mathieu_angular`, with the derivatives by `mp.diff`."""

    def residual(m: mp.mpf, q: mp.mpf, x: mp.mpf) -> mp.mpf:
        def dz(derivative: bool) -> mp.mpf:
            def g(t: mp.mpf) -> mp.mpf:
                return angular(odd, derivative)(m, q, t)[0]

            return mp.diff(g, x) * 180 / mp.pi

        v = characteristic(odd)(m, q) - 2 * q * mp.cos(2 * mp.radians(x))
        w = max(abs(v), 1)
        y, scale = angular(odd, False)(m, q, x)
        yp = angular(odd, True)(m, q, x)[0]
        return max(abs(dz(False) - yp) / mp.sqrt(w), abs(dz(True) + v * y) / w) / scale

    return residual


def km(a: tuple[float, ...]) -> float:
    """`qm(|q|) + m / 2` of `mathieu_km_too_large` (src/xsf/mathieu.rs), the number of Fourier
    coefficients that specfun predicts."""
    q = abs(a[1])
    s = math.sqrt(q)
    if q <= 1:
        qm = 7.5 + 56.1 * s - 134.7 * q + 90.7 * s * q
    else:
        qm = 17.0 + 3.1 * s - 0.126 * q + 0.0037 * s * q
    return qm + 0.5 * a[0]
