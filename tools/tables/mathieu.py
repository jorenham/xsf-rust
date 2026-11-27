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


def _value(v: mp.mpf | tuple[mp.mpf, mp.mpf]) -> mp.mpf:
    return v[0] if isinstance(v, tuple) else v


def radial(
    odd: int, kind: int, derivative: bool
) -> Callable[[mp.mpf, mp.mpf, mp.mpf], mp.mpf | tuple[mp.mpf, mp.mpf]]:
    """The modified Mathieu function of the given kind or its derivative."""

    def at(k: int, d: bool, m: mp.mpf, q: mp.mpf, x: mp.mpf):
        res = specs._mathieu_radial(odd, k, m, q, x, d)
        if k == 1 or not isinstance(res, tuple) or res[0] == 0:
            return res
        # The series is accurate relative to the scale. Where the value is much smaller (as
        # for large q at x = 0), more digits make it accurate relative to itself.
        extra = max(int(mp.log10(res[1] / abs(res[0]))), 0)
        if extra:
            with specs.wdps(extra + 5):
                res = specs._mathieu_radial(odd, k, m, q, x, d)
        return res

    def reference(m: mp.mpf, q: mp.mpf, x: mp.mpf) -> mp.mpf | tuple[mp.mpf, mp.mpf]:
        if x >= 0:
            return at(kind, derivative, m, q, x)
        # The first kind has definite parity: even for Mc, odd for Ms.
        flip = derivative != bool(odd)
        if kind == 1:
            res = at(1, derivative, m, q, -x)
            if not flip:
                return res
            return (-res[0], res[1]) if isinstance(res, tuple) else -res
        # Mc2(-z) is a combination of the solutions Mc1(z) and Mc2(z) of the same equation,
        # with the coefficient of Mc2 from the Wronskian (DLMF 28.20.22), and that of Mc1 from
        # z = 0. The series of Mc2 does not converge for large negative x, where Y_n(u) has a
        # tiny u.
        sign = 1 - 2 * odd
        c0 = _value(at(2, bool(odd), m, q, mp.mpf(0)))
        s0 = _value(at(1, bool(odd), m, q, mp.mpf(0)))
        a = (2 - 4 * odd) * c0 / s0
        y1 = _value(at(1, derivative, m, q, -x))
        y2 = at(2, derivative, m, q, -x)
        y = a * y1 - sign * _value(y2)
        y = -y if derivative else y
        if isinstance(y2, tuple):
            # the modulus of the first and second kind, where they oscillate
            return y, mp.hypot(-y1 if flip else y1, y)
        return y

    return reference


def radial1_residual(odd: int) -> Callable[[mp.mpf, mp.mpf, mp.mpf], mp.mpf]:
    """Check the cross-product reference against the single-Bessel expansions
    in DLMF 28.23.7, 9, 11, 13. The Bessel recurrence removes the coth singularity."""

    def residual(m: mp.mpf, q: mp.mpf, x: mp.mpf) -> mp.mpf:
        order = int(m)
        _, coef, _ = specs._mathieu_eig(odd, order, q)
        h = mp.sqrt(q)
        u = 2 * h * mp.sinh(x)
        weighted = odd == (order % 2 == 0)
        phase = (-1) ** (order // 2 + int(odd and order % 2 == 0))
        norm = mp.fsum(
            (-1) ** k * c * (n if weighted else 1) for k, (n, c) in enumerate(coef)
        )

        def coeff(n: int) -> tuple[mp.mpf, mp.mpf]:
            j = mp.besselj(n, u)
            return j, (mp.besselj(n - 1, u) - mp.besselj(n + 1, u)) / 2

        if weighted:
            pairs = [(c, coeff(n - 1), coeff(n + 1)) for n, c in coef]
            g = mp.fsum(c * (a[0] + b[0]) for c, a, b in pairs)
            dg = mp.fsum(c * (a[1] + b[1]) for c, a, b in pairs)
            y = phase * h * mp.cosh(x) * g / norm
            dy = phase * h * (mp.sinh(x) * g + 2 * h * mp.cosh(x) ** 2 * dg) / norm
        else:
            pairs = [(c, coeff(n)) for n, c in coef]
            y = phase * mp.fsum(c * a[0] for c, a in pairs) / norm
            dy = phase * 2 * h * mp.cosh(x) * mp.fsum(c * a[1] for c, a in pairs) / norm

        errors: list[mp.mpf] = []
        for j, actual in enumerate((y, dy)):
            if x == 0 and j != odd:
                errors.append(abs(actual))
                continue
            v = radial(odd, 1, bool(j))(m, q, x)
            expected, scale = v if isinstance(v, tuple) else (v, abs(v))
            errors.append(abs(actual - expected) / scale)
        return max(errors)

    return residual


def radial2_residual(odd: int) -> Callable[[mp.mpf, mp.mpf, mp.mpf], mp.mpf]:
    """The difference of the Wronskian of the first and second kind and its value 2 / pi
    (DLMF 28.20.22), relative to it. At x < 0, also that of the reference and the series of
    `specs._mathieu_radial`, relative to the reference or its scale."""

    def residual(m: mp.mpf, q: mp.mpf, x: mp.mpf) -> mp.mpf:
        y1, y2 = (
            [
                # the first kind has definite parity, so its value or derivative is 0 at x = 0
                mp.mpf(0)
                if x == 0 and kind == 1 and d == (not odd)
                else _value(radial(odd, kind, d)(m, q, x))
                for d in (False, True)
            ]
            for kind in (1, 2)
        )
        errors = [abs((y1[0] * y2[1] - y1[1] * y2[0]) * mp.pi / 2 - 1)]
        if x < 0:
            for d in (False, True):
                ref = radial(odd, 2, d)(m, q, x)
                y, scale = ref if isinstance(ref, tuple) else (ref, abs(ref))
                direct = specs._mathieu_radial(odd, 2, m, q, x, d)
                errors.append(abs(_value(direct) - y) / scale)
        return max(errors)

    return residual
