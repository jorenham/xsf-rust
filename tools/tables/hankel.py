"""Real-axis references for the unscaled Hankel tables."""

from collections.abc import Callable

import flint
import mpmath as mp


def _axis[R: (mp.mpf, flint.arb)](
    kind: int,
    v: R,
    x: R,
    j: R,
    y: R,
    phase: Callable[[int], tuple[R, R]],
) -> tuple[R, R]:
    """H1 or H2 from J and Y at the positive order and argument.

    Combine order reflection and upper-cut continuation before adding the two parts
    (DLMF 10.4.6 and 10.11.5). In particular, H1 of negative order on the upper cut is
    exactly -J + iY, even when Y is many orders of magnitude larger than J.
    """
    s, c = phase(1)
    if x > 0:
        re, im = (j, y) if v >= 0 else (c * j - s * y, s * j + c * y)
        return re, im if kind == 1 else -im
    if v >= 0:
        if kind == 1:
            return -c * j + s * y, s * j + c * y
        return 3 * c * j - s * y, s * j - c * y
    if kind == 1:
        return -j, y
    s2, c2 = phase(2)
    return (2 + c2) * j - s2 * y, -s2 * j - c2 * y


def axis(kind: int, v: mp.mpf, x: mp.mpf) -> mp.mpc:
    """H1 or H2 on the real axis, with the upper limit on the cut."""
    u = abs(v)
    j, y = mp.besselj(u, abs(x)), mp.bessely(u, abs(x))
    return mp.mpc(
        *_axis(kind, v, x, j, y, lambda m: (mp.sinpi(m * u), mp.cospi(m * u)))
    )


def _arb_jy(u: flint.arb, x: flint.arb) -> tuple[flint.arb, flint.arb]:
    """J and Y at positive x and nonnegative u. When Y fails at a large
    noninteger u, expand J_(-u) and reflect its gamma factor (DLMF 10.2.2, 10.2.3, 5.5.3).
    """
    j, y = x.bessel_j(u), x.bessel_y(u)
    if y.is_finite():
        return j, y
    delta = (u - u.floor()).min(u.ceil() - u)
    if not delta > 0:
        return j, y
    # For every integer k >= 1, k |k-u| >= u delta / 2. Split at k = u/2.
    ratio = x * x / (2 * u * delta)
    if not ratio < 0.5:
        return j, y
    q = -x * x / 4
    term = total = flint.arb(1)
    k = 0
    cutoff = flint.arb(2) ** -flint.ctx.prec
    while abs(term) > cutoff:
        k += 1
        term *= q / (k * (k - u))
        total += term
    total += flint.arb(0, (abs(term) * ratio / (1 - ratio)).upper())
    leading = ((2 / x).log() * u + u.lgamma() - flint.arb.pi().log()).exp()
    return j, j * u.cos_pi() / u.sin_pi() - leading * total


def axis_arb(kind: int, v: flint.arb, x: flint.arb) -> flint.acb:
    """Arb version of axis."""
    u = abs(v)
    j, y = _arb_jy(u, abs(x))
    return flint.acb(
        *_axis(kind, v, x, j, y, lambda m: ((m * u).sin_pi(), (m * u).cos_pi()))
    )
