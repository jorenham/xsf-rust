"""The references, conventions, and tolerance exceptions of the tables that `generate.py`
replaces: an `Entry` per table and signature. The other tables keep those of xsref.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from typing import Final

import flint
import mpmath as mp
import specs

type Pred = Callable[[tuple[float, ...]], bool]  # of the f64 inputs
type Value = mp.mpf | tuple[mp.mpf, mp.mpf]  # a reference, or (reference, scale)
type Arb = Callable[..., flint.arb | flint.acb | None]  # None where Arb can't check
type RealArb = (
    Callable[[flint.arb], flint.arb] | Callable[[flint.arb, flint.arb], flint.arb]
)


@dataclass(frozen=True)
class Rule:
    """A convention: the result is `value`, with its class and sign, where `when` holds."""

    when: Pred
    value: complex
    why: str


@dataclass(frozen=True)
class Override:
    """Where `when` holds, xsf is wrong or much less accurate than xsf-rust, so the tolerance
    is at most the largest error of xsf-rust there."""

    when: Pred
    why: str


@dataclass(frozen=True)
class Out:
    """A logical output: real, or complex with a reference per part."""

    func: str  # the name in `tools/accuracy/src/main.rs` and `harness.cpp`
    refs: tuple[Callable[..., Value], ...]  # of the inputs as `mpf`s
    arbs: tuple[Arb, ...] = ()  # of the inputs as balls
    rules: tuple[Rule, ...] = ()
    overrides: tuple[Override, ...] = ()


@dataclass(frozen=True)
class Identity:
    """For an output without an Arb version: `residual` (of the inputs as `mpf`s, relative
    to 1) vanishes where `when` holds."""

    when: Pred
    residual: Callable[..., mp.mpf]


@dataclass(frozen=True)
class Entry:
    outs: tuple[Out, ...]
    identity: Identity | None = None


NAN: Final = Rule(lambda a: any(map(math.isnan, a)), math.nan, "NaN propagates")
INF: Final = Rule(lambda a: a[0] == math.inf, math.inf, "the limit for x -> +inf")
NEG_INF: Final = Rule(lambda a: a[0] == -math.inf, math.nan, "no limit for x -> -inf")
POLE: Final = Rule(
    lambda a: a[0] < 0 and a[0].is_integer(),
    math.nan,
    "a pole whose one-sided limits have opposite signs",
)
GAMMA: Final = Out(
    "gamma",
    (specs.r_gamma,),
    (specs.ARB["gamma"],),
    rules=(
        NAN,
        INF,
        NEG_INF,
        Rule(
            lambda a: a[0] == 0 and math.copysign(1, a[0]) > 0,
            math.inf,
            "the limit for x -> 0+, like C's tgamma(+0)",
        ),
        POLE,
    ),
)


def _ln_abs_gamma(x: flint.arb) -> flint.arb:
    """ln |Gamma(x)|, also for x < 0, where `arb.lgamma` is NaN."""
    return flint.acb(x).lgamma().real


LN_GAMMA_POLE: Final = Rule(
    lambda a: a[0] <= 0 and a[0].is_integer(),
    math.inf,
    "ln |Gamma(x)| -> +inf at a pole",
)


def _gammasgn_arb(x: flint.arb) -> flint.arb:
    return x.gamma().sgn()


def _lgam1p(x: mp.mpf) -> mp.mpf:
    """`specs.r_lgam1p`, but ln |Gamma(1 + x)| for x < -1, like xsf."""
    return mp.re(specs.r_lgam1p(x))


def _lgam1p_arb(x: flint.arb) -> flint.arb:
    return _ln_abs_gamma(1 + x)


def _digamma(x: mp.mpf) -> Value:
    """`specs.r_digamma`, with an absolute error where |psi(x)| < 1 for x < -1 (the accuracy
    that Cephes states), where the terms of the reflection cancel near a zero."""
    psi = specs.r_digamma(x)
    return (psi, mp.mpf(1)) if x < -1 else psi


def _limits(neg: float, pos: float) -> tuple[Rule, ...]:
    """NaN propagates, and the limits for x -> -inf and x -> +inf are `neg` and `pos`."""
    return (
        NAN,
        Rule(lambda a: a[0] == -math.inf, neg, "the limit for x -> -inf"),
        Rule(lambda a: a[0] == math.inf, pos, "the limit for x -> +inf"),
    )


# Arb loses about log2(x^2) bits in exp(+-x^2), so up to 2048 for an f64 x (and up to 1024 in
# exp(+-x))
EXP_SQUARE_BITS: Final = 2048


def _exp_square(
    f: Callable[[flint.arb], flint.arb],
) -> Callable[[flint.arb], flint.arb]:
    """`f` with `EXP_SQUARE_BITS` more bits."""

    def arb(x: flint.arb) -> flint.arb:
        with flint.ctx.extraprec(EXP_SQUARE_BITS):
            return f(x)

    return arb


EXP_SQUARE_OVERRIDE: Final = Override(
    lambda a: abs(a[0]) > 26.64,
    "xsf returns 0 for erfc and inf for erfi where exp(x^2) overflows, for |x| > 26.64",
)


def _ndtr(x: mp.mpf) -> mp.mpf:
    """`specs.r_ndtr` with the 2 log10 |x| more digits that x / sqrt(2) needs in
    exp(-x^2 / 2)."""
    with specs.wdps(2 * specs.xd(x)):
        return specs.r_ndtr(x)


def _ndtr_arb(x: flint.arb) -> flint.arb:
    """Q(1/2, x^2 / 2) / 2 for x <= 0 (DLMF 7.11.2), whose argument is exact, unlike that of
    erfc(-x / sqrt(2)), whose error Arb bounds absolutely."""
    q = (x * x / 2).gamma_upper(0.5, regularized=1) / 2
    return q if x <= 0 else 1 - q


def _erfcx_arb(x: flint.arb) -> flint.arb:
    return x.erfc() * (x * x).exp()


def _dawsn_arb(x: flint.arb) -> flint.arb:
    return flint.arb.pi().sqrt() / 2 * (-x * x).exp() * x.erfi()


def _ndtri_arb(p: flint.arb) -> flint.arb:
    return -flint.arb(2).sqrt() * (2 * p).erfcinv()


def _it2struve0_residual(x: mp.mpf) -> mp.mpf:
    """`r_it2struve0(x) - (pi/2 - int_0^x H_0(t)/t dt)` by quadrature, which checks the
    derivative and the value at 0. The limit at infinity is checked at x = 1e4, where
    `r_it2struve0` is below 1e-4."""
    if abs(specs.r_it2struve0(mp.mpf(10**4))) > 1e-4:
        return mp.inf

    def integrand(t: mp.mpf) -> mp.mpf:
        return mp.struveh(0, t) / t

    integral = mp.quad(integrand, mp.linspace(0, x, 2 + int(abs(x))))
    return abs(specs.r_it2struve0(x) - mp.pi / 2 + integral)


KELVIN: Final = ("ber", "bei", "ker", "kei", "berp", "beip", "kerp", "keip")


def _kelvin(kind: int) -> Callable[[mp.mpf], Value]:
    """`specs._kelvin`, and the limits for x -> 0+ (where it raises)."""

    def ref(x: mp.mpf) -> Value:
        if x == 0:
            return mp.mpf((1, 0, mp.inf, -mp.pi / 4, 0, 0, -mp.inf, 0)[kind])
        return specs._kelvin(kind, x)

    return ref


def _kelvin_arb(kind: int) -> Callable[[flint.arb], flint.arb | None]:
    """ber x + i bei x = J_0(x e^(3 pi i/4)) and ker x + i kei x = K_0(x e^(pi i/4)) (DLMF
    10.61.1, 10.61.2), with the derivatives -e^(3 pi i/4) J_1 and -e^(pi i/4) K_1."""

    def arb(x: flint.arb) -> flint.arb | None:
        if x == 0:
            return None
        k, order = kind % 4 >= 2, kind // 4
        w = flint.acb(0.25 if k else 0.75).exp_pi_i()
        # ber and bei are even
        z = w * abs(x)
        f = z.bessel_k(order) if k else z.bessel_j(order)
        if order:
            f *= -w * x.sgn()
        return f.imag if kind % 2 else f.real

    return arb


KE_NEGATIVE: Final = "Ke has a branch cut along the negative real axis"
KELVIN_OVERRIDE: Final = Override(
    lambda a: 3 <= abs(a[0]) < 18,
    "xsf loses up to 6 digits for 3 <= |x| < 18, where xsf-rust differs (scipy/xsf#232)",
)
KER_AT_0: Final = Rule(lambda a: a[0] == 0, math.inf, "ker x ~ -ln x for x -> 0+")


def _kelvin_out(kind: int) -> Out:
    rules = (KER_AT_0,) if kind == 2 else ()
    if kind % 4 >= 2:
        rules += (Rule(lambda a: a[0] < 0, math.nan, KE_NEGATIVE),)
    return Out(
        KELVIN[kind],
        (_kelvin(kind),),
        (_kelvin_arb(kind),),
        rules,
        (KELVIN_OVERRIDE,),
    )


def _kelvin_complex_out(j: int) -> Out:
    """Be, Ke, Be', or Ke' (j = 0 to 3)"""
    parts = (2 * j, 2 * j + 1)
    nan = complex(math.nan, math.nan)
    return Out(
        f"kelvin.{j}",
        tuple(_kelvin(k) for k in parts),
        tuple(_kelvin_arb(k) for k in parts),
        (Rule(lambda a: a[0] < 0, nan, KE_NEGATIVE),) if j % 2 else (),
        (KELVIN_OVERRIDE,),
    )


AIRY_OVERRIDES: Final = (
    Override(
        lambda a: a[0] > 103,
        "xsf returns 0 or NaN for x > 103.09, where AMOS under- or overflows (scipy/xsf#310)",
    ),
    Override(
        lambda a: a[0] < -10,
        "xsf rounds the phase 2/3 |x|^(3/2) for x < -10, and returns NaN for x < -1.05e6",
    ),
)


def _airy_out(i: int) -> Out:
    """Ai, Ai', Bi, or Bi' (i = 0 to 3), like `specs._airy` but without the functions that
    neither it nor its scale needs."""
    d, f = i % 2, (mp.airyai, mp.airybi)[i // 2]

    def ref(x: mp.mpf) -> Value:
        with specs.wdps(specs.xd(x) + 10):
            if x >= 0:
                return f(x, derivative=d)
            ai, bi = mp.airyai(x, derivative=d), mp.airybi(x, derivative=d)
            return (ai, bi)[i // 2], mp.hypot(ai, bi)

    def arb(x: flint.arb) -> flint.arb:
        return (x.airy_ai, x.airy_bi)[i // 2](derivative=d)

    return Out(f"airy[{i}]", (ref,), (arb,), overrides=AIRY_OVERRIDES)


def _real(
    func: str,
    ref: Callable[[mp.mpf], Value] | Callable[[mp.mpf, mp.mpf], Value],
    arb: RealArb,
    rules: tuple[Rule, ...] = (),
    overrides: tuple[Override, ...] = (),
) -> Entry:
    """A table with a real output."""
    return Entry((Out(func, (ref,), (arb,), rules, overrides),))


def _log1pmx(x: mp.mpf) -> mp.mpf:
    """`specs.r_log1pmx`, whose series does not end at x = 0."""
    return x if x == 0 else specs.r_log1pmx(x)


EXPIT_OVERRIDE: Final = Override(
    lambda a: a[0] < -709.78,
    "xsf returns 0 where exp(-x) overflows, for x < -709.78",
)
# of the last input
LOG1P_DOMAIN: Final = Rule(lambda a: a[-1] < -1, math.nan, "log1p is not real below -1")
# like in xsf, also where the logarithm is infinite
X_ZERO: Final = Rule(lambda a: a[0] == 0, 0.0, "x log y = 0 for x = 0")


def _sign_at(x: float, n: float) -> float:
    """The sign of x times (-1)^n, for an integer n."""
    return math.copysign(1, x) * (-1 if n % 2 else 1)


def _zeros(at: Pred, sign: Callable[[float], float], why: str) -> tuple[Rule, Rule]:
    """The zeros where `at` holds, with the sign `sign(x)`."""
    return (
        Rule(lambda a: at(a) and sign(a[0]) > 0, 0.0, why),
        Rule(lambda a: at(a) and sign(a[0]) < 0, -0.0, why),
    )


def _degrees(f: Callable[[flint.arb], flint.arb]) -> Callable[[flint.arb], flint.arb]:
    """`f` of x / 180, so that pi x / 180 has no rounding error in pi."""
    return lambda x: f(x / 180)


def _radian_arb(d: flint.arb, m: flint.arb, s: flint.arb) -> flint.arb:
    return (d + m / 60 + s / 3600) * flint.arb.pi() / 180


ZETA_POLE: Final = Rule(
    lambda a: a[0] == 1,
    math.inf,
    "the pole at x = 1, as +inf like in xsf (the limit for x -> 1- is -inf)",
)


def _zetac_arb(x: flint.arb) -> flint.arb:
    """zeta(x, 2) = zeta(x) - 1, without the cancellation for large x."""
    return x.zeta(2)


def _round_arb(x: flint.arb) -> flint.arb:
    """Python's `round`, which also rounds half to even, of the exact f64 input."""
    return flint.arb(round(float(x)))


ODD_AT_ZERO: Final = Rule(
    lambda a: a[0] == 0 and math.copysign(1, a[0]) > 0, 0.0, "f(+0) = +0 for an odd f"
)


def _bessel(n: int) -> dict[tuple[str, str], Entry]:
    """The tables of J_n, Y_n, I_n, e^-|x| I_n, K_n, and e^x K_n."""

    def j(x: mp.mpf) -> Value:
        """J_n(x), scaled by the Bessel modulus like `specs.r_jv` (at |x|: J_n is even or odd)."""
        return specs._osc(mp.besselj(n, x), n, abs(x))

    def at_n(f: Callable[[int, mp.mpf], Value]) -> Callable[[mp.mpf], Value]:
        return lambda x: f(n, x)

    def ie(x: flint.arb) -> flint.arb:
        return x.bessel_i(n) * (-abs(x)).exp()

    def ke(x: flint.arb) -> flint.arb:
        return x.bessel_k(n, scaled=True)

    def cut(f: str) -> tuple[Rule, ...]:
        why = f"{f}_{n} has a branch cut along the negative real axis"
        return (Rule(lambda a: a[0] < 0, math.nan, why),)

    def entry(
        f: str,
        ref: Callable[[mp.mpf], Value],
        rules: tuple[Rule, ...],
        arb: Arb | None = None,
    ) -> tuple[tuple[str, str], Entry]:
        """The table of `f` ("j", "ie", ...), with the Arb version of `specs.ARB` by default."""
        func = f"bessel_{f[0]}{n}{f[1:]}"
        out = Out(func, (ref,), (arb or specs.ARB[func],), rules)
        return (f"cyl_{func}", "d-d"), Entry((out,))

    odd = (ODD_AT_ZERO,) if n else ()
    return dict(
        [
            entry("j", j, odd),
            entry("y", at_n(specs.r_yv), cut("Y")),
            entry("i", at_n(specs.r_iv), odd),
            entry("ie", at_n(specs.r_ie), odd, ie),
            entry("k", at_n(specs.r_kv), cut("K")),
            entry("ke", at_n(specs.r_ke), cut("K"), ke),
        ]
    )


def _of_abs(f: Callable[[flint.arb], flint.arb]) -> Callable[[flint.arb], flint.arb]:
    """`f` of |x| with `EXP_SQUARE_BITS` more bits, like `specs.r_ci` and `specs.r_chi`."""

    def arb(x: flint.arb) -> flint.arb:
        return f(abs(x))

    return _exp_square(arb)


# also at x = -0
LOG_AT_0: Final = Rule(
    lambda a: a[0] == 0, -math.inf, "the limit for x -> 0, like ln |x|"
)


# with these more bits, 1 - x is exact and x - n pi keeps its precision for any f64 x
F64_BITS: Final = 1100
# 1 - m sin^2 phi cancels down to cos^2(fl(pi/2)) = 3.7e-33 for f64 inputs
CANCEL_DIGITS: Final = 40


def _ellipkm1_arb(p: flint.arb) -> flint.arb:
    with flint.ctx.extraprec(F64_BITS):
        return flint.acb(1 - p).elliptic_k().real


def _beyond_real(a: tuple[float, ...]) -> bool:
    """|phi| > asin(m^(-1/2)) for m > 1, where the incomplete integrals are not real."""
    phi, m = a
    return m > 1 and abs(phi) > math.asin(1 / math.sqrt(m))


def _incomplete(
    func: str,
    ref: Callable[[mp.mpf, mp.mpf], mp.mpf],
    arb: Callable[[flint.acb, flint.acb], flint.acb],
    letter: str,
    rules: tuple[Rule, ...],
) -> Entry:
    """An incomplete elliptic integral `letter`(phi|m), whose `rules` follow those for NaN and
    m > 1, and precede those for phi = +-inf."""

    def ref_digits(phi: mp.mpf, m: mp.mpf) -> mp.mpf:
        with specs.wdps(CANCEL_DIGITS):
            return ref(phi, m)

    def arb_real(phi: flint.arb, m: flint.arb) -> flint.arb:
        with flint.ctx.extraprec(F64_BITS):
            return arb(flint.acb(phi), flint.acb(m)).real

    f = f"{letter}(phi|m)"
    return _real(
        func,
        ref_digits,
        arb_real,
        (
            NAN,
            Rule(
                _beyond_real,
                math.nan,
                f"{f} is not real for m > 1, |phi| > asin(m^(-1/2))",
            ),
            Rule(
                lambda a: a[1] > 1,
                math.nan,
                f"{f} is real for m > 1, |phi| <= asin(m^(-1/2)), but NaN like in xsf",
            ),
            *rules,
            Rule(lambda a: a[0] == math.inf, math.inf, "the limit for phi -> +inf"),
            Rule(lambda a: a[0] == -math.inf, -math.inf, "the limit for phi -> -inf"),
        ),
    )


def _kolmogc_arb(x: flint.arb, r: int = 0) -> flint.arb:
    """1 - Q(x) = theta_4(0, tau) with the nome e^(i pi tau) = e^(-2 x^2), or its r-th
    derivative in z."""
    tau = flint.acb(0, 2 * x * x / flint.arb.pi())
    return flint.acb(0).modular_theta(tau, r)[3].real


def _kolmogorov_arb(x: flint.arb) -> flint.arb:
    return 1 - _kolmogc_arb(x)


def _kolmogp_arb(x: flint.arb) -> flint.arb:
    """Q'(x) = -(x / pi^2) theta_4''(0, tau), by the heat equation of theta_4."""
    return -x / flint.arb.pi() ** 2 * _kolmogc_arb(x, 2)


def _smirnov_jmax_arb(k: int, x: flint.arb) -> int | None:
    """`specs._smirnov_jmax` of the ball x, or None where x straddles a breakpoint (where the
    sum changes)."""

    def exact(e: flint.arb) -> Fraction:
        man, exp = e.mid().man_exp()
        return Fraction(int(man)) * Fraction(2) ** int(exp)

    lo, hi = (specs._smirnov_jmax(k, exact(e)) for e in (x.lower(), x.upper()))
    return lo if lo == hi else None


def _smirnov_arb(n: flint.arb, x: flint.arb, rest: bool = False) -> flint.arb:
    """`specs._smirnov`, or with `rest` the rest of its sum (`specs._smirnovc`); NaN, which
    Arb does not confirm, where x straddles a breakpoint."""
    k = int(float(n))
    jmax = _smirnov_jmax_arb(k, x)
    if jmax is None:
        return flint.arb.nan()
    js = range(jmax + 1, k + 1) if rest else range(jmax + 1)
    # the terms binom(k, j) a^(j - 1) b^(k - j), with a = x + j/k and b = 1 - a
    return x * sum(
        math.comb(k, j)
        * ((j + k * x) / k) ** (j - 1)
        * ((k - j - k * x) / k) ** (k - j)
        for j in js
    )


def _smirnovc_arb(n: flint.arb, x: flint.arb) -> flint.arb:
    """`specs._smirnovc`: 1 - smirnov(n, x), or the rest of its sum where that cancels."""
    c = 1 - _smirnov_arb(n, x)
    return c if c > 1e-3 else _smirnov_arb(n, x, rest=True)


def _smirnovp(n: mp.mpf, x: mp.mpf) -> mp.mpf:
    """`specs.r_smirnovp`, the left derivative, but the mean of the one-sided derivatives
    (which differ by 1) at the kink x = 1/n for n > 1, like in xsf."""
    return specs.r_smirnovp(n, x) + (0.5 if n > 1 and n * x == 1 else 0)


def _smirnovp_arb(n: flint.arb, x: flint.arb) -> flint.arb:
    """`_smirnovp`, the derivative of the sum of `_smirnov_arb` term by term."""
    k = int(float(n))
    jmax = _smirnov_jmax_arb(k, x)
    if jmax is None:
        return flint.arb.nan()
    d = -k * (1 - x) ** (k - 1) + (0.5 if k > 1 and k * x == 1 else 0)
    for j in range(1, jmax + 1):
        a, b = (j + k * x) / k, (k - j - k * x) / k
        t = a ** (j - 1) * b ** (k - j) + x * (j - 1) * a ** (j - 2) * b ** (k - j)
        d += math.comb(k, j) * (t - x * (k - j) * a ** (j - 1) * b ** (k - j - 1))
    return d


def _smirnovi(n: mp.mpf, q: mp.mpf) -> mp.mpf:
    """`specs.r_smirnovi`, but 1 - q exactly for n = 1 (smirnov(1, d) = 1 - d): where that is
    halfway between two f64s, a root found to 50 digits may round to either."""
    return 1 - q if n == 1 else specs.r_smirnovi(n, q)


def _inverse_arb(
    f: RealArb, ref: Callable[[mp.mpf], mp.mpf] | Callable[[mp.mpf, mp.mpf], mp.mpf]
) -> RealArb:
    """The root x of f(*a, x) = y for the inputs (*a, y), for a monotonic f: a ball around
    `ref(*a, y)` at 50 digits, which Arb confirms by the signs of f - y at its ends."""

    def arb(*args: flint.arb) -> flint.arb:
        y = args[-1]
        with mp.workdps(50):
            x = flint.arb(mp.nstr(ref(*(mp.mpf(float(t)) for t in args)), 50)).mid()
        r = abs(x) / 10**40
        # one unpacking, whose length pyrefly does not count against either signature
        lo, hi = (f(*(*args[:-1], x + s * r)) - y for s in (-1, 1))
        # a ball that Arb cannot confirm, unless f - y changes its sign
        return flint.arb(x, r if lo < 0 < hi or hi < 0 < lo else abs(x))

    return arb


def _unit_domain(why: str) -> Rule:
    """NaN where the last input is outside of [0, 1]."""
    return Rule(lambda a: not 0 <= a[-1] <= 1, math.nan, why)


def _ends(at_0: float, at_1: float, why: str) -> tuple[Rule, Rule]:
    """The values where the last input is 0 and where it is 1."""
    return (
        Rule(lambda a: a[-1] == 0, at_0, why),
        Rule(lambda a: a[-1] == 1, at_1, why),
    )


SMIRNOV_N: Final = Rule(lambda a: a[0] < 1, math.nan, "there is no sample for n < 1")


ENTRIES: Final[dict[tuple[str, str], Entry]] = {
    ("gamma", "d-d"): Entry((GAMMA,)),
    ("gammaln", "d-d"): Entry(
        (
            Out(
                "gammaln",
                (specs.r_gammaln,),
                (_ln_abs_gamma,),
                (LN_GAMMA_POLE,),
            ),
        )
    ),
    ("loggamma", "d-d"): Entry(
        (
            Out(
                "loggamma",
                (specs.r_loggamma,),
                (_ln_abs_gamma,),
                (
                    Rule(
                        lambda a: a[0] < 0,
                        math.nan,
                        "the principal ln Gamma(x) is not real for x < 0",
                    ),
                    LN_GAMMA_POLE,
                ),
            ),
        )
    ),
    ("rgamma", "d-d"): Entry(
        (Out("rgamma", (specs.r_rgamma,), (specs.ARB["rgamma"],), (NAN,)),)
    ),
    ("gammasgn", "d-d"): Entry(
        (
            Out(
                "gammasgn",
                (specs.r_gammasgn,),
                (_gammasgn_arb,),
                (
                    Rule(lambda a: a[0] == math.inf, 1.0, "the limit for x -> +inf"),
                    NEG_INF,
                    Rule(
                        lambda a: a[0] == 0 and math.copysign(1, a[0]) < 0,
                        -1.0,
                        "the sign of Gamma(-0) = -inf",
                    ),
                    POLE,
                ),
            ),
        )
    ),
    ("lgam1p", "d-d"): Entry((Out("lgam1p", (_lgam1p,), (_lgam1p_arb,)),)),
    ("digamma", "d-d"): Entry(
        (
            Out(
                "digamma",
                (_digamma,),
                (specs.ARB["digamma"],),
                (
                    INF,
                    NEG_INF,
                    Rule(
                        lambda a: a[0] == 0 and math.copysign(1, a[0]) > 0,
                        -math.inf,
                        "the limit for x -> 0+",
                    ),
                    POLE,
                ),
            ),
        )
    ),
    ("erf", "d-d"): Entry(
        (Out("erf", (specs.r_erf,), (flint.arb.erf,), _limits(-1, 1)),)
    ),
    ("erfc", "d-d"): Entry(
        (
            Out(
                "erfc",
                (specs.r_erfc,),
                (_exp_square(flint.arb.erfc),),
                _limits(2, 0),
                (EXP_SQUARE_OVERRIDE,),
            ),
        )
    ),
    ("erfi", "d-d"): Entry(
        (
            Out(
                "erfi",
                (specs.r_erfi,),
                (_exp_square(flint.arb.erfi),),
                _limits(-math.inf, math.inf),
                (EXP_SQUARE_OVERRIDE,),
            ),
        )
    ),
    ("erfcx", "d-d"): Entry(
        (Out("erfcx", (specs.r_erfcx,), (_erfcx_arb,), _limits(math.inf, 0)),)
    ),
    ("dawsn", "d-d"): Entry(
        (Out("dawsn", (specs.r_dawsn,), (_dawsn_arb,), _limits(-0.0, 0.0)),)
    ),
    ("erfcinv", "d-d"): Entry(
        (
            Out(
                "erfcinv",
                (specs.r_erfcinv,),
                (specs.ARB["erfcinv"],),
                (
                    Rule(
                        lambda a: not 0 <= a[0] <= 2,
                        math.nan,
                        "erfc takes values in [0, 2]",
                    ),
                    Rule(lambda a: a[0] == 0, math.inf, "erfc(x) -> 0 for x -> +inf"),
                    Rule(lambda a: a[0] == 2, -math.inf, "erfc(x) -> 2 for x -> -inf"),
                ),
            ),
        )
    ),
    ("ndtr", "d-d"): Entry((Out("ndtr", (_ndtr,), (_exp_square(_ndtr_arb),)),)),
    ("ndtri", "d-d"): Entry(
        (
            Out(
                "ndtri",
                (specs.r_ndtri,),
                (_ndtri_arb,),
                (
                    Rule(
                        lambda a: not 0 <= a[0] <= 1,
                        math.nan,
                        "ndtr takes values in [0, 1]",
                    ),
                    Rule(lambda a: a[0] == 0, -math.inf, "ndtr(x) -> 0 for x -> -inf"),
                    Rule(lambda a: a[0] == 1, math.inf, "ndtr(x) -> 1 for x -> +inf"),
                ),
            ),
        )
    ),
    ("it2struve0", "d-d"): Entry(
        (Out("it2struve0", (specs.r_it2struve0,)),),
        Identity(lambda a: math.isfinite(a[0]), _it2struve0_residual),
    ),
    **{(name, "d-d"): Entry((_kelvin_out(kind),)) for kind, name in enumerate(KELVIN)},
    ("kelvin", "d-cd_cd_cd_cd"): Entry(tuple(_kelvin_complex_out(j) for j in range(4))),
    ("airy", "d-d_d_d_d"): Entry(tuple(_airy_out(i) for i in range(4))),
    ("expm1", "d-d"): _real(
        "expm1", specs.r_expm1, flint.arb.expm1, _limits(-1, math.inf)
    ),
    ("exprel", "d-d"): _real(
        "exprel",
        specs.r_exprel,
        lambda x: x.expm1() / x,
        (
            *_limits(0.0, math.inf),
            Rule(lambda a: a[0] == 0, 1.0, "the limit for x -> 0"),
        ),
    ),
    ("log1p", "d-d"): _real(
        "log1p",
        specs.r_log1p,
        flint.arb.log1p,
        (
            INF,
            LOG1P_DOMAIN,
            Rule(lambda a: a[0] == -1, -math.inf, "the limit for x -> -1+"),
        ),
    ),
    ("log1pmx", "d-d"): _real(
        "log1pmx", _log1pmx, lambda x: x.log1p() - x, (LOG1P_DOMAIN,)
    ),
    ("exp2", "d-d"): _real("exp2", specs.r_exp2, lambda x: flint.arb(2) ** x),
    ("exp10", "d-d"): _real("exp10", specs.r_exp10, lambda x: flint.arb(10) ** x),
    ("expit", "d-d"): _real(
        "expit",
        specs.r_expit,
        lambda x: 1 / (1 + (-x).exp()),
        overrides=(EXPIT_OVERRIDE,),
    ),
    ("log_expit", "d-d"): _real(
        "log_expit", specs.r_log_expit, lambda x: -(-x).exp().log1p()
    ),
    ("logit", "d-d"): _real(
        "logit",
        specs.r_logit,
        lambda p: p.log() - (-p).log1p(),
        (
            Rule(
                lambda a: not 0 <= a[0] <= 1, math.nan, "expit takes values in [0, 1]"
            ),
            Rule(lambda a: a[0] == 0, -math.inf, "expit(x) -> 0 for x -> -inf"),
            Rule(lambda a: a[0] == 1, math.inf, "expit(x) -> 1 for x -> +inf"),
        ),
    ),
    ("cbrt", "d-d"): _real("cbrt", specs.r_cbrt, lambda x: x.sgn() * abs(x).root(3)),
    ("cosm1", "d-d"): _real("cosm1", specs.r_cosm1, lambda x: x.cos() - 1),
    ("xlogy", "d_d-d"): _real(
        "xlogy",
        specs.r_xlogy,
        lambda x, y: x * y.log(),
        (
            NAN,
            X_ZERO,
            Rule(lambda a: a[1] < 0, math.nan, "log y is not real for y < 0"),
            Rule(
                lambda a: a[0] > 0 and a[1] == math.inf,
                math.inf,
                "the limit for y -> +inf, with x > 0",
            ),
        ),
    ),
    ("xlog1py", "d_d-d"): _real(
        "xlog1py",
        specs.r_xlog1py,
        lambda x, y: x * y.log1p(),
        (
            NAN,
            X_ZERO,
            LOG1P_DOMAIN,
            Rule(
                lambda a: a[0] > 0 and a[1] == -1, -math.inf, "the limit for y -> -1+"
            ),
            Rule(lambda a: a[0] < 0 and a[1] == -1, math.inf, "the limit for y -> -1+"),
        ),
    ),
    ("sinpi", "d-d"): _real(
        "sinpi",
        specs.r_sinpi,
        flint.arb.sin_pi,
        _zeros(
            lambda a: a[0].is_integer(),
            lambda x: _sign_at(x, x),
            "sin(pi n) = +-0 with the sign of (-1)^n n, like in xsf",
        ),
    ),
    ("cospi", "d-d"): _real(
        "cospi",
        specs.r_cospi,
        flint.arb.cos_pi,
        (Rule(lambda a: a[0] % 1 == 0.5, 0.0, "cos(pi (n + 1/2)) = +0, like in xsf"),),
    ),
    ("sindg", "d-d"): _real("sindg", specs.r_sindg, _degrees(flint.arb.sin_pi)),
    ("cosdg", "d-d"): _real(
        "cosdg",
        specs.r_cosdg,
        _degrees(flint.arb.cos_pi),
        (
            Rule(
                lambda a: abs(a[0]) % 360 == 90,
                -0.0,
                "cos(x) = -0 for |x| = 90 + 360 n, like in xsf",
            ),
        ),
    ),
    ("tandg", "d-d"): _real(
        "tandg",
        specs.r_tandg,
        _degrees(flint.arb.tan_pi),
        _zeros(
            lambda a: a[0] % 180 == 0,
            lambda x: _sign_at(x, x // 180),
            "tan(180 n) = +-0 with the sign of (-1)^n n, like in xsf",
        ),
    ),
    ("cotdg", "d-d"): _real(
        "cotdg",
        specs.r_cotdg,
        _degrees(flint.arb.cot_pi),
        _zeros(
            lambda a: a[0] % 180 == 90,
            lambda x: _sign_at(1, (x - 90) // 180),
            "cot(90 + 180 n) = +-0 with the sign of (-1)^n, like in xsf",
        ),
        (
            Override(
                lambda a: not 45 <= abs(a[0]) % 180 <= 135,
                "xsf loses digits near the poles at 180 n, where xsf-rust differs "
                "(scipy/xsf#363)",
            ),
        ),
    ),
    ("radian", "d_d_d-d"): Entry((Out("radian", (specs.r_radian,), (_radian_arb,)),)),
    ("riemann_zeta", "d-d"): _real(
        "riemann_zeta",
        specs.r_riemann_zeta,
        flint.arb.zeta,
        (
            NAN,
            Rule(lambda a: a[0] == math.inf, 1.0, "the limit for x -> +inf"),
            Rule(
                lambda a: a[0] < 0 and a[0] % 2 == 0,
                0.0,
                "the trivial zeros at x = -2, -4, ..., as +0 like in xsf",
            ),
        ),
    ),
    ("zetac", "d-d"): _real(
        "zetac",
        specs.r_zetac,
        _zetac_arb,
        (
            NEG_INF,
            Rule(lambda a: a[0] == math.inf, 0.0, "the limit for x -> +inf"),
            ZETA_POLE,
        ),
    ),
    ("zeta", "d_d-d"): _real(
        "zeta", specs.r_hurwitz_zeta, flint.arb.zeta, (ZETA_POLE,)
    ),
    ("round", "d-d"): _real("round", mp.nint, _round_arb),
    **_bessel(0),
    **_bessel(1),
    ("exp1", "d-d"): _real(
        "exp1",
        specs.r_exp1,
        _exp_square(lambda x: x.expint(1)),
        (
            Rule(lambda a: a[0] < 0, math.nan, "E1(x) is not real for x < 0"),
            Rule(lambda a: a[0] == 0, math.inf, "the limit for x -> 0+, at x = +-0"),
        ),
    ),
    ("expi", "d-d"): _real(
        "expi", specs.r_expi, _exp_square(flint.arb.ei), (LOG_AT_0,)
    ),
    ("expn", "p_d-d"): _real(
        "expn",
        specs.r_expn,
        lambda n, x: x.expint(n),
        (
            Rule(lambda a: a[0] < 0, math.nan, "E_n is NaN for n < 0, like in xsf"),
            Rule(
                lambda a: a[1] < 0,
                math.nan,
                "NaN for x < 0, where E_n is not real for n >= 1",
            ),
        ),
    ),
    ("scaled_exp1", "d-d"): _real(
        "scaled_exp1",
        specs.r_scaled_exp1,
        _exp_square(lambda x: x * x.exp() * x.expint(1)),
        (
            Rule(lambda a: a[0] == 0, 0.0, "the limit for x -> 0+, at x = +-0"),
            Rule(lambda a: a[0] == math.inf, 1.0, "the limit for x -> +inf"),
        ),
    ),
    ("shichi", "d-d_d"): Entry(
        (
            Out("shichi.0", (specs.r_shi,), (_exp_square(flint.arb.shi),)),
            Out(
                "shichi.1",
                (specs.r_chi,),
                (_of_abs(flint.arb.chi),),
                (LOG_AT_0,),
            ),
        )
    ),
    ("sici", "d-d_d"): Entry(
        (
            Out(
                "sici.0",
                (specs.r_si,),
                (_exp_square(flint.arb.si),),
                _limits(-math.pi / 2, math.pi / 2)[1:],
            ),
            Out(
                "sici.1",
                (specs.r_ci,),
                (_of_abs(flint.arb.ci),),
                (
                    LOG_AT_0,
                    Rule(lambda a: a[0] == math.inf, 0.0, "the limit for x -> +inf"),
                    Rule(
                        lambda a: a[0] == -math.inf,
                        math.nan,
                        "NaN for x = -inf, like in xsf, not the limit 0 of Ci(|x|)",
                    ),
                ),
            ),
        )
    ),
    ("ellipk", "d-d"): _real(
        "ellipk",
        mp.ellipk,
        lambda m: flint.acb(m).elliptic_k().real,
        (
            Rule(lambda a: a[0] > 1, math.nan, "K(m) is not real for m > 1"),
            Rule(lambda a: a[0] == 1, math.inf, "K(m) -> +inf for m -> 1-"),
        ),
    ),
    ("ellipkm1", "d-d"): _real(
        "ellipkm1",
        specs.r_ellipkm1,
        _ellipkm1_arb,
        (
            NAN,
            Rule(lambda a: a[0] < 0, math.nan, "K(1 - p) is not real for p < 0"),
            Rule(lambda a: a[0] == 0, math.inf, "K(1 - p) -> +inf for p -> 0+"),
            Rule(lambda a: a[0] == math.inf, 0.0, "the limit for p -> +inf"),
        ),
    ),
    ("ellipe", "d-d"): _real(
        "ellipe",
        mp.ellipe,
        lambda m: flint.acb(m).elliptic_e().real,
        (
            NAN,
            Rule(lambda a: a[0] > 1, math.nan, "E(m) is not real for m > 1"),
            Rule(lambda a: a[0] == -math.inf, math.inf, "the limit for m -> -inf"),
        ),
    ),
    ("ellipeinc", "d_d-d"): _incomplete(
        "ellipeinc",
        mp.ellipe,
        flint.acb.elliptic_e_inc,
        "E",
        (
            Rule(
                lambda a: a[1] == -math.inf and 0 < a[0] < math.inf,
                math.inf,
                "the limit for m -> -inf, with 0 < phi < inf",
            ),
        ),
    ),
    ("kolmogorov", "d-d"): _real(
        "kolmogorov",
        specs.r_kolmogorov,
        _kolmogorov_arb,
        (NAN, Rule(lambda a: a[0] <= 0, 1.0, "Q(x) = P(K >= x) = 1 for x <= 0")),
    ),
    ("kolmogc", "d-d"): _real(
        "kolmogc",
        specs.r_kolmogc,
        _kolmogc_arb,
        (Rule(lambda a: a[0] <= 0, 0.0, "1 - Q(x) = P(K < x) = 0 for x <= 0"),),
    ),
    ("kolmogp", "d-d"): _real(
        "kolmogp",
        specs.r_kolmogp,
        _kolmogp_arb,
        (NAN, Rule(lambda a: a[0] <= 0, -0.0, "Q'(x) = -0 for x <= 0, like in xsf")),
    ),
    ("kolmogi", "d-d"): _real(
        "kolmogi",
        specs.r_kolmogi,
        _inverse_arb(_kolmogorov_arb, specs.r_kolmogi),
        (
            NAN,
            _unit_domain("kolmogorov takes values in [0, 1]"),
            Rule(lambda a: a[0] == 1, 0.0, "kolmogorov(x) -> 1 for x -> 0+"),
        ),
    ),
    ("kolmogci", "d-d"): _real(
        "kolmogci",
        specs.r_kolmogci,
        _inverse_arb(_kolmogc_arb, specs.r_kolmogci),
        (NAN, Rule(lambda a: a[0] == 0, 0.0, "kolmogc(x) -> 0 for x -> 0+")),
    ),
    ("smirnov", "p_d-d"): _real(
        "smirnov",
        specs.r_smirnov,
        _smirnov_arb,
        (
            NAN,
            SMIRNOV_N,
            _unit_domain("NaN for d outside of [0, 1], like in xsf"),
            *_ends(1.0, 0.0, "P(D_n^+ >= d) is 1 at d = 0 and 0 at d = 1"),
        ),
    ),
    ("smirnovc", "p_d-d"): _real(
        "smirnovc",
        specs.r_smirnovc,
        _smirnovc_arb,
        (NAN, *_ends(0.0, 1.0, "P(D_n^+ < d) is 0 at d = 0 and 1 at d = 1")),
    ),
    ("smirnovp", "p_d-d"): _real(
        "smirnovp",
        _smirnovp,
        _smirnovp_arb,
        (
            NAN,
            Rule(lambda a: a[1] == 0, -1.0, "the right derivative at d = 0 is -1"),
            Rule(
                lambda a: a[0] >= 2 and a[1] == 1,
                -0.0,
                "the left derivative at d = 1 is -0 for n >= 2, like in xsf",
            ),
        ),
    ),
    ("smirnovi", "p_d-d"): _real(
        "smirnovi",
        _smirnovi,
        _inverse_arb(_smirnov_arb, _smirnovi),
        (
            NAN,
            SMIRNOV_N,
            _unit_domain("smirnov takes values in [0, 1]"),
            *_ends(1.0, 0.0, "smirnov(n, d) is 0 only at d = 1 and 1 only at d = 0"),
        ),
    ),
    ("smirnovci", "p_d-d"): _real(
        "smirnovci",
        specs.r_smirnovci,
        _inverse_arb(_smirnovc_arb, specs.r_smirnovci),
        (
            NAN,
            *_ends(0.0, 1.0, "smirnovc(n, d) is 0 only at d = 0 and 1 only at d = 1"),
        ),
    ),
}
