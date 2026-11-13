"""The references, conventions, and tolerance exceptions of the tables that `generate.py`
replaces: an `Entry` per table and signature. The other tables keep those of xsref.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Final

import flint
import mpmath as mp
import specs

type Pred = Callable[[tuple[float, ...]], bool]  # of the f64 inputs
type Value = mp.mpf | tuple[mp.mpf, mp.mpf]  # a reference, or (reference, scale)
type Arb = Callable[..., flint.arb | flint.acb | None]  # None where Arb can't check


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


# Arb loses about log2(x^2) bits in exp(+-x^2), so up to 2048 for an f64 x
EXP_SQUARE_BITS: Final = 2048


def _exp_square(f: Callable[[flint.arb], flint.arb]) -> Arb:
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
    arb: Callable[[flint.arb], flint.arb] | Callable[[flint.arb, flint.arb], flint.arb],
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
}
