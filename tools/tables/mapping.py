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


NAN: Final = Rule(lambda a: math.isnan(a[0]), math.nan, "NaN propagates")
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
        Rule(lambda a: a[0] == math.inf, math.inf, "the limit for x -> +inf"),
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
TAYLOR_RADIUS: Final = 0.2  # `LOGGAMMA_TAYLOR_RADIUS` of `src/xsf/loggamma.rs`
LN_GAMMA_ZEROS: Final = Override(
    lambda a: abs(a[0] - 1) < TAYLOR_RADIUS or 2 - TAYLOR_RADIUS < a[0] < 2,
    "xsf loses digits near the zeros at 1 and 2, where xsf-rust differs (scipy/xsf#529)",
)


def _gammasgn_arb(x: flint.arb) -> flint.arb:
    return x.gamma().sgn()


def _lgam1p(x: mp.mpf) -> mp.mpf:
    """`specs.r_lgam1p`, but ln |Gamma(1 + x)| for x < -1, like xsf."""
    return mp.re(specs.r_lgam1p(x))


def _lgam1p_arb(x: flint.arb) -> flint.arb:
    return _ln_abs_gamma(1 + x)


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


ENTRIES: Final[dict[tuple[str, str], Entry]] = {
    ("gamma", "d-d"): Entry((GAMMA,)),
    ("gammaln", "d-d"): Entry(
        (
            Out(
                "gammaln",
                (specs.r_gammaln,),
                (_ln_abs_gamma,),
                (LN_GAMMA_POLE,),
                (LN_GAMMA_ZEROS,),
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
                (LN_GAMMA_ZEROS,),
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
    ("it2struve0", "d-d"): Entry(
        (Out("it2struve0", (specs.r_it2struve0,)),),
        Identity(lambda a: math.isfinite(a[0]), _it2struve0_residual),
    ),
    **{(name, "d-d"): Entry((_kelvin_out(kind),)) for kind, name in enumerate(KELVIN)},
    ("kelvin", "d-cd_cd_cd_cd"): Entry(tuple(_kelvin_complex_out(j) for j in range(4))),
    ("airy", "d-d_d_d_d"): Entry(tuple(_airy_out(i) for i in range(4))),
}
