"""The references, conventions, and tolerance exceptions of the tables that `generate.py`
replaces: an `Entry` per table and signature. The other tables keep those of xsref.
"""

import functools
import math
import sys
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from typing import Final, Protocol

import flint
import mpmath as mp
import specs
import spheroidal

import accuracy

type Pred = Callable[[tuple[float, ...]], bool]  # of the f64 inputs
type Value = mp.mpf | tuple[mp.mpf, mp.mpf]  # a reference, or (reference, scale)
type Ball = flint.arb | flint.acb
type Arb = Callable[..., Ball | None]  # None where Arb can't check
type RealArb = (
    Callable[[flint.arb], flint.arb]
    | Callable[[flint.arb, flint.arb], flint.arb]
    | Callable[[flint.arb, flint.arb, flint.arb], flint.arb]
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
    is at most the largest error of xsf-rust there, and a convention need not pass xsref's
    test."""

    when: Pred
    why: str


@dataclass(frozen=True)
class Floor:
    """Where `when` holds, the rounding errors of xsf-rust can exceed those of xsf by more than
    the test allows, because it uses another method, or because they depend on the platform's
    libm, so the tolerance is at least `tol`."""

    when: Pred
    tol: float
    why: str


@dataclass(frozen=True)
class Out:
    """A logical output: real, or complex with a reference per part."""

    func: str  # the name in `tools/accuracy/src/main.rs` and `harness.cpp`
    refs: tuple[Callable[..., Value], ...]  # of the inputs as `mpf`s
    arbs: tuple[Arb, ...] = ()  # of the inputs as balls
    rules: tuple[Rule, ...] = ()
    overrides: tuple[Override, ...] = ()
    floors: tuple[Floor, ...] = ()


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
COMPLEX_NAN: Final = Rule(NAN.when, complex(math.nan, math.nan), NAN.why)
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


def _nonpos_int(x: float) -> bool:
    return math.isfinite(x) and x <= 0 and x.is_integer()


COMPLEX_POLE: Final = Rule(
    lambda a: a[1] == 0 and _nonpos_int(a[0]),
    complex(math.nan, math.nan),
    "a pole at z = 0, -1, -2, ..., whose limit depends on the direction",
)


def _sum_nonpos_int(a: tuple[float, ...]) -> bool:
    """Whether the exact sum of the two inputs is a non-positive integer."""
    if not all(map(math.isfinite, a)):
        return False
    s = sum(map(Fraction, a))
    return s <= 0 and s.denominator == 1


def _gamma_pole(a: tuple[float, ...]) -> bool:
    """Whether Gamma(a) or Gamma(b) is a pole that Gamma(a + b) does not cancel, or both are."""
    poles = tuple(map(_nonpos_int, a))
    return all(poles) or any(poles) and not _sum_nonpos_int(a)


def _sums_exact(
    f: Callable[[flint.arb, flint.arb], flint.arb],
) -> Callable[[flint.arb, flint.arb], flint.arb]:
    """`f` with the bits that make sums of its arguments exact, even of very different size."""

    def arb(a: flint.arb, b: flint.arb) -> flint.arb:
        with flint.ctx.extraprec(2 * max(abs(math.frexp(float(x))[1]) for x in (a, b))):
            return f(a, b)

    return arb


def _gamma_ratio(a: flint.arb, m: flint.arb) -> flint.arb:
    """Gamma(a + m) / Gamma(a), also where both are poles, by reflection (DLMF 5.5.3)."""
    s = a + m
    if a.is_integer() and a <= 0 and s.is_integer() and s <= 0:
        return (-1) ** (int(float(m)) % 2) * (1 - a).gamma() * (1 - s).rgamma()
    return s.gamma() * a.rgamma()


@_sums_exact
def _beta_arb(a: flint.arb, b: flint.arb) -> flint.arb:
    """Gamma(a) Gamma(b) / Gamma(a + b), with the first factor that is not a pole in front."""
    if a.is_integer() and a <= 0:
        a, b = b, a
    return a.gamma() * _gamma_ratio(a + b, -a)


def _betaln_arb(a: flint.arb, b: flint.arb) -> flint.arb:
    return abs(_beta_arb(a, b)).log()


@_sums_exact
def _binom_arb(n: flint.arb, k: flint.arb) -> flint.arb:
    return (k + 1).rgamma() * _gamma_ratio(n + 1 - k, k)


BETA_POLES: Final = (
    Rule(
        lambda a: all(map(_nonpos_int, a)) and a[1] % 2 == 0,
        math.inf,
        "both arguments are poles, and xsf returns (-1)^b inf",
    ),
    Rule(
        lambda a: all(map(_nonpos_int, a)) and a[1] % 2 == 1,
        -math.inf,
        "both arguments are poles, and xsf returns (-1)^b inf",
    ),
    Rule(_gamma_pole, math.inf, "a pole whose one-sided limits have opposite signs"),
)


# the Lanczos approximation (g = 6.024680040776729583740234375) of Gamma, scaled by exp(g): the
# coefficients of the rational function in descending powers, as in xsf
LANCZOS_NUM: Final = (
    6.061842346248907e-3,
    5.098416655656676e-1,
    19.519927882476175,
    449.9445569063168,
    6955.999602515376,
    75999.29304014542,
    601859.6171681099,
    3481712.154980646,
    14605578.087685067,
    43338889.32467614,
    86363131.2881386,
    103794043.11634454,
    56906521.913471565,
)
LANCZOS_DENOM: Final = (
    1.0,
    66.0,
    1925.0,
    32670.0,
    357423.0,
    2637558.0,
    13339535.0,
    45995730.0,
    105258076.0,
    150917976.0,
    120543840.0,
    39916800.0,
    0.0,
)


# `MAXGAM` of src/utils.rs, above which Gamma overflows, `MAXLGM` of gamma.rs, above which
# lgam does, and `BETA_ASYMP_FACTOR` of beta.rs
MAXGAM: Final = 171.6243769563027
MAXLGM: Final = 2.556348e305
ASYMP_FACTOR: Final = 1e6


SMALL_B: Final = (
    "for |b| < 1 and a > 4, xsf-rust computes ln Gamma(a) - ln Gamma(a + b) with lgam_ratio "
    "(scipy/xsf#357), which avoids the rounding error of a + b that xsf amplifies"
)
STIRLING: Final = (
    "for a, b > 0 and a + b > 171.6, xsf-rust uses Stirling's formula instead of "
    "ln Gamma(a) + ln Gamma(b) - ln Gamma(a + b) (scipy/xsf#357), and so rounds differently"
)


def _by_size(args: tuple[float, ...]) -> tuple[float, float]:
    """The arguments of `beta` and `lbeta`, the larger |a| first."""
    a, b = sorted(args, key=abs, reverse=True)
    return a, b


def _asymptotic(a: float, b: float) -> bool:
    """Whether `lbeta_asymp` applies."""
    return a > ASYMP_FACTOR * abs(b) and a > ASYMP_FACTOR


def _small_b_beta(args: tuple[float, ...]) -> bool:
    """Whether `beta` takes the `small_b` path: |b| < 1 and a > 4 for the larger |a|."""
    a, b = _by_size(args)
    return abs(b) < 1 and 4 < a < math.inf


def _small_b(args: tuple[float, ...]) -> bool:
    """Whether `lbeta` takes the `small_b` path, which the asymptotic expansion precedes."""
    return _small_b_beta(args) and not _asymptotic(*_by_size(args))


def _stirling(args: tuple[float, ...]) -> bool:
    """Whether `beta` takes the Stirling branch of `lbeta_ratio`: a, b > 0, b >= a / 8, and
    a + b > MAXGAM, after `small_b` and the asymptotic expansion."""
    a, b = _by_size(args)
    return (
        b > 0
        and b >= a / 8
        and a + b > MAXGAM
        and not _small_b_beta(args)
        and not _asymptotic(a, b)
    )


def _beta_huge_negint(args: tuple[float, ...]) -> bool:
    """Whether xsf's `beta` takes an integer below -2^31 as a pole, which does not fit in an
    `int`, although the limit is finite: for the first argument that is a non-positive integer
    (xsf checks a before b), and the other one an integer 0 < y <= -x."""
    a, b = args
    x, y = (a, b) if _nonpos_int(a) else (b, a)
    return _nonpos_int(x) and x < -(2.0**31) and y.is_integer() and 0 < y <= -x


def _no_pole(args: tuple[float, ...]) -> bool:
    return not any(map(_nonpos_int, args))


def _beta_float_sum_pole(args: tuple[float, ...]) -> bool:
    """Whether the rounded a + b is a pole of Gamma, but their exact sum is not, and neither
    a nor b is (xsf then has 1 / Gamma(a + b) = 0, or ln Gamma(a + b) = inf)."""
    return (
        _no_pole(args) and _nonpos_int(args[0] + args[1]) and not _sum_nonpos_int(args)
    )


BETA_OVERRIDES: Final = (
    Override(
        _beta_huge_negint,
        "xsf takes the integers below -2^31, which do not fit in an `int`, as poles "
        "(scipy/xsf#338)",
    ),
    Override(
        lambda a: _stirling(a) and max(a) > MAXLGM,
        "xsf has ln Gamma(a) = inf and ln Gamma(a + b) = inf for the larger a > 2.6e305, and "
        "so NaN (scipy/xsf#357)",
    ),
    Override(
        _beta_float_sum_pole,
        "xsf takes a + b as a pole where only its rounded value is one (scipy/xsf#357)",
    ),
)
BETA_FLOORS: Final = (
    Floor(_small_b_beta, 5e-16, SMALL_B),
    Floor(_stirling, 8e-15, STIRLING),
)


def _binom_kx(args: tuple[float, ...]) -> float | None:
    """The `kx` of the integer k case of `binom`, after the symmetry for an integer n > 0."""
    n, k = args
    if not (math.isfinite(k) and k.is_integer() and (abs(n) > 1e-8 or n == 0)):
        return None
    return n - k if n.is_integer() and n > 0 and k > n / 2 else k


def _binom_late(args: tuple[float, ...]) -> bool:
    """Whether `binom` gets past the integer k (kx < 20), `n >= 1e10 k`, `k > 1e8 |n|`, and
    `|k| > 2e4` branches, to the evaluation with `a = 1 + n - k`."""
    n, k = args
    kx = _binom_kx(args)
    return not (
        (kx is not None and kx < 20)
        or (n >= 1e10 * k and k > 0)
        or k > 1e8 * abs(n)
        or (math.isfinite(k) and abs(k) > 2e4 and abs(n) < 170)
    )


def _binom_reflects(args: tuple[float, ...]) -> bool:
    """Whether `binom` reflects Gamma(1 + n) and Gamma(1 + n - k), for a < 0 and a fractional n."""
    n, k = args
    return _binom_late(args) and 1 + n - k < 0 and n % 1 != 0


def _binom_beta(args: tuple[float, ...]) -> tuple[float, float] | None:
    """The arguments of the `beta` at the end of `binom`, if it gets there."""
    n, k = args
    return (
        (1 + n - k, 1 + k) if _binom_late(args) and not _binom_reflects(args) else None
    )


def _in_beta(pred: Pred) -> Pred:
    """`pred` of the arguments of the `beta` at the end of `binom`."""

    def when(args: tuple[float, ...]) -> bool:
        ab = _binom_beta(args)
        return ab is not None and pred(ab)

    return when


def _binom_big_zero(args: tuple[float, ...]) -> bool:
    """Whether xsf's `binom` gets NaN where `kx < 0` after the symmetry: k > 1e8 |n| takes its
    expansion, where Gamma(1 + n) = inf for 1 + n > MAXGAM times sinpi(...) = 0 is NaN."""
    n, k = args
    kx = _binom_kx(args)
    return kx is not None and kx < 0 and k > 1e8 * abs(n) and 1 + n > MAXGAM


def _reflects_pole(args: tuple[float, ...]) -> bool:
    """Whether `binom` reflects, where xsf's `beta(a, 1 + k)` gets a = 1 + n - k that is a pole
    of Gamma only after rounding (1 + k is not a pole)."""
    n, k = args
    exact = Fraction(n) + 1 - Fraction(k)
    return (
        _binom_reflects(args)
        and _nonpos_int(1 + n - k)
        and not _nonpos_int(1 + k)
        and exact.denominator != 1
    )


def _xsf_gamma_mag(x: float) -> float:
    """|Gamma(x)| of xsf for -MAXGAM < x < -33, from the reflection formula and Stirling's."""
    q = -x
    z = q - math.floor(q)
    z = min(z, 1 - z)
    return math.pi / abs(q * math.sin(math.pi * z) * math.gamma(q))


def _reflects_underflow(args: tuple[float, ...]) -> bool:
    """Whether `binom` reflects, where xsf's `beta(a, 1 + k)` evaluates Gamma(a) (|a|, |1 + k|,
    |a + 1 + k| <= MAXGAM) for a = 1 + n - k < -33, which is a subnormal with fewer than 26 bits
    (< 2^-1049), so that its relative error exceeds 2^-26."""
    n, k = args
    a, b = 1 + n - k, 1 + k
    return (
        _binom_reflects(args)
        and -MAXGAM < a < -33
        and abs(a + b) <= MAXGAM
        and abs(b) <= MAXGAM
        and _xsf_gamma_mag(a) < 2.0**-1049
    )


BINOM_OVERRIDES: Final = (
    Override(
        _binom_big_zero,
        "xsf has Gamma(1 + n) = inf for 1 + n > 171.6 in its expansion for k > 1e8 |n|, and so "
        "NaN where the result is 0",
    ),
    Override(
        _reflects_pole,
        "xsf takes 1 + n - k as a pole where only its rounded value is one, and so returns 0 "
        "(scipy/xsf#352)",
    ),
    Override(
        _reflects_underflow,
        "xsf has Gamma(1 + n - k) as a subnormal with fewer than 26 bits, or 0, for "
        "-171.6 < 1 + n - k < -33 (scipy/xsf#352)",
    ),
)
BINOM_FLOORS: Final = (
    Floor(
        _binom_reflects,
        2e-14,
        "xsf-rust reflects Gamma(1 + n) and Gamma(1 + n - k) for 1 + n - k < 0 and a fractional "
        "n (scipy/xsf#352), and so rounds differently",
    ),
    Floor(_in_beta(_stirling), 2e-14, STIRLING),
    Floor(_in_beta(_small_b_beta), 3e-16, SMALL_B),
)


def _xsf_poch_loops(a: float, m: float) -> tuple[float, float] | None:
    """The `r` and the reduced `m` of the loops of xsf's `poch`, or None if they take > 100 steps
    (so that r is still finite and nonzero)."""
    r = 1.0
    for _ in range(100):
        if m >= 1 and a + m != 1:
            m -= 1
            r *= a + m
        elif m <= -1 and a + m != 0:
            r /= a + m
            m += 1
        else:
            return r, m
        if not math.isfinite(r) or r == 0:
            return r, m
    return None


def _xsf_pole(x: float) -> bool:
    """`is_nonpos_int` of xsf's `poch`, which is false for |x| >= 1e13."""
    return _nonpos_int(x) and abs(x) < 1e13


def _poch_final_nan(args: tuple[float, ...]) -> str | None:
    """How the last line of xsf's `poch`, `r * exp(lgam(a + m) - lgam(a)) * gammasgn(a + m) *
    gammasgn(a)` with the reduced m of its loops, is NaN: "pole" for `gammasgn` at a pole (a or
    a + m is a negative integer, both below 1e13 or one of them not), "lgam" for r = 0 or inf
    times exp(+-inf) or 0, or inf - inf (lgam is inf for x > 2.6e305), else None."""
    a, m = args
    if not all(map(math.isfinite, args)) or (loops := _xsf_poch_loops(a, m)) is None:
        return None
    r, m = loops
    s = a + m
    if m == 0 or (a > 1e4 and abs(m) <= 1):
        return None
    if _xsf_pole(s) and not _xsf_pole(a) and s != m:
        return None
    if _xsf_pole(a) and not _xsf_pole(s):
        return None
    if (a < 0 and a.is_integer()) or (s < 0 and s.is_integer()):
        return "pole"
    inf_s, inf_a = s > MAXLGM, a > MAXLGM
    if (inf_s and inf_a) or (r == 0 and inf_s) or (math.isinf(r) and inf_a):
        return "lgam"
    return None


def _poch_float_sum_pole(args: tuple[float, ...]) -> bool:
    """Whether xsf's `poch` returns inf from a rounded a + m' that is a pole (a is not) for the
    reduced m' of its loops, but the exact a + m is not a pole."""
    a, m = args
    if not all(map(math.isfinite, args)) or (loops := _xsf_poch_loops(a, m)) is None:
        return False
    s = a + loops[1]
    return (
        loops[1] != 0
        and not (a > 1e4 and abs(loops[1]) <= 1)
        and _xsf_pole(s)
        and not _xsf_pole(a)
        and s != loops[1]
        and not _sum_nonpos_int(args)
    )


POCH_OVERRIDES: Final = (
    Override(
        lambda a: _poch_final_nan(a) == "pole",
        "xsf has gammasgn(x) = NaN for a negative integer x = a or a + m (after the reduction "
        "of m), where the exact values need not be poles (scipy/xsf#357)",
    ),
    Override(
        lambda a: _poch_final_nan(a) == "lgam",
        "xsf has NaN from inf - inf in lgam(a + m) - lgam(a), or from r = 0 or inf (its loops "
        "stopped on an under- or overflow) times exp(+-inf) (scipy/xsf#357)",
    ),
    Override(
        _poch_float_sum_pole,
        "xsf takes a + m as a pole where only its rounded value is one (scipy/xsf#357)",
    ),
)


def _lanczos[T: (mp.mpf, flint.arb)](x: T) -> T:
    """The rational function of `lanczos_sum_expg_scaled`, exactly."""

    def horner(coef: tuple[float, ...]) -> T:
        acc = x * 0 + coef[0]
        for c in coef[1:]:
            acc = acc * x + c
        return acc

    return horner(LANCZOS_NUM) / horner(LANCZOS_DENOM)


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


def _exp_square[**P](f: Callable[P, flint.arb]) -> Callable[P, flint.arb]:
    """`f` with `EXP_SQUARE_BITS` more bits."""

    def arb(*args: P.args, **kwargs: P.kwargs) -> flint.arb:
        with flint.ctx.extraprec(EXP_SQUARE_BITS):
            return f(*args, **kwargs)

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


def _erfcx_arb[B: (flint.arb, flint.acb)](x: B) -> B:
    return x.erfc() * (x * x).exp()


def _dawsn_arb[B: (flint.arb, flint.acb)](x: B) -> B:
    return flint.arb.pi().sqrt() / 2 * (-x * x).exp() * x.erfi()


def _zero_signs(
    when: Callable[[tuple[float, ...], float], bool],
    value: Callable[[float], complex],
    why: str,
) -> tuple[Rule, Rule]:
    """The `Rule`s for both signs s = +-1 of a zero part of the result."""
    return (
        Rule(lambda a: when(a, 1.0), value(1.0), why),
        Rule(lambda a: when(a, -1.0), value(-1.0), why),
    )


def _wofz_arb(z: flint.acb) -> flint.acb:
    """w(z) = e^(-z^2) erfc(-i z)."""
    return (-z * z).exp() * (-flint.acb(0, 1) * z).erfc()


# w(x) = e^(-x^2) + 2i D(x) / sqrt(pi) for a real x, with Dawson's D(x) ~ 1 / (2x)
WOFZ_AT_INF: Final = (
    Rule(
        lambda a: a[0] == math.inf and a[1] == 0,
        0j,
        "w(z) -> 0 + 0i for z -> +inf along Im z = +-0",
    ),
    Rule(
        lambda a: a[0] == -math.inf and a[1] == 0,
        complex(0, -0.0),
        "w(z) -> 0 - 0i for z -> -inf along Im z = +-0",
    ),
    *_zero_signs(
        lambda a, s: a[0] == 0 and math.copysign(1, a[0]) == s and math.isnan(a[1]),
        lambda s: complex(math.nan, s * 0.0),
        "w(z) = NaN + 0i for z = 0 + NaN i, with the sign of Re z, as w(iy) = erfcx(y) is real",
    ),
)


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
    return Out(
        f"kelvin.{j}",
        tuple(_kelvin(k) for k in parts),
        tuple(_kelvin_arb(k) for k in parts),
        (Rule(lambda a: a[0] < 0, COMPLEX_NAN.value, KE_NEGATIVE),) if j % 2 else (),
        (KELVIN_OVERRIDE,),
    )


MODIFIED_FRESNEL_OVERRIDE: Final = Override(
    lambda a: a[0] < 0,
    "xsf subtracts the updated real part from the imaginary part of F for x < 0, where "
    "xsf-rust differs (scipy/xsf#367)",
)


def _modified_fresnel(k: int, part: int) -> Callable[[mp.mpf], mp.mpf]:
    """Part (0 re, 1 im) of F+, K+, F-, K- (k = 0 to 3)."""

    def ref(x: mp.mpf) -> mp.mpf:
        if x == 0 and k % 2 and part:
            return mp.mpf(0)  # K(0) = 1/2, where the phases cancel
        z = specs.c_modified_fresnel(x)[k]
        return z.imag if part else z.real

    return ref


def _modified_fresnel_arb(k: int, part: int) -> Callable[[flint.arb], flint.arb | None]:
    """F+(x) = sqrt(pi/2) (1/2 - C(u) + i (1/2 - S(u))) with u = x sqrt(2/pi), K+(x) =
    exp(-i (x^2 + pi/4)) F+(x) / sqrt(pi), and the conjugates for F- and K-."""

    def arb(x: flint.arb) -> flint.arb | None:
        if x == 0 and k % 2 and part:
            return None  # the ball of the exact 0 cannot be told from 0
        u = x * (2 / flint.arb.pi()).sqrt()
        f = (flint.arb.pi() / 2).sqrt() * flint.acb(
            0.5 - u.fresnel_c(normalized=True), 0.5 - u.fresnel_s(normalized=True)
        )
        if k % 2:
            f *= (
                -flint.acb(0, 1) * (x * x + flint.arb.pi() / 4)
            ).exp() / flint.arb.pi().sqrt()
        if k >= 2:
            f = f.conjugate()
        return f.imag if part else f.real

    return arb


def _modified_fresnel_out(sign: int, i: int) -> Out:
    """F (i = 0) or K (i = 1) of the plus (sign = 0) or minus (sign = 1) functions."""
    k = 2 * sign + i
    return Out(
        f"modified_fresnel_{('plus', 'minus')[sign]}.{i}",
        tuple(_modified_fresnel(k, part) for part in range(2)),
        tuple(_modified_fresnel_arb(k, part) for part in range(2)),
        overrides=(MODIFIED_FRESNEL_OVERRIDE,) if i == 0 else (),
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


AIRYE_NAN: Final = Rule(
    lambda a: a[0] < 0,
    math.nan,
    "the factor e^zeta of Ai and Ai' is not real for x < 0 (zeta = 2/3 x^(3/2) is imaginary)",
)


def _airy_out(i: int, scaled: bool = False) -> Out:
    """Ai, Ai', Bi, or Bi' (i = 0 to 3), like `specs._airy` but without the functions that
    neither it nor its scale needs; with `scaled`, Ai and Ai' times e^zeta, and Bi and Bi' times
    e^-zeta (for x >= 0, and Bi, Bi' unscaled for x < 0, where Re zeta = 0)."""
    d, f = i % 2, (mp.airyai, mp.airybi)[i // 2]
    sign = 1 - i // 2 * 2

    def ref(x: mp.mpf) -> Value:
        with specs.wdps(specs.xd(x) + 10):
            if x >= 0:
                y = f(x, derivative=d)
                return y * mp.exp(sign * 2 * x * mp.sqrt(x) / 3) if scaled else y
            ai, bi = mp.airyai(x, derivative=d), mp.airybi(x, derivative=d)
            return (ai, bi)[i // 2], mp.hypot(ai, bi)

    def arb(x: flint.arb) -> flint.arb:
        y = (x.airy_ai, x.airy_bi)[i // 2](derivative=d)
        return y * (sign * 2 * x * x.sqrt() / 3).exp() if scaled and x >= 0 else y

    if not scaled:
        return Out(f"airy[{i}]", (ref,), (arb,), overrides=AIRY_OVERRIDES)
    rules = (NAN, AIRYE_NAN) if i < 2 else (NAN,)
    return Out(f"airy_scaled[{i}]", (ref,), (arb,), rules)


def _real(
    func: str,
    ref: Callable[[mp.mpf], Value] | Callable[[mp.mpf, mp.mpf], Value],
    arb: RealArb,
    rules: tuple[Rule, ...] = (),
    overrides: tuple[Override, ...] = (),
) -> Entry:
    """A table with a real output."""
    return Entry((Out(func, (ref,), (arb,), rules, overrides),))


def _complex_out(
    func: str,
    ref: Callable[..., mp.mpc],
    arb: Callable[[flint.acb], flint.acb] | Callable[[flint.acb, flint.acb], flint.acb],
    rules: tuple[Rule, ...] = (),
    overrides: tuple[Override, ...] = (),
) -> Out:
    """A complex output of complex inputs, which are (real, imaginary) pairs of f64s: `ref`
    takes the parts as `mpf`s and `arb` the inputs as `acb`s."""

    def ref_part(part: int) -> Callable[..., mp.mpf]:
        def value(*a: mp.mpf) -> mp.mpf:
            return (mp.re, mp.im)[part](ref(*a))

        return value

    def arb_part(part: int) -> Arb:
        def ball(*a: flint.arb) -> flint.arb:
            # the bits that exp(z) and a periodic reduction of z lose to the magnitude of z
            bits = max(math.frexp(float(x))[1] for x in a)
            with flint.ctx.extraprec(max(bits, 0)):
                z = arb(*(flint.acb(*a[k : k + 2]) for k in range(0, len(a), 2)))
            return z.imag if part else z.real

        return ball

    refs, arbs = (ref_part(0), ref_part(1)), (arb_part(0), arb_part(1))
    return Out(func, refs, arbs, rules, overrides)


def _complex(
    func: str,
    ref: Callable[..., mp.mpc],
    arb: Callable[[flint.acb], flint.acb] | Callable[[flint.acb, flint.acb], flint.acb],
    rules: tuple[Rule, ...] = (),
    overrides: tuple[Override, ...] = (),
) -> Entry:
    """A table with the complex output `_complex_out(func, ref, arb, rules, overrides)`."""
    return Entry((_complex_out(func, ref, arb, rules, overrides),))


def _fresnel_arb(i: int) -> Callable[[flint.arb], flint.arb]:
    """S (i = 0) or C (i = 1) with the normalization of scipy."""
    return _exp_square(lambda x: (x.fresnel_s, x.fresnel_c)[i](normalized=True))


FRESNEL_RULES: Final = (
    NAN,
    Rule(lambda a: a[0] == math.inf, 0.5, "the limit of S and C for x -> +inf"),
    Rule(lambda a: a[0] == -math.inf, -0.5, "the limit of S and C for x -> -inf"),
)


OWENS_T_H_INF: Final = (
    Rule(
        lambda a: abs(a[0]) == math.inf and a[1] >= 0,
        0.0,
        "T(h, a) -> 0 for |h| -> inf and a >= 0 (+-0 included)",
    ),
    Rule(
        lambda a: abs(a[0]) == math.inf and a[1] < 0,
        -0.0,
        "T(h, a) -> -0 for |h| -> inf and a < 0 (odd in a)",
    ),
)
# above this H = h^2 / 2 (or a^2 H for a > 1, after the argument reduction), Owen's series
# returns 0 instead of the tiny result, because it needs too many digits for its cancellation
OWENS_T_SERIES_H: Final = 2000
# for a > 1 the series drops T(a h, 1 / a) ~ exp(-a^2 H), which is negligible relative to
# T ~ exp(-H) from (a^2 - 1) H of this on (exp(-100) < 1e-43)
OWENS_T_DROPPED_H: Final = 100


def _owens_t(h: mp.mpf, a: mp.mpf) -> mp.mpf:
    """`specs.r_owens_t`, also for a = +-inf (T = sign(a) Phi(-|h|) / 2), and by quadrature
    where its series would return 0 although the result is not negligible."""
    if abs(a) == mp.inf:
        return mp.sign(a) * mp.ncdf(-abs(h)) / 2
    big_h = h * h / 2
    dropped = big_h * (a * a - 1)
    if big_h * max(a * a, 1) <= OWENS_T_SERIES_H or dropped >= OWENS_T_DROPPED_H:
        return specs.r_owens_t(h, a)
    # T = exp(-H) int_0^(a sqrt(H)) exp(-s^2) / (1 + s^2 / H) ds / (2 pi sqrt(H)), where
    # exp(-s^2) is negligible beyond s = 40
    top = min(abs(a) * mp.sqrt(big_h), mp.mpf(40))

    def integrand(s: mp.mpf) -> mp.mpf:
        return mp.exp(-s * s) / (1 + s * s / big_h)

    def scaled(u: mp.mpf) -> mp.mpf:
        return integrand(top * u)

    # (over [0, 1] for top < 1, where the absolute tolerance of `mp.quad` would stop it early)
    f = (
        top * mp.quad(scaled, [0, 1])
        if top < 1
        else mp.quad(integrand, [0, *range(1, int(top))[::2], top])
    )
    return mp.sign(a) * mp.exp(-big_h) * f / (2 * mp.pi * mp.sqrt(big_h))


def _owens_t_arb(h: flint.arb, a: flint.arb) -> flint.arb:
    """T(h, a) = exp(-H) int_0^a exp(-H t^2) / (1 + t^2) dt / (2 pi) for H = h^2 / 2, whose
    integral has no tiny absolute tolerance to miss for a large H."""
    if not a.is_finite():
        return a.sgn() * (abs(h) / flint.arb(2).sqrt()).erfc() / 4
    big_h = h * h / 2

    def integrand(t: flint.acb, _analytic: bool) -> flint.acb:
        return (-big_h * t * t).exp() / (1 + t * t)

    return (
        (-big_h).exp() * flint.acb.integral(integrand, 0, a).real / (2 * flint.arb.pi())
    )


def _log1pmx(x: mp.mpf) -> mp.mpf:
    """`specs.r_log1pmx`, whose series does not end at x = 0."""
    return x if x == 0 else specs.r_log1pmx(x)


# the largest f64 and f32 whose exp is finite: ln of the largest one, rounded down (to the
# multiples of 2^-17, the ulp of the f32s in [64, 128))
LN_MAX: Final = math.log(sys.float_info.max)
LN_MAX_F32: Final = math.floor(math.log((2 - 2.0**-23) * 2.0**127) * 2**17) / 2**17


def _logistic(sig: str, suffix: str, ln_max: float) -> dict[tuple[str, str], Entry]:
    """The `expit`, `log_expit`, and `logit` tables of a float type, whose exp(-x) overflows
    for x < -`ln_max`."""
    return {
        ("expit", sig): _real(
            f"expit{suffix}",
            specs.r_expit,
            lambda x: 1 / (1 + (-x).exp()),
            overrides=(
                Override(
                    lambda a: a[0] < -ln_max,
                    f"xsf returns 0 where exp(-x) overflows, for x < -{ln_max}",
                ),
            ),
        ),
        ("log_expit", sig): _real(
            f"log_expit{suffix}", specs.r_log_expit, lambda x: -(-x).exp().log1p()
        ),
        ("logit", sig): _real(
            f"logit{suffix}",
            specs.r_logit,
            lambda p: p.log() - (-p).log1p(),
            (
                Rule(
                    lambda a: not 0 <= a[0] <= 1,
                    math.nan,
                    "expit takes values in [0, 1]",
                ),
                Rule(lambda a: a[0] == 0, -math.inf, "expit(x) -> 0 for x -> -inf"),
                Rule(lambda a: a[0] == 1, math.inf, "expit(x) -> 1 for x -> +inf"),
            ),
        ),
    }


# of the last input
LOG1P_DOMAIN: Final = Rule(lambda a: a[-1] < -1, math.nan, "log1p is not real below -1")
# like in xsf, also where the logarithm is infinite
X_ZERO: Final = Rule(lambda a: a[0] == 0, 0.0, "x log y = 0 for x = 0")
COMPLEX_X_ZERO: Final = Rule(lambda a: a[0] == a[1] == 0, 0j, X_ZERO.why)


def _log1p_at_inf(phase: float, name: str) -> Rule:
    """The limit of ln(1 + z) for z -> inf with arg z = `phase`, like C's clog."""
    return Rule(
        lambda a: (
            any(map(math.isinf, a))
            and not any(map(math.isnan, a))
            and math.atan2(a[1], a[0]).hex() == phase.hex()
        ),
        complex(math.inf, phase),
        f"ln(1 + z) -> inf + {name} i for an infinite z with arg z = {name}, like C's clog",
    )


LOG1P_AT_INF: Final = (
    *(
        _log1p_at_inf(k * math.pi / 4, name)
        for k, name in enumerate(("0", "pi/4", "pi/2", "3 pi/4", "pi"))
    ),
    Rule(
        lambda a: any(map(math.isinf, a)) and any(map(math.isnan, a)),
        complex(math.inf, math.nan),
        "ln(1 + z) -> inf + NaN i for an infinite z with a NaN part, like C's clog",
    ),
)


def _expm1(re: mp.mpf, im: mp.mpf) -> mp.mpc:
    """e^z - 1 = expm1(x) cos y - 2 sin^2(y/2) + i e^x sin y, since the complex `mp.expm1`
    loses the real part where e^x cos y - 1 cancels."""
    with specs.wdps(specs.xd(re, im)):
        real = mp.expm1(re) * mp.cos(im) - 2 * mp.sin(im / 2) ** 2
        return mp.mpc(real, mp.exp(re) * mp.sin(im))


def _expm1_at_inf(re: float, im: float) -> Rule:
    """The limit of e^z - 1 = e^x (cos y + i sin y) - 1 for x -> +inf, where cos y and sin y
    have the signs of `re` and `im`."""
    value = complex(re * math.inf, im * math.inf)
    return Rule(
        lambda a: (
            a[0] == math.inf
            and math.isfinite(a[1])
            and a[1] != 0
            and (math.copysign(1, math.cos(a[1])), math.copysign(1, math.sin(a[1])))
            == (re, im)
        ),
        value,
        f"e^z - 1 -> {value} for Re z -> +inf, with a finite y = Im z != 0 where cos y and"
        " sin y have these signs, like C's cexp",
    )


EXPM1_AT_INF: Final = (
    Rule(
        lambda a: a[0] == math.inf and not math.isfinite(a[1]),
        complex(math.inf, math.nan),
        "e^z - 1 -> inf + NaN i for Re z -> +inf, with an infinite or NaN Im z, like C's cexp",
    ),
    Rule(
        lambda a: a[0] == -math.inf and not math.isfinite(a[1]),
        complex(-1, 0),
        "e^z - 1 -> -1 + 0i for Re z -> -inf, with an infinite or NaN Im z, like C's cexp",
    ),
    *_zero_signs(
        lambda a, s: a[0] == -math.inf and math.copysign(1, math.sin(a[1])) == s,
        lambda s: complex(-1, s * 0.0),
        "e^z - 1 -> -1 + 0i for Re z -> -inf, with a finite Im z and the sign of sin(Im z),"
        " like C's cexp",
    ),
    *_zero_signs(
        lambda a, s: math.isnan(a[0]) and a[1] == 0 and math.copysign(1, a[1]) == s,
        lambda s: complex(math.nan, s * 0.0),
        "e^z - 1 = NaN + 0i for z = NaN + 0i, with the sign of Im z, like C's cexp",
    ),
    COMPLEX_NAN,
    Rule(
        lambda a: math.isfinite(a[0]) and math.isinf(a[1]),
        COMPLEX_NAN.value,
        "e^z - 1 has no limit for Im z -> +-inf, with a finite Re z",
    ),
    *_zero_signs(
        lambda a, s: a[0] == math.inf and a[1] == 0 and math.copysign(1, a[1]) == s,
        lambda s: complex(math.inf, s * 0.0),
        "e^z - 1 -> inf + 0i for z -> +inf along Im z = 0, with its sign, like C's cexp",
    ),
    *(_expm1_at_inf(re, im) for re in (1, -1) for im in (1, -1)),
)


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


def _cut(f: str) -> Rule:
    """NaN for x < 0, where `f` has a branch cut."""
    return Rule(
        lambda a: a[0] < 0,
        math.nan,
        f"{f} has a branch cut along the negative real axis",
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
            entry("y", at_n(specs.r_yv), (_cut(f"Y_{n}"),)),
            entry("i", at_n(specs.r_iv), odd),
            entry("ie", at_n(specs.r_ie), odd, ie),
            entry("k", at_n(specs.r_kv), (_cut(f"K_{n}"),)),
            entry("ke", at_n(specs.r_ke), (_cut(f"K_{n}"),), ke),
        ]
    )


def _j_pole(v: float) -> int:
    """The sign of the limit of J_v(x) for x -> 0+, 0 if it is finite: J_v(x) ~ (x / 2)^v /
    Gamma(v + 1), where Gamma(v + 1) has the sign (-1)^floor(-v) for a non-integer v < 0."""
    return 0 if v >= 0 or v.is_integer() else -1 if math.floor(-v) % 2 else 1


def _y_pole(v: float) -> int:
    """The sign of the limit of Y_v(x) for x -> 0+, 0 if it is finite: Y_v(x) ~ -Gamma(v) / pi
    (2 / x)^v for v > 0, and cos(pi v) times that of -v for v < 0, which is 0 for a
    half-integer v."""
    if v >= 0:
        return -1
    f = Fraction(v) % 2
    return 0 if f.denominator == 2 else -1 if f < 0.5 or f > 1.5 else 1


def _poles_at_0(name: str, pole: Callable[[float], int]) -> tuple[Rule, ...]:
    return tuple(
        Rule(
            lambda a, s=s: a[1] == 0 and pole(a[0]) == s,
            math.copysign(math.inf, s),
            f"{name} has the limit {s * math.inf} for x -> 0+, at x = +-0",
        )
        for s in (1, -1)
    )


# the orders above which all f64 are even integers
EVEN_ORDERS: Final = 2.0**53


# the terms and bits that the series of mpmath need for x up to 1e8, above its defaults
SERIES: Final = {"maxterms": 10**8, "maxprec": 20000}


def _yv(v: mp.mpf, x: mp.mpf) -> Value:
    """`specs.r_yv` with `SERIES` (also for its Bessel modulus), or for |v| > `EVEN_ORDERS` and
    x^2 < |v|, where mpmath does not converge for large x, Y_|v|(x) (DLMF 10.4.1) by the finite
    sum of DLMF 10.8.1, as the rest is smaller by a factor of about (x / 2)^(2 |v|) /
    Gamma(|v|)^2."""
    n = abs(v)
    with specs.wdps(specs.xd(v, x)):
        if n <= EVEN_ORDERS or x * x >= n:
            y = mp.bessely(v, x, **SERIES)
            return (y, mp.hypot(mp.besselj(v, x, **SERIES), y)) if x > n else y
        s, t, k = mp.mpf(0), mp.mpf(1), 0
        while t > mp.eps * s:
            s, k = s + t, k + 1
            t *= x * x / 4 / (k * (n - k))
        return -mp.gamma(n) / mp.pi * (2 / x) ** n * s


def _iv(v: mp.mpf, x: mp.mpf) -> mp.mpf:
    """`specs.r_iv` with `SERIES`."""
    with specs.wdps(specs.xd(v, x)):
        return mp.besseli(v, x, **SERIES)


# the orders above which AMOS returns NaN
AMOS_ORDERS: Final = 2.0**51
LOG_MAX: Final = math.log(sys.float_info.max)
CEPHES_MACHEP: Final = 2.0**-53


def _inf(f: Callable[[], float]) -> bool:
    """Whether `f()` overflows, like in C."""
    try:
        return math.isinf(f())
    except OverflowError:
        return True


def _inputs(vs: tuple[float, ...], xs: tuple[float, ...]) -> Pred:
    """Whether (v, x) is one of `vs` x `xs`, the inputs where xsf fails, checked one by one."""
    rows = {(v, x) for v in vs for x in xs}
    return lambda a: a in rows


def _jv_nan(a: tuple[float, ...]) -> bool:
    """Whether cephes `jv` is NaN for x > 0, as pow(x / 2, n) / Gamma(n + 1) of its power
    series for x^2 < |n + 1| eps / 2 is inf * 0 (with n = |v| for an integer v)."""
    v, x = a
    n = abs(v) if v.is_integer() else v
    return (
        x > 0
        and x * x < abs(n + 1) * CEPHES_MACHEP
        and _inf(lambda: math.pow(x / 2, n))
    )


def _asymptotic_nan(a: tuple[float, ...]) -> bool:
    """Whether I_v(x) of xsf is NaN instead of inf, as 2 pi x of its large-x expansion (for
    |v| <= 50) or 1 + (x / v)^2 of its uniform asymptotic expansion (for |v| > 50) overflows."""
    v, x = a
    n = abs(v)
    y = (abs(x) / n) * (abs(x) / n) if n > 50 else 2 * math.pi * abs(x)
    return not (x < 0 and not v.is_integer()) and math.isinf(y)


def _uniform_overflow(a: tuple[float, ...]) -> bool:
    """Whether e^(v eta) (or for a non-integer v < 0, e^(-v eta) of its K_v term) overflows in
    the uniform asymptotic expansion of I_v(x) of xsf for |v| > 50, but I_v(x) does not."""
    v, x = a
    n, k_term = abs(v), v < 0 and not v.is_integer()
    if n <= 50 or x == 0 or x < 0 and not v.is_integer() or _asymptotic_nan(a):
        return False
    z = abs(x) / n
    # like `ikv_asymptotic_uniform` of xsf
    t = 1 / math.sqrt(1 + z * z)
    e = n * (math.sqrt(1 + z * z) + math.log(z / (1 + 1 / t)))
    # I_v(x) ~ e^e sqrt(t / (2 pi n)), which is inf beyond this bound of e
    i_overflows = (
        _inf(lambda: math.exp(e)) and e < LOG_MAX + 1 - math.log(t / (math.tau * n)) / 2
    )
    k_overflows = k_term and _inf(lambda: math.exp(-e))
    return (i_overflows or k_overflows) and math.isfinite(float(_iv(*map(mp.mpf, a))))


J_OVERRIDES: Final = (
    Override(
        lambda a: a[1] < 0 and a[0].is_integer() and not -(2**31) <= a[0] < 2**31,
        "xsf returns NaN for x < 0 and an integer v outside of the range of int, as its integer"
        " test casts v to int",
    ),
    Override(
        lambda a: abs(a[0]) > AMOS_ORDERS and _jv_nan(a),
        "xsf returns NaN for |v| > 2^51 (the limit of AMOS) where its cephes fallback is NaN",
    ),
    Override(
        _inputs((10.0, -10.0), (1e-30, -1e-30)),
        "xsf returns 0 where AMOS underflows, which J_10(1e-30) = 2.7e-310 does"
        " (scipy/xsf#402)",
    ),
)
Y_OVERRIDES: Final = (
    Override(
        lambda a: abs(a[0]) > AMOS_ORDERS and a[0].is_integer() and a[1] > 0,
        "xsf returns NaN for an integer |v| > 2^51 (the limit of AMOS) and x > 0, as its cephes"
        " fallback casts v to int",
    ),
    Override(
        _inputs((10.0, -10.0), (1e-30,)),
        "xsf returns -inf where AMOS overflows, which Y_10(1e-30) = -1.2e308 does"
        " (scipy/xsf#372)",
    ),
)
I_OVERRIDES: Final = (
    Override(
        _uniform_overflow,
        "xsf returns inf for |v| > 50 where e^(+-v eta) of its uniform asymptotic expansion"
        " overflows, before I_v(x) does",
    ),
    Override(
        _asymptotic_nan,
        "xsf returns NaN instead of inf where its asymptotic expansions overflow in an"
        " intermediate",
    ),
)
# AMOS's `rl` and `fnul`, computed as in `MACH` of `amos.rs`
AMOS_DIG: Final = min(52 * math.log10(2), 18)
AMOS_RL: Final = 1.2 * AMOS_DIG + 3
AMOS_FNUL: Final = 10 + 6 * (AMOS_DIG - 3)


def _amos_not_j(a: tuple[float, ...]) -> bool:
    """Whether J_v(x) isn't computed by `besj`: NaN for a non-integer v and x < 0, `besh` for
    `v <= fnul < x` (`amos_j`), or `ierr = 4` above `0.5 / eps`."""
    v, x = a
    nan = x < 0 and not v.is_integer()
    return nan or abs(v) <= AMOS_FNUL < x or max(abs(v), abs(x)) > 2**51


# `ZUNHJ_GAMA` of `amos.rs`: the power series of `zeta / w2` in `w2 = 1 - (x / n)^2` of `unhj`
AMOS_GAMA: Final = (
    6.29960524947436582e-01, 2.51984209978974633e-01, 1.54790300415655846e-01,
    1.10713062416159013e-01, 8.57309395527394825e-02, 6.97161316958684292e-02,
    5.86085671893713576e-02, 5.04698873536310685e-02, 4.42600580689154809e-02,
    3.93720661543509966e-02, 3.54283195924455368e-02, 3.21818857502098231e-02,
    2.94646240791157679e-02, 2.71581677112934479e-02, 2.51768272973861779e-02,
    2.34570755306078891e-02, 2.19508390134907203e-02, 2.06210828235646240e-02,
    1.94388240897880846e-02, 1.83810633800683158e-02, 1.74293213231963172e-02,
    1.65685837786612353e-02, 1.57865285987918445e-02, 1.50729501494095594e-02,
    1.44193250839954639e-02, 1.38184805735341786e-02, 1.32643378994276568e-02,
    1.27517121970498651e-02, 1.22761545318762767e-02, 1.18338262398482403e-02,
)  # fmt: skip


def _amos_airy_arg(x: float, n: float) -> tuple[float, float]:
    """The Airy argument `n^(2/3) zeta` of `unhj` for J_n(x) with x > 0, rounded as in AMOS (its
    `C64` operations reduce to these, as the imaginary parts are 0 or negligible)."""
    zb = x * (1 / n)
    w2 = 1 - zb * zb
    if abs(w2) <= 0.25:
        p, suma, ap = 1.0, AMOS_GAMA[0], 1.0
        if abs(w2) >= 2**-52:
            for g in AMOS_GAMA[1:]:
                p *= w2
                suma += p * g
                ap *= abs(w2)
                if ap < 2**-52:
                    break
        zeta = (w2 * suma, 0.0)
    elif w2 > 0:
        return (1.0, 0.0)  # positive, so `airy` doesn't call `acai`
    else:
        # `w = t i`, `za = (1 + w) / zb`, and `zth = 3/2 (ln(za) - w)`, whose real part (from
        # ln|za| ~ eps, clamped to >= 0) doesn't change `|zth| = -zth.im`, so its angle is 3 pi/2
        t = math.sqrt(-w2)
        zci = min(max(math.atan2(t / zb, 1 / zb), 0.0), 1.57079632679489662)
        pp, ang = ((t - zci) * 1.5) ** (2 / 3), 4.71238898038468986 * (2 / 3)
        zeta = (pp * math.cos(ang), max(pp * math.sin(ang), 0.0))
    fn13 = n ** (1 / 3)
    return (zeta[0] * (fn13 * fn13), zeta[1] * (fn13 * fn13))


def _amos_acai_mlri(arg: tuple[float, float]) -> bool:
    """Whether `airy` of `arg` calls `acai`, and that `mlri`, for Ai (order 1/3) or Ai' (2/3)."""
    zr, zi = arg
    if zr > 0 or math.hypot(zr, zi) <= 1:
        return False
    # `zta = 2/3 z sqrt(z)`, with `csqrt(z)` for `zr < 0`, and `|zta| = |zta.im|` where `zi = 0`
    t = math.sqrt((-zr + math.hypot(zr, zi)) * 0.5)
    rx, ry = abs(zi) / (2 * t), math.copysign(t, zi)
    zta = ((zr * rx - zi * ry) * (2 / 3), (zr * ry + zi * rx) * (2 / 3))
    az = abs(zta[1]) if zi == 0 else math.hypot(zta[0], zta[1])
    return any(2 < az < AMOS_RL and az * az * 0.25 > f + 1 for f in (1 / 3, 2 / 3))


def _amos_airy_mlri(a: tuple[float, ...]) -> bool:
    """Whether `besj` computes J_v(x) with Airy functions (`binu` -> `buni` -> `uni2`), whose
    `acai` uses `mlri`."""
    v, x = abs(a[0]), abs(a[1])
    # `binu` uses `buni` for `max(v, x) > fnul`, unless `asyi` for `2 x >= v^2`
    if _amos_not_j(a) or max(v, x) <= AMOS_FNUL or 2 * x >= v * v:
        return False
    # `buni` evaluates `uni2` at the orders `v + nui` and `v + nui + 1` if `nui > 0`
    nui = max(int(AMOS_FNUL - v) + 1, 0)
    orders = (v + nui, v + nui + 1) if nui else (v,)
    return any(_amos_acai_mlri(_amos_airy_arg(x, n)) for n in orders)


def _amos_binu_mlri(a: tuple[float, ...]) -> bool:
    """Whether `besj` computes J_v(x) with `mlri` in `binu`: outside of its power series
    (`seri`), for `x < rl` and v <= fnul (else `buni`), or at `x = rl` unless `asyi`."""
    v, x = abs(a[0]), abs(a[1])
    if _amos_not_j(a) or x <= 2 or x * x / 4 <= v + 1 or v > AMOS_FNUL:
        return False
    return x < AMOS_RL or x == AMOS_RL and 2 * x < v * v


J_AIRY_FLOOR: Final = Floor(
    _amos_airy_mlri,
    2e-15,
    "AMOS's Miller algorithm (`mlri`) for Airy functions is more accurate than in xsf, but"
    " rounds differently (7.1e-15 at (442, 481.23), where xsf has 1.8e-15)",
)
J_MILLER_FLOOR: Final = Floor(
    _amos_binu_mlri,
    5e-16,
    "AMOS's Miller algorithm (`mlri`) is more accurate than in xsf, but its rounding depends"
    " on the platform's libm (1.05e-15 at (6.67, 6.67) on macOS and Windows)",
)
I_FLOOR: Final = Floor(
    lambda a: a[0] < 0 and not a[0].is_integer(),
    3e-16,
    "for a non-integer v < 0, the error of I_|v| + 2 / pi sin(pi |v|) K_|v| depends on the"
    " platform's libm (1.15e-15 at (-48.6, 0.0551) on macOS and Windows)",
)
Y_FLOOR: Final = Floor(
    lambda a: a[0] < 0 and not a[0].is_integer(),
    4e-16,
    "for a non-integer v < 0, the error of cos(pi |v|) Y_|v| + sin(pi |v|) J_|v| depends on"
    " the platform's libm (1.37e-15 at (-7.37, 2.22) on macOS and Windows)",
)
K_OVERRIDE: Final = Override(
    _inputs((-20.0, -10.0, -1.0, -0.5, 0.0, 1.0, 12.0, 12.49), (700.6,)),
    "xsf returns 0 where AMOS underflows, which K_v(700.6) of about 3e-306 does",
)


def _bessel_arb(f: str) -> Arb:
    """Arb's J_v, Y_v, I_v, or K_v (`f`), at -v for an integer v < 0 (DLMF 10.4.1, 10.27.1)
    or for K_v (DLMF 10.27.3), with up to 4 times the precision until its ball can confirm a
    reference, else None."""
    g = _precise(specs.ARB[f"bessel_{f}"])

    def arb(v: flint.arb, x: flint.arb) -> Ball | None:
        n = float(v)
        sign = 1
        if n < 0 and (n.is_integer() or f == "k"):
            v, sign = -v, -1 if f in "jy" and n % 2 else 1
        for k in range(3):
            with flint.ctx.extraprec(flint.ctx.prec * (2**k - 1)):
                if (y := g(v, x)) is not None:
                    return sign * y
        return None

    return arb


def _bessel_v() -> dict[tuple[str, str], Entry]:
    """The tables of J_v, Y_v, I_v, K_v, and their exponentially scaled versions, of real order
    v (J_v and Y_v are not scaled for real x)."""

    def j(v: mp.mpf, x: mp.mpf) -> Value:
        """(-1)^v J_v(|x|) for x < 0 and an integer v, with the Bessel modulus as error scale
        like `specs.r_jv`."""
        with specs.wdps(specs.xd(v, x)):
            j = mp.besselj(v, abs(x)) * (-1 if x < 0 and v % 2 else 1)
            return specs._osc(j, v, abs(x))

    def ie(v: flint.arb, x: flint.arb) -> flint.arb:
        return x.bessel_i(v) * (-abs(x)).exp()

    def ke(v: flint.arb, x: flint.arb) -> flint.arb:
        return x.bessel_k(v, scaled=True)

    def entry(
        f: str,
        ref: Callable[[mp.mpf, mp.mpf], Value],
        arb: Arb,
        *rules: Rule,
        overrides: tuple[Override, ...] = (),
        floors: tuple[Floor, ...] = (),
    ) -> tuple[tuple[str, str], Entry]:
        func = f"bessel_{f}"
        out = Out(func, (ref,), (arb,), rules, overrides, floors)
        return (f"cyl_{func}", "d_d-d"), Entry((out,))

    def not_real(f: str) -> Rule:
        return Rule(
            lambda a: a[1] < 0 and not a[0].is_integer(),
            math.nan,
            f"{f} is not real for x < 0 and a non-integer v",
        )

    def cut(f: str) -> Rule:
        return Rule(
            lambda a: a[1] < 0,
            math.nan,
            f"{f} has a branch cut along the negative real axis",
        )

    k0 = Rule(lambda a: a[1] == 0, math.inf, "the limit for x -> 0+, at x = +-0")
    return dict(
        [
            entry(
                "j",
                j,
                _bessel_arb("j"),
                not_real("J_v"),
                *_poles_at_0("J_v", _j_pole),
                overrides=J_OVERRIDES,
                floors=(J_AIRY_FLOOR, J_MILLER_FLOOR),
            ),
            entry("je", j, specs.ARB["bessel_j"], not_real("J_v")),
            entry(
                "y",
                _yv,
                _bessel_arb("y"),
                cut("Y_v"),
                *_poles_at_0("Y_v", _y_pole),
                overrides=Y_OVERRIDES,
                floors=(Y_FLOOR,),
            ),
            entry("ye", specs.r_yv, specs.ARB["bessel_y"], cut("Y_v")),
            entry(
                "i",
                _iv,
                _bessel_arb("i"),
                not_real("I_v"),
                overrides=I_OVERRIDES,
                floors=(I_FLOOR,),
            ),
            entry("ie", specs.r_ie, ie, not_real("I_v")),
            entry(
                "k",
                specs.r_kv,
                _bessel_arb("k"),
                cut("K_v"),
                k0,
                overrides=(K_OVERRIDE,),
            ),
            entry("ke", specs.r_ke, ke, cut("K_v"), k0),
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

COMPLEX_EI_AT_0: Final = Rule(
    lambda a: a[0] == a[1] == 0,
    complex(-math.inf, -0.0),
    "Ei(z) -> -inf for z -> 0, as -inf - 0i like in xsf",
)
COMPLEX_LOG_AT_0: Final = Rule(
    lambda a: a[0] == a[1] == 0,
    complex(-math.inf, math.nan),
    "-inf + NaN i at z = 0, like in xsf, as the limit of the imaginary part depends on arg z",
)


def _e1_overflow(re_w: Callable[[tuple[float, ...]], float]) -> Override:
    """Where xsf's E1(w) overflows, with Re w = `re_w` of the inputs."""
    return Override(
        lambda a: re_w(a) < -LN_MAX,
        "xsf multiplies the continued fraction of E1(w) by exp(-w) first, which overflows for "
        "Re w < -709.78, and so has NaN parts where xsf-rust has the infinities (#1671)",
    )


def _off_axis(f: Callable[[mp.mpc], mp.mpc]) -> Callable[[mp.mpf, mp.mpf], mp.mpc]:
    """`specs._c(f)` with |log10 |Re z / Im z|| more digits, which a part that vanishes on an
    axis loses near it."""
    ref = specs._c(f)

    def value(re: mp.mpf, im: mp.mpf) -> mp.mpc:
        with specs.wdps(abs(int(mp.log10(abs(re / im)))) if re and im else 0):
            return ref(re, im)

    return value


# Shi and Chi from Ei(z) = -E1(-z) +- pi i and Ei(-z) = -E1(z) +- pi i
SHICHI_OVERFLOW: Final = (_e1_overflow(lambda a: -abs(a[0])),)


def _expn(n: mp.mpf, x: mp.mpf) -> mp.mpf:
    """`specs.r_expn` of the integer part of n, to which xsf and xsf-rust truncate n (like
    scipy's `expn` of a float n), but mpmath's E_n for x > 1e15, where the points
    1 + 1/(x + n) of the quadrature of `specs.r_expn` lose their digits."""
    n = mp.floor(n)
    return mp.expint(int(n), x) if x > 1e15 else specs.r_expn(n, x)


EXPN: Final = _real(
    "expn",
    _expn,
    _exp_square(lambda n, x: x.expint(n.floor())),
    (
        Rule(lambda a: a[0] < 0, math.nan, "E_n is NaN for n < 0, like in xsf"),
        Rule(
            lambda a: a[1] < 0,
            math.nan,
            "NaN for x < 0, where E_n is not real for n >= 1",
        ),
        Rule(
            lambda a: a[0] in (0, 1) and a[1] == 0,
            math.inf,
            "E_n(x) -> +inf for x -> 0+ and n = 0, 1",
        ),
    ),
    (
        Override(
            lambda a: a[1] > 709.78,
            "xsf returns 0 for x > 709.78, where E_n(x) ~ e^-x / x can still be subnormal "
            "(scipy/xsf#536)",
        ),
    ),
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
    overrides: tuple[Override, ...] = (),
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
        overrides,
    )


def _jacobi_arb(u: flint.arb, m: flint.arb) -> tuple[flint.arb, flint.arb, flint.arb]:
    """sn, cn, and dn for m < 1, by theta functions with the nome of tau = i K(1 - m) / K(m),
    and for m < 0 by sn(u|m) = sd(u k'|-m/k'^2) / k', with k' = sqrt(1 - m)."""
    if m < 0:
        k1 = (1 - m).sqrt()
        sn, cn, dn = _jacobi_arb(u * k1, -m / (1 - m))
        return sn / (k1 * dn), cn / dn, 1 / dn
    if m == 0:
        return u.sin(), u.cos(), flint.arb(1)
    k = flint.acb(m).elliptic_k()
    tau = 1j * flint.acb(1 - m).elliptic_k() / k
    t1, t2, t3, t4 = flint.acb.modular_theta(u / (2 * k), tau)
    _, a2, a3, a4 = flint.acb.modular_theta(flint.acb(0), tau)
    return (a3 / a2 * t1 / t4).real, (a4 / a2 * t2 / t4).real, (a4 / a3 * t3 / t4).real


ELLIPJ_RULES: Final = (
    NAN,
    Rule(
        lambda a: a[1] > 1,
        math.nan,
        "sn, cn, dn, and am are real for m > 1, but NaN like in xsf",
    ),
)
ELLIPJ_OVERRIDE: Final = Override(
    lambda a: a[1] == 1 and abs(a[0]) > 710.47,
    "xsf returns NaN for m = 1 where cosh(u) overflows, for |u| > 710.47 (scipy/xsf#375)",
)
ELLIPJ_FLOOR: Final = Floor(
    lambda a: a[1] < 0,
    1e-15,
    "for m < 0, the error depends on the platform's libm (1.15e-15 at (-1, -10) on macOS "
    "and Windows)",
)
AM_OVERRIDE: Final = Override(
    lambda a: a[1] >= 0.9999999999 and abs(a[0]) < 1,
    "xsf's am = 2 atan(e^u) - pi/2 for m >= 0.9999999999 cancels for |u| < 1 (scipy/xsf#375)",
)


def _ellipj_out(i: int) -> Out:
    """sn, cn, dn, or am (i = 0 to 3) of `specs._ellipj`, exact at u = 0, with the limits for
    u -> +-inf with m = 1 (where Arb can't check), and an absolute error for sn = sin am and
    cn = cos am where they oscillate (m < 1)."""

    def ref(u: mp.mpf, m: mp.mpf) -> Value:
        if u == 0:
            f = mp.mpf((0, 1, 1, 0)[i])
        elif mp.isinf(u):
            f = (mp.sign(u), mp.mpf(0), mp.mpf(0), mp.sign(u) * mp.pi / 2)[i]
        else:
            f = specs._ellipj(u, m)[i]
        return (f, mp.mpf(1)) if i < 2 and m < 1 else f

    def arb(u: flint.arb, m: flint.arb) -> flint.arb | None:
        if not u.is_finite():
            return None
        if m == 1:
            sech = 1 / u.cosh()
            return (u.tanh(), sech, sech, u.sinh().atan())[i]
        sn, cn, dn = _jacobi_arb(u, m)
        if i < 3:
            return (sn, cn, dn)[i]
        k = flint.acb(m).elliptic_k().real
        n = math.floor(float((u + k) / (2 * k)))
        s = -1 if n % 2 else 1
        return n * flint.arb.pi() + flint.arb.atan2(s * sn, s * cn)

    overrides = (ELLIPJ_OVERRIDE, AM_OVERRIDE) if i == 3 else (ELLIPJ_OVERRIDE,)
    return Out(f"ellipj.{i}", (ref,), (arb,), ELLIPJ_RULES, overrides, (ELLIPJ_FLOOR,))


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
    f: RealArb,
    ref: Callable[[mp.mpf], mp.mpf]
    | Callable[[mp.mpf, mp.mpf], mp.mpf]
    | Callable[[mp.mpf, mp.mpf, mp.mpf], mp.mpf],
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


def _struve_arb[T: (flint.arb, flint.acb)](v: float, x: T, sign: int) -> T:
    """H_v(x) (sign -1) or L_v(x) (sign 1) by DLMF 11.2.1 or 11.2.2, also for x < 0 with an
    integer v."""
    # in Arb, since v + 1 and v + 3/2 round in f64 for tiny v
    w = flint.arb(v)
    p = int(v) + 1 if v.is_integer() else w + 1
    return (x / 2) ** p * (sign * x * x / 4).hypgeom(
        [1], [1.5, w + 1.5], regularized=True
    )


def _struve(sign: int) -> Entry:
    """The table of H_v (sign -1) or L_v (sign 1)."""
    name = "LH"[sign < 0]

    def ref(v: mp.mpf, x: mp.mpf) -> Value:
        """`specs.r_struve_h` or `specs.r_struve_l`, and for H_v at x < 0 the scale at |x|
        (H_n is even or odd)."""
        if sign > 0:
            return specs.r_struve_l(v, x)
        h = specs.r_struve_h(v, x)
        return h if x >= 0 else specs._osc(h, v, -x)

    def arb(v: flint.arb, x: flint.arb) -> flint.arb:
        # the terms of H_v cancel up to about e^|x|, and for v << 0 Arb needs about |v| / 10
        # more bits (measured), so 2 (|x| + |v|) bits cover both
        with flint.ctx.extraprec(2 * int(abs(float(x)) + abs(float(v)))):
            return _struve_arb(float(v), x, sign)

    def gamma_sign(v: float) -> int:
        """The sign of Gamma(v + 3/2) for v < -1, or 0 at a pole."""
        return 0 if (v + 0.5).is_integer() else -1 if math.ceil(-v - 1.5) % 2 else 1

    def pole(s: int) -> Rule:
        return Rule(
            lambda a: a[1] == 0 and a[0] < -1 and gamma_sign(a[0]) == s,
            s * math.inf,
            f"the limit for x -> 0+, for v < -1 with Gamma(v + 3/2) {'<>'[s > 0]} 0",
        )

    rules = (
        Rule(lambda a: a[1] == 0 and a[0] > -1, 0.0, f"{name}_v(0) = 0 for v > -1"),
        pole(1),
        pole(-1),
    )
    # only the table of H_v has such rows
    if sign < 0:
        why = "H_v(x) is not real for x < 0 and a non-integer v"
        rules += (Rule(lambda a: a[1] < 0 and not a[0].is_integer(), math.nan, why),)
    return _real(f"struve_{name.lower()}", ref, arb, rules)


def _integral_arb(
    f: Callable[[flint.acb], flint.acb],
) -> Callable[[flint.arb], flint.arb]:
    """int_0^x f(t) dt, for an entire f."""
    return lambda x: flint.acb.integral(lambda t, _: f(t), 0, x).real


def _int_c0_arb(c: Callable[[flint.arb, int], flint.arb], sign: int) -> Arb:
    """int_0^x Y_0(t) dt (sign -1) or int_0^x K_0(t) dt (sign 1), like `specs.r_ity0` and
    `specs.r_itk0`, where `c` is Y or K."""

    def arb(x: flint.arb) -> flint.arb:
        c0, c1 = c(x, 0), c(x, 1)
        s0, s1 = (_struve_arb(v, x, sign) for v in (0.0, 1.0))
        return x * c0 + flint.arb.pi() * x / 2 * (c1 * s0 + sign * c0 * s1)

    return arb


INTEGRAL_AT_0: Final = Rule(lambda a: a[0] == 0, 0.0, "int_0^0 f(t) dt = +0")


def _it2k0_residual(x: mp.mpf) -> mp.mpf:
    """`r_it2k0(x) - int_x^inf K_0(t)/t dt` by quadrature."""

    def integrand(t: mp.mpf) -> mp.mpf:
        return mp.besselk(0, t) / t

    return abs(specs.r_it2k0(x) - mp.quad(integrand, [x, mp.inf]))


def _it2i0_arb(x: flint.arb) -> flint.arb:
    """Like `specs.r_it2i0`."""
    return x * x / 8 * (x * x / 4).hypgeom([1, 1], [2, 2, 2])


def _mid_x_override(issue: int) -> Override:
    """Where xsf-rust uses another method than xsf, for 0.5 <= |x| < 40."""
    return Override(
        lambda a: 0.5 <= abs(a[0]) < 40,
        f"xsf loses digits for 0.5 <= |x| < 40, where xsf-rust differs (scipy/xsf#{issue})",
    )


ITAIRY_OVERRIDE: Final = Override(
    lambda a: 2 < abs(a[0]) <= 20,
    "xsf loses digits for 2 < |x| <= 20, where xsf-rust differs (scipy/xsf#448)",
)


def _itairy_out(i: int) -> Out:
    """int_0^x Ai(t), Bi(t), Ai(-t), or Bi(-t) dt (i = 0 to 3)."""

    def ref(x: mp.mpf) -> mp.mpf:
        return specs._airy_int(x)[i]

    def f(t: flint.acb) -> flint.acb:
        return (flint.acb.airy_ai, flint.acb.airy_bi)[i % 2](-t if i // 2 else t)

    arbs = (_integral_arb(f),)
    return Out(f"airy_integrals[{i}]", (ref,), arbs, overrides=(ITAIRY_OVERRIDE,))


# from here on, Temme's expansion (with a relative error of about a^-3, checked against Arb)
# replaces Arb, whose Q(a, x) loses all digits near x = a from about a = 1e17
TEMME_A: Final = 1e15


def _temme(a: mp.mpf, x: mp.mpf, upper: bool) -> mp.mpf:
    """Q(a, x) (or P(a, x)) by Temme's uniform expansion (DLMF 8.12.3, 8.12.4, 8.12.7) to the
    c_2(eta) term, with c_0 and c_1 of DLMF 8.12.9, c_2 by the recursion of DLMF 8.12.10, and
    their limits at eta = 0."""
    mu = (x - a) / a
    # c_2(eta) has a removable singularity of order 5 at mu = 0
    with specs.wdps(5 * specs.xd(1 / mu) if mu else 0):
        mu = (x - a) / a
        c = (mp.mpf(-1) / 3, mp.mpf(-1) / 540, mp.mpf(25) / 6048)
        eta = mp.mpf(0)
        if mu:
            eta = mp.sign(mu) * mp.sqrt(2 * (mu - mp.log1p(mu)))
            c = (
                1 / mu - 1 / eta,
                1 / eta**3 - 1 / mu**3 - 1 / mu**2 - 1 / (12 * mu),
                (1 + mu) * (3 / mu**5 + 2 / mu**4 + 1 / (12 * mu**3))
                - 3 / eta**5
                + 1 / (288 * mu),
            )
        s = c[0] + c[1] / a + c[2] / a**2
        s *= mp.exp(-a * eta**2 / 2) / mp.sqrt(2 * mp.pi * a)
        e = mp.erfc((eta if upper else -eta) * mp.sqrt(a / 2)) / 2
        return e + s if upper else e - s


def _gammainc(upper: bool) -> Callable[[mp.mpf, mp.mpf], mp.mpf]:
    """Q(a, x) (or P(a, x)) of `specs`, or by Temme's expansion for large a."""
    return lambda a, x: (
        _temme(a, x, upper) if a >= TEMME_A else specs._gammainc_pq(a, x, upper)
    )


# the relative accuracy of an Arb ball that can confirm a reference
ARB_BITS: Final = 128


def _gammaincc_arb(a: flint.arb, x: flint.arb) -> flint.arb:
    return x.gamma_upper(a, regularized=1)


def _gammainc_arb(a: flint.arb, x: flint.arb) -> flint.arb:
    """1 - Q(a, x), or where P(a, x) is too small for that, P(a, x) (which is slow for large a)."""
    p = 1 - _gammaincc_arb(a, x)
    if p.rel_accuracy_bits() < ARB_BITS and a < TEMME_A:
        p = x.gamma_lower(a, regularized=1)
    return p


def _precise[*A](f: Callable[[*A], Ball]) -> Callable[[*A], Ball | None]:
    """`f`, but None where its ball is too wide to confirm a reference."""

    def arb(*args: *A) -> Ball | None:
        y = f(*args)
        return y if y.rel_accuracy_bits() >= ARB_BITS else None

    return arb


def _gammainc_rules(upper: bool) -> tuple[Rule, ...]:
    """xsf's conventions for Q(a, x) (or P(a, x)), which is `x_inf` for x -> +inf and `x_0` at
    x = 0."""
    x_inf, x_0 = (0.0, 1.0) if upper else (1.0, 0.0)
    f = "Q" if upper else "P"
    return (
        NAN,
        Rule(lambda a: a[0] < 0 or a[1] < 0, math.nan, "NaN outside of a >= 0, x >= 0"),
        Rule(
            lambda a: math.isinf(a[0]) and math.isinf(a[1]),
            math.nan,
            "no limit for a, x -> +inf",
        ),
        Rule(lambda a: math.isinf(a[0]), x_0, "the limit for a -> +inf, with x finite"),
        Rule(
            lambda a: math.isinf(a[1]), x_inf, "the limit for x -> +inf, with a finite"
        ),
        Rule(lambda a: a[0] == a[1] == 0, math.nan, "no limit for a, x -> 0+"),
        Rule(lambda a: a[0] == 0, x_inf, "the limit for a -> 0+, with x > 0"),
        Rule(lambda a: a[1] == 0, x_0, f"{f}(a, 0) = {x_0:g}, with a > 0"),
    )


def _gammaincinv_rules(upper: bool, at_one: float) -> tuple[Rule, ...]:
    """xsf's conventions for the inverse of Q(a, x) (or P(a, x)), which is `at_one` at 1."""
    f = "Q" if upper else "P"
    return (
        Rule(
            lambda a: a[0] < 0 or not 0 <= a[1] <= 1,
            math.nan,
            f"NaN outside of a >= 0, 0 <= {f.lower()} <= 1",
        ),
        Rule(
            lambda a: a[0] > 0 and a[1] == 1,
            at_one,
            f"{f}(a, {at_one:g}) = 1, with a > 0",
        ),
    )


RGAMMA: Final = (
    "xsf-rust computes 1 / Gamma(a) with rgamma(a) instead of exp(-lgam(a)) for a < 1 "
    "(scipy/xsf#456), and x^a e^-x / Gamma(a) in igam_fac with pow_rgamma instead of "
    "exp(a ln(x) - x - lgam(a))"
)


@functools.cache
def _rust() -> list[str]:
    """The command of xsf-rust's dump program."""
    return [*accuracy.build_rust(), "xsf"]


@functools.cache
def _rgamma(a: float) -> float:
    """rgamma(a) of xsf-rust."""
    return accuracy.run(_rust(), [("rgamma", (a,))])[0].real


def _normal(y: float) -> bool:
    return math.isfinite(y) and abs(y) >= sys.float_info.min


def _igam_fac_rgamma(a: float, x: float) -> bool:
    """Whether `igam_fac(a, x)` uses rgamma(a): for a < 1, or where `pow_rgamma(x, a, 0, a, x)`
    succeeds (the same operations, with the same libm)."""
    if not (abs(a - x) > 0.4 * a and a <= 200):
        return False
    if a < 1:
        return True
    try:
        p = math.pow(x, a / 2) * math.exp(-x / 2)
    except OverflowError:
        return False
    return _normal(p) and _normal(p * _rgamma(a) * 1.0 * p)


def _uses_rgamma(a: float, x: float, upper: bool) -> bool:
    """Whether xsf-rust's Q(a, x) (or P(a, x)) for 0 < a < 1 and x > 0 reaches a method that
    differs from that of xsf (scipy/xsf#456): `igam_fac` or `igamc_series` with rgamma(a), or
    P = 1 - Q where P = 1 - O(a), by the branches of `igamc` (or `igam`)."""
    fac = _igam_fac_rgamma(a, x)
    if not upper and not (x > 1 and x > a):
        return fac or a * max(abs(math.log(x)), 1) <= 0.01
    if x > 1.1:
        return fac
    # otherwise igamc_series, or 1 - igam_series
    return fac or (-0.4 / math.log(x) >= a if x <= 0.5 else x * 1.1 >= a)


def _gammainc_out(upper: bool, floor: float) -> Out:
    """Q(a, x) (or P(a, x)), with the tolerance floor `floor` where `igam_fac` uses rgamma(a),
    so that each branch of `igamc` (or `igam`) does."""
    return Out(
        "gammaincc" if upper else "gammainc",
        (_gammainc(upper),),
        (_precise(_gammaincc_arb if upper else _gammainc_arb),),
        _gammainc_rules(upper),
        floors=(Floor(lambda a: _igam_fac_rgamma(*a), floor, RGAMMA),),
    )


def _gammaincinv_out(upper: bool, floor: float) -> Out:
    """The inverse of Q(a, x) (or P(a, x)), with the tolerance floor `floor` where its Halley
    steps reach rgamma(a) at the root x (by `specs`, with 20 digits): in `igam_fac`, or in Q
    (or P), which they evaluate instead of P (or Q) for y > 0.9."""
    ref = specs.r_gammainccinv if upper else specs.r_gammaincinv

    def when(args: tuple[float, ...]) -> bool:
        a, y = args
        if not 0 < a <= 200:
            return False
        with mp.workdps(20):
            x = float(ref(mp.mpf(a), mp.mpf(y)))
        return (
            _igam_fac_rgamma(a, x) or a < 1 and _uses_rgamma(a, x, upper != (y > 0.9))
        )

    return Out(
        "gammainccinv" if upper else "gammaincinv",
        (ref,),
        (_inverse_arb(_gammaincc_arb if upper else _gammainc_arb, ref),),
        _gammaincinv_rules(upper, 0.0 if upper else math.inf),
        floors=(Floor(when, floor, RGAMMA),),
    )


class ToArgs(Protocol):
    """The (a, x) of Q(a, x), or P(a, x), of the inputs as `float`s, `mpf`s, or balls."""

    def __call__[T: (float, mp.mpf, flint.arb)](
        self, a: tuple[T, ...], /
    ) -> tuple[T, T]: ...


def _rules_at(rules: tuple[Rule, ...], to_args: ToArgs, sub: str) -> tuple[Rule, ...]:
    """The rules of Q(a, x), P(a, x), or their inverses, for the function of other inputs with
    (a, x) (or (a, y)) = `sub`."""
    return tuple(
        Rule(
            lambda a, r=r: r.when(to_args(a)), r.value, f"{r.why}, with (a, x) = {sub}"
        )
        for r in rules
    )


def _floors_at(floors: tuple[Floor, ...], to_args: ToArgs) -> tuple[Floor, ...]:
    return tuple(Floor(lambda a, f=f: f.when(to_args(a)), f.tol, f.why) for f in floors)


def _half[T: (float, mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T]:
    return a[0] / 2, a[1] / 2


def _floor[T: (float, mp.mpf, flint.arb)](x: T) -> T:
    """floor(x), which is x for an infinity or NaN."""
    if isinstance(x, flint.arb):
        return x.floor()
    return type(x)(math.floor(x)) if math.isfinite(x) else x


def _poisson[T: (float, mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T]:
    return _floor(a[0]) + 1, a[1]


def _gamma_dist[T: (float, mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T]:
    return a[1], a[0] * a[2]


def _distribution(
    func: str,
    upper: bool,
    to_args: ToArgs,
    sub: str,
    extra: tuple[Rule, ...] = (),
    floors: tuple[Floor, ...] = (),
) -> Out:
    """A distribution function that is Q(a, x) (or P(a, x)) at (a, x) = `sub`, after the rules
    `extra`."""
    g = _gammainc(upper)
    arb = _precise(_gammaincc_arb if upper else _gammainc_arb)

    def ref(*a: mp.mpf) -> mp.mpf:
        return g(*to_args(a))

    def ball(*a: flint.arb) -> Ball | None:
        return arb(*to_args(a))

    return Out(
        func,
        (ref,),
        (ball,),
        extra + _rules_at(_gammainc_rules(upper), to_args, sub),
        floors=_floors_at(floors, to_args),
    )


X_NEG: Final = Rule(lambda a: a[2] < 0, math.nan, "NaN for x < 0")
K_M_NEG: Final = Rule(
    lambda a: a[0] < 0 or a[1] < 0, math.nan, "NaN outside of k >= 0, m >= 0"
)
PDTRI_DOMAIN: Final = Rule(
    lambda a: a[0] < 0 or not 0 <= a[1] <= 1,
    math.nan,
    "NaN outside of k >= 0, 0 <= y <= 1",
)


def _chdtri_args[T: (float, mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T]:
    return a[0] / 2, a[1]


def _chdtri() -> Out:
    """The x with Q(v / 2, x / 2) = y, which is twice the root of Q(v / 2, .) = y."""
    return Out(
        "chdtri",
        (specs.r_chdtri,),
        (_inverse_arb(lambda v, x: _gammaincc_arb(v / 2, x / 2), specs.r_chdtri),),
        _rules_at(_gammaincinv_rules(True, 0.0), _chdtri_args, "(v / 2, y)"),
    )


def _pdtri() -> Out:
    """The m with Q(k + 1, m) = y, for the integer k >= 0 and 0 <= y <= 1."""
    return Out(
        "pdtri",
        (specs.r_pdtri,),
        (_inverse_arb(lambda k, m: _gammaincc_arb(k + 1, m), specs.r_pdtri),),
        (PDTRI_DOMAIN, Rule(lambda a: a[1] == 1, 0.0, "pdtri(k, 1) = 0")),
        (
            Override(
                lambda a: a[1] == 1, "xsf returns NaN for y = 1, inside the domain"
            ),
        ),
    )


def _floor_k(a: tuple[float, ...]) -> float:
    """floor(k) of the first input."""
    return _floor(a[0])


BIN_P: Final = Rule(
    lambda a: not 0 <= a[2] <= 1, math.nan, "NaN outside of 0 <= p <= 1"
)
BIN_N: Final = Rule(lambda a: a[1] < _floor_k(a), math.nan, "NaN for n < floor(k)")
BDTR_RULES: Final = (
    BIN_P,
    BIN_N,
    Rule(lambda a: a[0] < 0, math.nan, "NaN for k < 0"),
    Rule(lambda a: a[1] == _floor_k(a), 1.0, "bdtr(k, n, p) = 1 for floor(k) = n"),
)
BDTRC_RULES: Final = (
    BIN_P,
    BIN_N,
    Rule(lambda a: a[0] < 0, 1.0, "bdtrc(k, n, p) = 1 for k < 0"),
    Rule(lambda a: a[1] == _floor_k(a), 0.0, "bdtrc(k, n, p) = 0 for floor(k) = n"),
)
BDTRI_RULES: Final = (
    Rule(
        lambda a: not 0 <= a[2] <= 1 or a[0] < 0 or a[1] <= _floor_k(a),
        math.nan,
        "NaN outside of 0 <= y <= 1, k >= 0, n > floor(k)",
    ),
    Rule(lambda a: a[2] == 1, 0.0, "bdtri(k, n, 1) = 0, since bdtr(k, n, 0) = 1"),
)
NBDTR_RULES: Final = (
    Rule(
        lambda a: a[0] < 0 or not 0 <= a[2] <= 1,
        math.nan,
        "NaN outside of k >= 0, 0 <= p <= 1",
    ),
    Rule(lambda a: a[1] <= 0, math.nan, "NaN for n <= 0"),
)
FDTR_DOMAIN: Final = Rule(
    lambda a: a[0] <= 0 or a[1] <= 0 or a[2] < 0,
    math.nan,
    "NaN outside of a, b > 0, x >= 0",
)
FDTRI_RULES: Final = (
    Rule(
        lambda a: a[0] <= 0 or a[1] <= 0 or not 0 <= a[2] <= 1,
        math.nan,
        "NaN outside of a, b > 0, 0 <= y <= 1",
    ),
    Rule(lambda a: a[2] == 1, math.inf, "the limit for y -> 1"),
)


def _betainc_arb(a: flint.arb, b: flint.arb, x: flint.arb, y: flint.arb) -> flint.arb:
    """I_x(a, b) with y = 1 - x, as I_x(a, b) or 1 - I_y(b, a), by the one that is more
    precise (starting on the side of the mean of the beta distribution, as `specs._betainc`)."""
    low = x < (a + 1) / (a + b + 2)
    r = (
        x.beta_lower(a, b, regularized=1)
        if low
        else 1 - y.beta_lower(b, a, regularized=1)
    )
    if r.rel_accuracy_bits() >= ARB_BITS:
        return r
    return (
        1 - y.beta_lower(b, a, regularized=1)
        if low
        else x.beta_lower(a, b, regularized=1)
    )


class BetaArgs(Protocol):
    """The (a, b, x, y = 1 - x) of I_x(a, b) of the inputs as `mpf`s, or balls."""

    def __call__[T: (mp.mpf, flint.arb)](
        self, a: tuple[T, ...], /
    ) -> tuple[T, T, T, T]: ...


def _bdtr_args[T: (mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T, T, T]:
    k, n, p = _floor(a[0]), a[1], a[2]
    return n - k, k + 1, 1 - p, p


def _bdtrc_args[T: (mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T, T, T]:
    k, n, p = _floor(a[0]), a[1], a[2]
    return k + 1, n - k, p, 1 - p


def _nbdtr_args[T: (mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T, T, T]:
    k, n, p = a
    return n, k + 1, p, 1 - p


def _nbdtrc_args[T: (mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T, T, T]:
    k, n, p = a
    return k + 1, n, 1 - p, p


def _fdtr_args[T: (mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T, T, T]:
    """w = a x / (b + a x) and 1 - w = b / (b + a x), both without cancellation."""
    d1, d2, x = a
    s = d2 + d1 * x
    return d1 / 2, d2 / 2, d1 * x / s, d2 / s


def _fdtrc_args[T: (mp.mpf, flint.arb)](a: tuple[T, ...], /) -> tuple[T, T, T, T]:
    d1, d2, x = a
    s = d2 + d1 * x
    return d2 / 2, d1 / 2, d2 / s, d1 * x / s


def _beta_dist(func: str, to_args: BetaArgs, rules: tuple[Rule, ...]) -> Out:
    """A distribution function that is I_x(a, b) at (a, b, x) = `to_args`, after `rules`."""
    arb = _precise(_betainc_arb)

    def ref(*a: mp.mpf) -> mp.mpf:
        return specs._betainc(*to_args(a))

    def ball(*a: flint.arb) -> Ball | None:
        return arb(*to_args(a))

    return Out(func, (ref,), (ball,), (NAN, *rules))


def _bdtr_ref(k: mp.mpf, n: mp.mpf, y: mp.mpf) -> mp.mpf:
    """The p with bdtr(floor(k), n, p) = y."""
    return specs.r_bdtri(_floor(k), n, y)


def _bdtri() -> Out:
    return Out(
        "bdtri",
        (_bdtr_ref,),
        (
            _inverse_arb(
                lambda k, n, p: _betainc_arb(*_bdtr_args((k, n, p))), _bdtr_ref
            ),
        ),
        (NAN, *BDTRI_RULES),
    )


def _fdtri() -> Out:
    return Out(
        "fdtri",
        (specs.r_fdtri,),
        (
            _inverse_arb(
                lambda a, b, x: _betainc_arb(*_fdtr_args((a, b, x))), specs.r_fdtri
            ),
        ),
        (NAN, *FDTRI_RULES),
    )


def _besselpoly_order(v: float) -> tuple[float, int]:
    """(|v|, (-1)^v) for a negative integer v, else (v, 1)."""
    if v < 0 and v.is_integer():
        return -v, -1 if -v % 2 else 1
    return v, 1


def _besselpoly_pole(a: tuple[float, ...]) -> int:
    """The sign of the integral at a pole lambda + |nu| + 1 = 0, -2, ..., else 0."""
    x, lmbda, v = a
    if x == 0 or not all(map(math.isfinite, a)) or (x < 0 and not v.is_integer()):
        return 0
    n, sign = _besselpoly_order(v)
    e = Fraction(lmbda) + Fraction(n) + 1
    if e > 0 or e.denominator != 1 or e % 2:
        return 0
    lead = 1 if n + 1 > 0 else (-1) ** math.floor(n + 1)  # the sign of Gamma(n + 1)
    return lead * sign * (-1 if x < 0 and v % 2 else 1)


BESSELPOLY_RULES: Final = tuple(
    Rule(
        lambda a, s=s: _besselpoly_pole(a) == s,
        math.copysign(math.inf, s),
        "the integral diverges at a pole",
    )
    for s in (1, -1)
)

BESSELPOLY_OVERRIDES: Final = (
    Override(
        lambda a: a[0] < 0 and a[2].is_integer(),
        "xsf returns NaN for a < 0",
    ),
    Override(
        lambda a: a[0] > 0 and _besselpoly_pole(a) != 0,
        "xsf returns a partial sum at a pole",
    ),
)


def _besselpoly(a: mp.mpf, lmbda: mp.mpf, v: mp.mpf) -> mp.mpf:
    """`specs.r_besselpoly`, with J_{-n} = (-1)^n J_n."""
    n, sign = _besselpoly_order(float(v))
    return sign * specs.r_besselpoly(a, lmbda, mp.mpf(n))


def _besselpoly_arb(a: flint.arb, lmbda: flint.arb, v: flint.arb) -> flint.arb:
    """`_besselpoly` as a 1F2 in Arb."""
    order, sign = _besselpoly_order(float(v))
    n = flint.arb(order)
    c = (lmbda + n + 1) / 2
    return (
        sign
        * a**n
        / ((lmbda + n + 1) * (n + 1).gamma())
        * (-a * a).hypgeom([c], [n + 1, c + 1])
    )


def _voigt(x: mp.mpf, sigma: mp.mpf, gamma: mp.mpf) -> mp.mpf:
    """`specs.r_voigt_profile`, and its limits for sigma = 0 (Cauchy) and gamma = 0 (normal)."""
    if sigma == 0:
        return gamma / (mp.pi * (x * x + gamma * gamma))
    if gamma == 0:
        return mp.exp(-((x / sigma) ** 2) / 2) / (abs(sigma) * mp.sqrt(2 * mp.pi))
    return specs.r_voigt_profile(x, sigma, gamma)


VOIGT_RULES: Final = (
    NAN,
    Rule(
        lambda a: a[1] == a[2] == 0 and a[0] == 0,
        math.inf,
        "the density of a point mass at x = 0, for sigma = gamma = 0",
    ),
    Rule(
        lambda a: a[1] == a[2] == 0 and a[0] != 0,
        0.0,
        "the density of a point mass at x != 0, for sigma = gamma = 0",
    ),
    Rule(
        lambda a: math.isinf(a[0]) and math.isfinite(a[1]) and math.isfinite(a[2]),
        0.0,
        "the limit for x -> +-inf",
    ),
)


def _wright_bessel_arb(a: flint.arb, b: flint.arb, x: flint.arb) -> flint.arb:
    """Phi(a, b; x) for a, b, x >= 0: exp(x) / Gamma(b) for a = 0, else the power series with
    a ball for the tail.

    The ratio of consecutive terms decreases, since Gamma(y + a) / Gamma(y) increases in y,
    so the tail after a term t with a ratio r <= 1/2 is at most t r / (1 - r) <= t."""
    if x == 0:
        return b.rgamma()
    if a == 0:
        return x.exp() * b.rgamma()
    prec = flint.ctx.prec
    s, power, k = flint.arb(0), flint.arb(1), 0
    while True:
        term = power * (a * k + b).rgamma()
        s += term
        k += 1
        power *= x / k
        if k > 1 and term < s * 2 ** (1 - prec):
            ratio = power * (a * k + b).rgamma() / term
            if ratio < 0.5:
                return s + flint.arb(0, term)


def _log_wright_bessel_arb(a: flint.arb, b: flint.arb, x: flint.arb) -> flint.arb:
    return _wright_bessel_arb(a, b, x).log()


def _log_wright_bessel(a: mp.mpf, b: mp.mpf, x: mp.mpf) -> mp.mpf:
    """`specs.r_log_wright_bessel`, with ln Phi = x - ln Gamma(b) for a = 0."""
    if a == 0 and x > 0:
        return x + mp.log(mp.rgamma(b))
    return specs.r_log_wright_bessel(a, b, x)


def _wright_bessel(a: mp.mpf, b: mp.mpf, x: mp.mpf) -> mp.mpf:
    """`specs.r_wright_bessel`, with exp(x) / Gamma(b) for a = 0."""
    if a == 0 and x > 0:
        return mp.exp(x) * mp.rgamma(b)
    return specs.r_wright_bessel(a, b, x)


WRIGHT_DOMAIN: Final = Rule(
    lambda a: min(a) < 0, math.nan, "only a, b, x >= 0 are supported, like in xsf"
)


def _voigt_arb(x: flint.arb, sigma: flint.arb, gamma: flint.arb) -> flint.arb:
    """Re w(z) / (sigma sqrt(2 pi)), with z = (x + i gamma) / (sigma sqrt(2)); the limits for
    sigma = 0 and gamma = 0 like `_voigt`."""
    if sigma == 0:
        return gamma / (flint.arb.pi() * (x * x + gamma * gamma))
    if gamma == 0:
        return (-((x / sigma) ** 2) / 2).exp() / (
            abs(sigma) * (2 * flint.arb.pi()).sqrt()
        )
    w = _wofz_arb(flint.acb(x, gamma) / (sigma * flint.arb(2).sqrt()))
    return w.real / (sigma * (2 * flint.arb.pi()).sqrt())


def _iv_ratio_rules(complement: bool) -> tuple[Rule, ...]:
    """The conventions of iv_ratio(v, x) = I_v(x) / I_{v-1}(x), or of its complement."""
    small, large = (
        (1.0, 0.0) if complement else (0.0, 1.0)
    )  # the ratio's limits 0 and 1
    return (
        NAN,
        Rule(
            lambda a: a[0] < 0.5 or a[1] < 0,
            math.nan,
            "NaN outside of v >= 0.5, x >= 0",
        ),
        Rule(
            lambda a: math.isinf(a[0]) and math.isinf(a[1]),
            math.nan,
            "no limit for v, x -> +inf",
        ),
        Rule(lambda a: a[1] == 0, small, "the limit for x -> 0+, with v >= 0.5"),
        Rule(
            lambda a: math.isinf(a[0]),
            small,
            "the limit for v -> +inf, with x finite",
        ),
        Rule(
            lambda a: math.isinf(a[1]),
            large,
            "the limit for x -> +inf, with v finite",
        ),
    )


def _iv_ratio_ball(v: flint.arb, x: flint.arb) -> flint.arb:
    """x / (2 v) 0F1(; v + 1; x^2 / 4) / 0F1(; v; x^2 / 4), which also works for a huge v."""
    z = x * x / 4
    return x / (2 * v) * z.hypgeom([], [v + 1]) / z.hypgeom([], [v])


def _iv_ratio_c_ball(v: flint.arb, x: flint.arb) -> flint.arb:
    return 1 - _iv_ratio_ball(v, x)


def _wright_asymptotic(a: float, b: float, x: float) -> bool:
    """Whether xsf evaluates the asymptotic expansion, which needs a < 5 here."""
    return (a * x) ** (1 / (1 + a)) >= 14 + b * b / (2 * (1 + a))


WRIGHT_NAN: Final = Override(
    lambda a: (
        0.5 <= a[0] <= 1.8
        and 100 <= a[1] < math.inf
        and 1e5 <= a[2] < math.inf
        and not _wright_asymptotic(*a)
    ),
    "xsf returns NaN for 1/2 <= a <= 9/5, b >= 100, x >= 10^5, where the asymptotic expansion "
    "does not apply",
)


def _wright_series(a: float, b: float, x: float) -> bool:
    """Whether `wright_bessel` in xsf-rust gets past its special cases and its series for a
    small a or x <= 2 (the `rgamma_zero` and `EXP_INF` constants of xsf)."""
    return (
        0 < a < 178.47241115886637
        and 0 <= b < 178.47241115886637
        and 2 < x < math.inf
        and not (a <= 1e-3 and b <= 50 and x <= 9)
        and not (a <= 1e-4 and b <= 70 and x <= 100)
        and not (a <= 1e-5 and b <= 170 and x < 709.7827128933841)
    )


def _wright_peak(a: float, b: float, x: float) -> bool:
    """Whether `wright_bessel` in xsf-rust sums the terms around the maximum term for a < 5."""
    if not (_wright_series(a, b, x) and a < 5):
        return False
    z = (a * x) ** (1 / (1 + a))
    return (
        z / a <= 1e8
        and z < 14 * (1 + a) + 4 * b * b / (1 + a)
        and (b >= 80 or (z >= 10 and x >= 200))
    )


def _wright_integral(a: float, b: float, x: float) -> bool:
    """Whether `wright_bessel` in xsf-rust evaluates the integral representation."""
    if not (_wright_series(a, b, x) and a < 5) or _wright_peak(a, b, x):
        return False
    z = (a * x) ** (1 / (1 + a))
    beyond = z / a > 1e8 or z >= 14 * (1 + a) + 4 * b * b / (1 + a)
    return not (_wright_asymptotic(a, b, x) and beyond)


WRIGHT_FLOORS: Final = (
    Floor(
        lambda a: _wright_series(*a) and a[0] >= 5,
        8e-16,
        "for a >= 5, xsf-rust centers the series on the maximum term including b, and sizes it "
        "by the spread of the terms",
    ),
    Floor(
        lambda a: _wright_peak(*a),
        5e-14,
        "for large b, or large Z and x, or Z near xsf's threshold of the asymptotic expansion, "
        "xsf-rust sums the terms around the maximum term instead",
    ),
    Floor(
        lambda a: _wright_integral(*a),
        7e-16,
        "the error of the integral representation depends on the platform's libm (2.66e-15 at "
        "(4.999, 1e-10, 9.999999999999778) on macOS)",
    ),
)


def _hypu_poly(m: int, b: Fraction, x: Fraction) -> Fraction:
    """U(-m, b, x), a polynomial (DLMF 13.2.7), exactly."""
    terms = (
        math.comb(m, s) * math.prod(b + s + j for j in range(m - s)) * (-x) ** s
        for s in range(m + 1)
    )
    return (-1) ** m * sum(terms, Fraction(0))


def _mpf(q: Fraction) -> mp.mpf:
    return mp.mpf(q.numerator) / q.denominator


def _hypu_terms(a: mp.mpf, x: mp.mpf) -> mp.mpf:
    """The size of the terms of U(a, 1, x) = -1 / Gamma(a) sum_k (a)_k x^k / k!^2 (ln x +
    psi(a + k) - 2 psi(1 + k)) (DLMF 13.2.9), with the three parts apart, since ln x and
    psi(a) + 2 gamma cancel near a zero of U."""
    t, psi_a, psi_1, ln_x = mp.rgamma(a), mp.digamma(a), -mp.euler, abs(mp.log(x))
    s, k = mp.mpf(0), 0
    while True:
        term = abs(t) * (ln_x + abs(psi_a) + 2 * abs(psi_1))
        s += term
        if k > abs(a) and term <= mp.eps * s:
            return s
        t *= (a + k) * x / (k + 1) ** 2
        psi_a += 1 / (a + k)
        psi_1 += mp.mpf(1) / (k + 1)
        k += 1


def _hypu(a: mp.mpf, b: mp.mpf, x: mp.mpf) -> Value:
    """`specs.r_hypu`, but exactly where a or a - b + 1 is a non-positive integer (DLMF 13.2.7
    and 13.2.40, where mpmath fails at the zeros), and U(a, b, 0) for b < 1 (DLMF 13.2.20-22).
    For b = 1 and 0 < x < 1, the error scale is the size of the terms that cancel: those of
    DLMF 13.2.9, and of DLMF 13.3.7 where xsf's `hypu` uses it."""
    fa, fb, fx = (Fraction(float(v)) for v in (a, b, x))
    if _nonpos_int(float(a)):
        return _mpf(_hypu_poly(-int(fa), fb, fx))
    if (aa := fa - fb + 1) <= 0 and aa.denominator == 1:
        return x ** (1 - b) * _mpf(_hypu_poly(-int(aa), 2 - fb, fx))
    if x == 0:
        return mp.gamma(1 - b) * mp.rgamma(_mpf(aa))
    u = specs.r_hypu(a, b, x)
    if b != 1 or x >= 1:
        return u
    scale = _hypu_terms(a, x)
    if -0.25 < a < 0.3:
        # U(a, 1, x) = (x + 2a + 1) U(a + 1, 1, x) - (a + 1)^2 U(a + 2, 1, x)
        rec = abs(x + 2 * a + 1) * abs(specs.r_hypu(a + 1, b, x))
        scale = max(scale, rec + (a + 1) ** 2 * abs(specs.r_hypu(a + 2, b, x)))
    return u, scale


def _hypu_arb(a: flint.arb, b: flint.arb, x: flint.arb) -> flint.arb:
    return x.hypgeom_u(a, b)


def _chgul_digits(a: float, b: float, x: float) -> int:
    """The digits that xsf-rust's chgul estimates for its DLMF 13.7.3 sum (a and a - b + 1 not
    0, -1, ...), which `chgu` returns from 12 digits, before it tries chgubi."""
    hu = r = rmax = 1.0
    r0 = ra = 0.0
    for k in range(1, 26):
        r = -r * (a + k - 1) * (a - b + k) / (k * x)
        ra = abs(r)
        if k > 5 and ra >= r0 or ra < 1e-15:
            break
        r0 = ra
        hu += r
        rmax = max(rmax, ra)
    return int(-math.log10(max(ra, rmax * 1e-15) / min(abs(hu), 1.0)))


def _hypu_x0_sign(a: tuple[float, ...]) -> int:
    """The sign of U(a, b, x) ~ 1 / Gamma(a) for x -> 0+ with b >= 1, or 0 where it is finite
    (a = 0, -1, ...) or x != 0."""
    if a[2] != 0 or a[1] < 1 or _nonpos_int(a[0]):
        return 0
    return 1 if a[0] > 0 else (-1) ** (math.floor(-a[0]) + 1)


HYPU_RULES: Final = (
    NAN,
    Rule(lambda a: a[2] < 0, math.nan, "NaN for x < 0, outside of the domain of xsf"),
    *(
        Rule(
            lambda a, s=s: _hypu_x0_sign(a) == s,
            math.copysign(math.inf, s),
            "the limit for x -> 0+ with b >= 1 and a != 0, -1, ... (DLMF 13.2.16-19)",
        )
        for s in (1, -1)
    ),
)

HYPU_OVERRIDES: Final = (
    Override(
        lambda a: (
            _hypu_x0_sign(a) == -1 or a[2] == 0 and a[1] > 1 and _nonpos_int(a[0])
        ),
        "xsf returns +inf at x = 0 with b >= 1, also where U is finite or -inf",
    ),
    Override(
        lambda a: (
            a[2] == 0
            and a[1] < 1
            and _nonpos_int(1 - a[1] + a[0])
            and (Fraction(a[0]) - Fraction(a[1])).denominator != 1
        ),
        "xsf returns 0 at x = 0 for b < 1 where its 1 - b + a rounds onto a pole of Gamma",
    ),
    Override(
        lambda a: a in {(1.0, -10.0, 10.0), (1.0, -1.0, 10.0), (1.0, 1.0, 10.0)},
        "xsf's chgubi loses digits for U(1, b, 10) with b = -10, -1, 1 (scipy/xsf#380)",
    ),
)


def _wright_wrong(why: str, *points: tuple[float, float, float]) -> Override:
    """At these inputs (found from the tables), xsf is wrong."""
    return Override(lambda a: a in points, why)


WRIGHT_OVERRIDES: Final = (
    WRIGHT_NAN,
    _wright_wrong(
        "xsf's quadrature of the integral representation has a relative error of 5e-8",
        (0.5, 100.0, 1000.0),
    ),
    _wright_wrong(
        "xsf's asymptotic expansion has a relative error of 8e-8 (its source says up to 1e-7)",
        (0.5, 10.0, 709.7827128933841),
    ),
)

LOG_WRIGHT_OVERRIDES: Final = (
    WRIGHT_NAN,
    _wright_wrong(
        "xsf returns NaN, where its quadrature of the integral representation fails",
        (0.001, 100.0, 1e5),
    ),
    _wright_wrong(
        "xsf's quadrature of the integral representation has a relative error of 4e-4",
        (0.001, 50.0, 1e5),
    ),
    _wright_wrong(
        "xsf takes the ln of a series that rounds to about 1 + 1e-10, without ln_1p",
        (10.0, 1.0, 1e-3),
        (10.0, 2.0, 1e-3),
    ),
    _wright_wrong(
        "xsf's Taylor series around the maximum term is off by up to 1e-4 for a = 5, b = 100, "
        "x >= 1e12",
        (5.0, 100.0, 1e12),
        (5.0, 100.0, 1e20),
    ),
)


def _pbdv(v: mp.mpf, x: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """D_v(x) and D_v'(x) = x / 2 D_v(x) - D_{v+1}(x) (DLMF 12.8.3) from the power series of
    U(-v - 1/2, x) at 0 (DLMF 12.4), since `mp.pcfd` does not converge for some integers v."""

    def u(a: mp.mpf) -> mp.mpf:
        c = mp.sqrt(mp.pi) * mp.mpf(2) ** (-a / 2)
        u0 = c / mp.mpf(2) ** (mp.mpf(1) / 4) * mp.rgamma(a / 2 + mp.mpf(3) / 4)
        u1 = -c * mp.mpf(2) ** (mp.mpf(1) / 4) * mp.rgamma(a / 2 + mp.mpf(1) / 4)
        s = x * x / 2
        even = mp.hyp1f1(a / 2 + mp.mpf(1) / 4, mp.mpf(1) / 2, s, zeroprec=mp.mp.prec)
        odd = mp.hyp1f1(a / 2 + mp.mpf(3) / 4, mp.mpf(3) / 2, s, zeroprec=mp.mp.prec)
        return mp.exp(-x * x / 4) * (u0 * even + u1 * x * odd)

    with specs.wdps(30 + int(x * x / 2)):
        a = -v - mp.mpf(1) / 2
        d = u(a)
        return d, x / 2 * d - u(a - 1)


def _pcf_u(a: flint.acb, z: flint.acb) -> flint.acb:
    """U(a, z) from its power series at 0 (DLMF 12.4)."""
    quarter = flint.acb(1) / 4
    u0 = flint.acb.pi().sqrt() * (-(a / 2 + quarter) * flint.acb(2).log()).exp()
    u0 *= (a / 2 + 3 * quarter).rgamma()
    u1 = -flint.acb.pi().sqrt() * (-(a / 2 - quarter) * flint.acb(2).log()).exp()
    u1 *= (a / 2 + quarter).rgamma()
    s = z * z / 2
    even = s.hypgeom_1f1(a / 2 + quarter, flint.acb(1) / 2)
    odd = s.hypgeom_1f1(a / 2 + 3 * quarter, flint.acb(3) / 2)
    return (-z * z / 4).exp() * (u0 * even + u1 * z * odd)


def _pcf_du(a: flint.acb, z: flint.acb) -> flint.acb:
    """U'(a, z) = z / 2 U(a, z) - U(a - 1, z) (DLMF 12.8.3)."""
    return z / 2 * _pcf_u(a, z) - _pcf_u(a - 1, z)


def _pbdv_arb(i: int) -> Callable[[flint.arb, flint.arb], flint.arb]:
    """D_v(x) = U(-v - 1/2, x) (i = 0) or its derivative (i = 1)."""

    def arb(v: flint.arb, x: flint.arb) -> flint.arb:
        a, z = flint.acb(-v - flint.arb(1) / 2), flint.acb(x)
        return (_pcf_du if i else _pcf_u)(a, z).real

    return arb


def _pair_ref(
    ref: Callable[[mp.mpf, mp.mpf], tuple[mp.mpf, mp.mpf]], i: int
) -> Callable[[mp.mpf, mp.mpf], mp.mpf]:
    """Output `i` of a reference with two outputs."""
    return lambda a, b: ref(a, b)[i]


def _pbwa_arb(i: int) -> Callable[[flint.arb, flint.arb], flint.arb]:
    """W(a, x) (i = 0) or its derivative (i = 1) from U(i a, x e^(-i pi / 4)) (DLMF 12.14)."""

    def arb(a: flint.arb, x: flint.arb) -> flint.arb:
        ia = flint.acb(0, a)
        rho = flint.arb.pi() / 8 + (flint.acb(1) / 2 + ia).lgamma().imag / 2
        k = (1 + (2 * flint.arb.pi() * a).exp()).sqrt() - (flint.arb.pi() * a).exp()
        phase = flint.acb(-flint.arb(1) / 4).exp_pi_i()
        z = x * phase
        u = (_pcf_du(ia, z) * phase) if i else _pcf_u(ia, z)
        e = flint.acb(0, rho).exp()
        return (k / 2).sqrt() * (flint.arb.pi() * a / 4).exp() * 2 * (e * u).real

    return arb


def _pbdv_recurrence(a: tuple[float, ...]) -> bool:
    """Whether xsf-rust computes D_v(x) by recurrence from an order <= -20, like its `pbdv`."""
    v, x = a
    lo, hi = (1.0, 36.5) if v < 0 else (4.5, 53.0)
    return abs(v) < 20 and (lo < x < 12 or hi <= x < 60)


PBDV_FLOORS: Final = (
    Floor(
        _pbdv_recurrence,
        2e-15,
        "xsf-rust recurs up from an order <= -20, where xsf sums `dvsa` or `dvla`",
    ),
    Floor(
        lambda a: -10 <= a[1] < -5.8 and -30 <= a[0] < -1e-15,
        2e-15,
        "xsf-rust sums `dvsa` (and D_{v-1}) instead of `dvla`",
    ),
)


def _pmv(m: mp.mpf, v: mp.mpf, x: mp.mpf) -> mp.mpf:
    """The Ferrers function P_v^m(x) for an integer m (DLMF 14.3.1), as the product of
    ((1 + x) / (1 - x))^(m / 2) and the regularized 2F1(v + 1, -v; 1 - m; (1 - x) / 2), which
    for m > 0 is (v + 1)_m (-v)_m / m! z^m 2F1(v + m + 1, m - v; m + 1; z); for an integer v
    and x < 0 with P_n^m(-x) = (-1)^(n + m) P_n^m(x) for its degree n (DLMF 14.9.10)."""
    with specs.wdps(specs.xd(v) + 20):
        if x < 0 and v == int(v):
            n = int(v) if v >= 0 else int(-v) - 1
            return (-1) ** (n + int(m)) * _pmv(m, v, -x)
        z = (1 - x) / 2
        if m >= 0:
            c = mp.rf(v + 1, m) * mp.rf(-v, m) / mp.factorial(m)
            f = mp.hyp2f1(v + m + 1, m - v, m + 1, z, zeroprec=mp.mp.prec)
            return (1 - x * x) ** (m / 2) / mp.mpf(2) ** m * c * f
        f = mp.hyp2f1(v + 1, -v, 1 - m, z, zeroprec=mp.mp.prec)
        return ((1 - x) / (1 + x)) ** (-m / 2) * f / mp.factorial(-m)


def _pmv_scaled(m: mp.mpf, v: mp.mpf, x: mp.mpf) -> Value:
    """P_v^m(x), with the scale M = hypot(P, P_theta / (n + 1/2)) for x = cos(theta) near the
    zero P_n^m(0) = 0 of an integer degree n with n + m odd and positive (DLMF 14.5.1), where
    the sum of xsf cancels down to an error of about M. M is the local amplitude of its
    oscillation in theta, as for large n, P_n^m(cos theta) is J_m((n + 1/2) theta) times a
    slowly varying factor (DLMF 14.15.11), whose zeros are about pi / (n + 1/2) apart: "near" is
    within that of theta = pi / 2. P_theta = -sin(theta) P' with DLMF 14.10.4."""
    with specs.wdps(specs.xd(v) + 20):
        p = _pmv(m, v, x)
        n = _pmv_degree(float(v))
        zero = n.is_integer() and n + m > 0 and (n + m) % 2 == 1
        if not zero or abs(x) >= 1 or abs(mp.asin(x)) >= mp.pi / (n + mp.mpf(1) / 2):
            return p
        dp = -((m - v - 1) * _pmv(m, v + 1, x) + (v + 1) * x * p) / mp.sqrt(1 - x * x)
        return p, mp.hypot(p, dp / (n + mp.mpf(1) / 2))


def _pmv_arb(m: flint.arb, v: flint.arb, x: flint.arb) -> flint.arb | None:
    """The same from the regularized 2F1 of Arb, or None at x = +-1, where its factor is 0 or
    infinite."""
    if abs(x) == 1:
        return None
    f = ((1 - x) / 2).hypgeom_2f1(v + 1, -v, 1 - m, regularized=True)
    return ((1 + x) / (1 - x)) ** (m / 2) * f


def _pmv_degree(v: float) -> float:
    """The degree n of P_v = P_(-v-1) (DLMF 14.9.5) that xsf uses, v for v >= 0."""
    return -v - 1 if v < 0 else v


PMV_RULES: Final = (
    NAN,
    Rule(
        lambda a: math.isfinite(a[0]) and not a[0].is_integer(),
        math.nan,
        "xsf takes the order m as a double, and it is NaN for a non-integer m",
    ),
    Rule(
        lambda a: abs(a[2]) > 1 and a[0] != 0,
        math.nan,
        "xsf has NaN for |x| > 1 and m != 0, where the Ferrers function has its sqrt(1 - x^2)",
    ),
    Rule(
        lambda a: (
            a[0] < 0
            and _pmv_degree(a[1]).is_integer()
            and _pmv_degree(a[1]) + a[0] + 1 <= 0
        ),
        math.nan,
        "xsf has no formula for a negative order -k and an integer degree n < k, where "
        "DLMF 14.9.3 does not help",
    ),
)


PBWA_NAN: Final = Rule(
    lambda a: abs(a[0]) > 5 or abs(a[1]) > 5,
    math.nan,
    "xsf and xsf-rust only compute W(a, x) for |a| <= 5 and |x| <= 5",
)


def _hyp2f1_degree(a: float, b: float, c: float) -> int | None:
    """The degree of 2F1 if a or b is a non-positive integer (the one closer to 0 ends it),
    unless c is a pole before that."""
    ends = [-int(v) for v in (a, b) if _nonpos_int(v)]
    m = min(ends, default=None)
    return None if m is None or _nonpos_int(c) and -m < c else m


def _hyp2f1_poly(args: tuple[float, ...], m: int) -> Fraction:
    """The polynomial 2F1 of degree m, also for c = -m - l with l >= 0 (DLMF 15.2.5), exactly."""
    a, b, c, x = map(Fraction, args)
    s = t = Fraction(1)
    for k in range(m):
        t *= (a + k) * (b + k) / ((c + k) * (k + 1)) * x
        s += t
    return s


def _hyp2f1(a: mp.mpf, b: mp.mpf, c: mp.mpf, x: mp.mpf) -> mp.mpf:
    """`specs.r_hyp2f1`, but a polynomial exactly, and the real part for x > 1, where the rules
    leave only real values (whose imaginary part is a rounding error)."""
    args = tuple(map(float, (a, b, c, x)))
    if (m := _hyp2f1_degree(*args[:3])) is not None:
        f = _hyp2f1_poly(args, m)
        return mp.mpf(f.numerator) / f.denominator
    return mp.re(specs.r_hyp2f1(a, b, c, x))


def _hyp2f1_arb(a: flint.arb, b: flint.arb, c: flint.arb, x: flint.arb) -> flint.arb:
    """The real part of 2F1 of Arb, or NaN, which Arb does not confirm, if it is not real.
    But a polynomial exactly, where Arb can be NaN (for c = a or c = b) or too wide (for large
    degrees) and has no exact zeros."""
    args = tuple(map(float, (a, b, c, x)))
    if (m := _hyp2f1_degree(*args[:3])) is None:
        f = flint.acb(x).hypgeom_2f1(a, b, c)
        return f.real if 0 in f.imag else flint.arb.nan()
    f = _hyp2f1_poly(args, m)
    return flint.arb(flint.fmpq(f.numerator, f.denominator))


def _hyp2f1_pole(a: tuple[float, ...]) -> bool:
    """Whether c is a non-positive integer before which no a or b ends the series, at x != 0."""
    return _nonpos_int(a[2]) and _hyp2f1_degree(*a[:3]) is None and a[3] != 0


def _hyp2f1_complex(a: tuple[float, ...]) -> bool:
    """Whether xsf-rust has NaN for x > 1, where 2F1 is generally complex (its fallback): unless
    a, b, c - a or c - b (rounded) is a non-positive integer up to `HYP2F1_EPS`, where it is a
    polynomial, or one times a power of 1 - x (DLMF 15.8.1), and else +inf, or c is near a pole
    (+inf, relative to it as in xsf-rust, so c = -2^-55 isn't), or c - a or c - b rounds to
    c > 0 (+inf). The table has none of these +inf rows, nor rows where a parameter is only near
    an integer, as Arb confirms neither for the inputs as given. Some others are real too, such
    as 2F1(2, 2; 4; 2) = -3, whose logarithmic term vanishes."""
    p, q, c, x = a
    poly = any(_near_nonpos_int(t) for t in (p, q, c - p, c - q))
    pole = c <= 0 and abs(c - round(c)) <= HYP2F1_EPS * min(abs(round(c)), 1)
    lost = c > 0 and c in (c - p, c - q)
    return x > 1 and not (poly or pole or lost)


def _hyp2f1_diverges(a: tuple[float, ...]) -> bool:
    """Whether 2F1 diverges to +inf for x -> 1-, for a, b, c > 0 and c - a - b <= 0."""
    p, q, c = map(Fraction, a[:3])
    return a[3] == 1 and min(a[:3]) > 0 and c - p - q <= 0


HYP2F1_EPS: Final = 1e-13  # as in xsf: up to this, a parameter counts as an integer


def _near_nonpos_int(t: float) -> bool:
    """Whether t is a non-positive integer up to `HYP2F1_EPS`, as in xsf."""
    return t <= 0 and abs(t - round(t)) < HYP2F1_EPS


def _xsf_c_eq_wrong(a: tuple[float, ...]) -> bool:
    """Whether c = a = -m or c = b = -m with no a or b ending the series earlier, at x != 0,
    except where xsf has a polynomial: for a = -1, and by `hyp2f1_neg_c_equal_bc` for b = c and
    |x| < 1 or x = -1."""
    p, q, c, x = a
    own = p == -1 or q == c and (abs(x) < 1 or x == -1)
    return (
        _nonpos_int(c)
        and c in (p, q)
        and _hyp2f1_degree(p, q, c) == -c
        and x != 0
        and not own
    )


def _xsf_degree_1_zero(a: tuple[float, ...]) -> bool:
    """Whether the polynomial 1 - a x / c for b = -1 (but a != -1) and an integer c < 0 is 0."""
    p, q, c, x = a
    neg_int_c = c < 0 and c.is_integer()
    return q == -1 and p != -1 and neg_int_c and Fraction(p) * Fraction(x) == c


def _xsf_euler_rounds(a: tuple[float, ...]) -> bool:
    """Whether xsf takes the Euler transformation (for x != 0, c - a - b <= -1, x <= 1 or an
    integer c - a - b, and a, b not polynomials) with c - a or c - b near a negative integer,
    where only its rounded value is one."""
    p, q, c = map(Fraction, a[:3])
    d = a[2] - a[0] - a[1]  # rounded, as in xsf

    def rounds(t: Fraction) -> bool:
        return round(t) < 0 and _near_nonpos_int(float(t)) and t.denominator != 1

    euler = a[3] != 0 and d <= -1 and (a[3] <= 1 or abs(d - round(d)) <= HYP2F1_EPS)
    poly = _near_nonpos_int(a[0]) or _near_nonpos_int(a[1])
    return euler and not poly and (rounds(c - p) or rounds(c - q))


def _hyp2f1_returns_early(a: tuple[float, ...]) -> bool:
    """Whether xsf-rust returns before it tries `hyp2f1_recur` (in `hyp2f1_impl`): for x = 0,
    a = 0 or b = 0, c = inf and x <= 1, a = -1 or b = -1 and c != 0, or an integer c <= 0 with
    c = a or c = b, or at x = 1 with an integer a or b < 0 above c."""
    p, q, c, x = a
    ends = any(_nonpos_int(t) and t > c for t in (p, q))
    int_c = _nonpos_int(c) and (c in (p, q) or x == 1 and ends)
    degree_1 = -1 in (p, q) and c != 0
    return x == 0 or 0 in (p, q) or (c == math.inf and x <= 1) or degree_1 or int_c


def _hyp2f1_recur_first(a: tuple[float, ...]) -> bool:
    """Whether xsf-rust first tries `hyp2f1_recur` (`in_recur_regime`): for max(|a|, |b|, |c|)
    >= 10, c - a - b (x >= 0.9) or b - a (x <= -9) within 1e-3 of an integer, and a, b not near
    a non-positive integer, a, b, c - a or c - b near (but not at) a negative integer, or a tiny
    negative a or b. Unless it returns before (`_hyp2f1_returns_early`), x >= 1 for a c that
    isn't a non-positive integer (where `hyp2f1_recur` has no result), c is near (but not at) a
    non-positive integer and not tiny, or a, b, c - a or c - b (rounded) is 0 or a negative
    subnormal."""
    p, q, c, x = a
    fp, fq, fc = map(Fraction, (p, q, c))
    exact = (fp, fq, fc - fp, fc - fq)
    params = tuple(map(float, exact))

    def tiny(t: float) -> bool:
        return sys.float_info.min <= -t < HYP2F1_EPS

    def near_poly(t: Fraction) -> bool:
        return float(t) < -0.5 and _near_nonpos_int(float(t)) and t.denominator != 1

    diff, w = (q - p, x / (x - 1)) if x < 0 else (c - p - q, x)
    poly = _near_nonpos_int(p) or _near_nonpos_int(q)
    near_log = w >= 0.9 and abs(diff - round(diff)) <= 1e-3 and not poly
    large = max(map(abs, a[:3])) >= 10
    regime = large or near_log or any(map(near_poly, exact)) or tiny(p) or tiny(q)
    at_pole = _near_nonpos_int(c) and not (tiny(c) or c.is_integer())
    degenerate = any(-sys.float_info.min < t <= 0 for t in params)
    no_result = x >= 1 and not _nonpos_int(c)
    skipped = at_pole or degenerate or no_result or _hyp2f1_returns_early(a)
    return regime and not skipped


def _hyp2f1_recur_large(a: tuple[float, ...]) -> bool:
    """Whether xsf-rust first tries `hyp2f1_recur` for max(|a|, |b|, |c|) >= 10."""
    return max(map(abs, a[:3])) >= 10 and _hyp2f1_recur_first(a)


def _hyp2f1_neg_c_series(a: tuple[float, ...]) -> bool:
    """Whether xsf-rust first tries the series of the Pfaff transformations for -1 < x < 0, or the
    power series for 0 < x < 1, for c - a - b <= -1 and a c < 0 that isn't an integer, unless
    it first tries `hyp2f1_recur`."""
    p, q, c, x = a
    pole = abs(c - round(c)) <= HYP2F1_EPS
    poly = _near_nonpos_int(p) or _near_nonpos_int(q)
    series = c - p - q <= -1 and c < 0 and not pole and not poly and 0 < abs(x) < 1
    return series and not _hyp2f1_recur_first(a)


HYP2F1_RULES: Final = (
    Rule(
        _hyp2f1_pole,
        math.inf,
        "+inf at a pole of c that no a or b ends the series before, for x != 0",
    ),
    Rule(
        _hyp2f1_diverges,
        math.inf,
        "the limit for x -> 1- for a, b, c > 0 and c - a - b <= 0 (DLMF 15.4.21, 15.4.23)",
    ),
    Rule(
        _hyp2f1_complex,
        math.nan,
        "NaN for x > 1 without a real closed form, where 2F1 is generally complex (DLMF 15.2.3)",
    ),
)
HYP2F1_OVERRIDES: Final = (
    Override(_hyp2f1_complex, "xsf has inf for x > 1 where 2F1 is complex (#1531)"),
    Override(
        _xsf_c_eq_wrong,
        "xsf has inf, or (1 - x)^(-b) of DLMF 15.2.6 for c = a and |x| < 1, instead of the "
        "polynomial of DLMF 15.2.5 for c = a = -m or c = b = -m (#1574)",
    ),
    Override(
        _xsf_degree_1_zero,
        "xsf has NaN where the polynomial 1 - a x / c for b = -1 and an integer c < 0 is 0 "
        "(#1573)",
    ),
    Override(
        _xsf_euler_rounds,
        "xsf takes c - a or c - b in its Euler transformation as a negative integer where only "
        "its rounded value is one, and so the polynomial of that degree",
    ),
)
HYP2F1_FLOORS: Final = (
    Floor(
        _hyp2f1_recur_large,
        3e-15,
        "xsf-rust first tries `hyp2f1_recur` for large parameters, and so rounds differently "
        "(1.0e-14 at (24.17, 2.23, 1.49, 0.822), where xsf has 7.0e-16)",
    ),
    Floor(
        _hyp2f1_neg_c_series,
        2e-15,
        "xsf-rust first tries a series instead of the Euler transformation (scipy/xsf#464), "
        "and so rounds differently (5.8e-15 at (1.8, 1, -1.8, 0.95), where xsf has 7.8e-16)",
    ),
)


SW_MODE: Final = Rule(
    lambda a: not 0 <= a[0] <= a[1], math.nan, "NaN outside of 0 <= m <= n"
)
SW_ANGLE: Final = Rule(
    lambda a: not -1 < a[-1] < 1, math.nan, "NaN outside of -1 < x < 1"
)
SW_PROLATE_X: Final = Rule(lambda a: not a[-1] > 1, math.nan, "NaN outside of x > 1")


def _spheroidal(
    name: str, rules: tuple[Rule, ...], cvs: tuple[bool, ...] = (False, True)
) -> dict[tuple[str, str], Entry]:
    """The tables of a spheroidal wave function (`prolate_aswfa`, ..., `oblate_radial2`) and
    its derivative: of (m, n, c, x) (`_nocv`), and of (m, n, c, cv, x) for `True` in `cvs`."""
    shape, kind = name.split("_")
    kd = 1 if shape == "prolate" else -1
    angular = kind == "aswfa"
    if angular:
        refs, residual = spheroidal.angular, spheroidal.angular_residual(kd)
    else:
        refs = functools.partial(spheroidal.radial, int(kind[-1]))
        residual = spheroidal.wronskian_residual(kd)

    def entry(cv: bool) -> Entry:
        outs = (
            Out(
                f"{name}{'' if cv else '_nocv'}.{i}",
                (spheroidal.with_cv(refs(kd, i), kd, angular) if cv else refs(kd, i),),
                rules=rules,
            )
            for i in range(2)
        )
        return Entry(tuple(outs), Identity(lambda a: True, residual))

    return {
        (f"{name}{'' if cv else '_nocv'}", f"d_d_d_d{'_d' * cv}-d_d"): entry(cv)
        for cv in cvs
    }


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
    ("beta", "d_d-d"): Entry(
        (
            Out(
                "beta",
                (specs.r_beta,),
                (_beta_arb,),
                (NAN, *BETA_POLES),
                BETA_OVERRIDES,
                BETA_FLOORS,
            ),
        )
    ),
    ("binom", "d_d-d"): Entry(
        (
            Out(
                "binom",
                (specs.r_binom,),
                (_binom_arb,),
                (
                    NAN,
                    Rule(
                        lambda a: a[0] < 0 and a[0].is_integer(),
                        math.nan,
                        "xsf leaves binom(n, k) undefined for a negative integer n",
                    ),
                ),
                BINOM_OVERRIDES,
                BINOM_FLOORS,
            ),
        )
    ),
    ("poch", "d_d-d"): _real(
        "pow_rising",
        specs.r_poch,
        _sums_exact(_gamma_ratio),
        (
            NAN,
            Rule(
                lambda a: _nonpos_int(a[0]) and not _sum_nonpos_int(a),
                0.0,
                "1 / Gamma(a) = 0 at a pole of Gamma(a) that Gamma(a + m) does not share",
            ),
            Rule(
                lambda a: _sum_nonpos_int(a) and not _nonpos_int(a[0]),
                math.inf,
                "a pole whose one-sided limits have opposite signs",
            ),
        ),
        POCH_OVERRIDES,
    ),
    ("betaln", "d_d-d"): Entry(
        (
            Out(
                "betaln",
                (specs.r_betaln,),
                (_betaln_arb,),
                (
                    NAN,
                    Rule(_gamma_pole, math.inf, "ln |B(a, b)| -> +inf at a pole"),
                    Rule(
                        lambda a: _sum_nonpos_int(a) and not any(map(_nonpos_int, a)),
                        -math.inf,
                        "B(a, b) = 0 if a + b is a pole of Gamma but a and b are not",
                    ),
                ),
                floors=(Floor(_small_b, 8e-16, SMALL_B),),
            ),
        )
    ),
    ("lanczos_sum_expg_scaled", "d-d"): _real(
        "lanczos_sum_expg_scaled", _lanczos, _lanczos, (NAN,)
    ),
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
    ("digamma", "cd-cd"): _complex(
        "digamma", specs._c(mp.digamma), flint.acb.digamma, (COMPLEX_POLE,)
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
    ("erfcx", "cd-cd"): _complex("erfcx", specs._c(specs.c_erfcx), _erfcx_arb),
    ("dawsn", "d-d"): Entry(
        (Out("dawsn", (specs.r_dawsn,), (_dawsn_arb,), _limits(-0.0, 0.0)),)
    ),
    ("dawsn", "cd-cd"): _complex("dawsn", specs._c(specs.c_dawsn), _dawsn_arb),
    ("wofz", "cd-cd"): _complex(
        "wofz", specs._c(specs.c_wofz), _wofz_arb, (*WOFZ_AT_INF, COMPLEX_NAN)
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
    **{
        (f"modified_fresnel_{name}", "d-cd_cd"): Entry(
            (_modified_fresnel_out(sign, 0), _modified_fresnel_out(sign, 1))
        )
        for sign, name in enumerate(("plus", "minus"))
    },
    ("airy", "d-d_d_d_d"): Entry(tuple(_airy_out(i) for i in range(4))),
    ("airye", "d-d_d_d_d"): Entry(tuple(_airy_out(i, scaled=True) for i in range(4))),
    ("fresnel", "d-d_d"): Entry(
        (
            Out("fresnel.0", (specs.r_fresnel_s,), (_fresnel_arb(0),), FRESNEL_RULES),
            Out("fresnel.1", (specs.r_fresnel_c,), (_fresnel_arb(1),), FRESNEL_RULES),
        )
    ),
    ("owens_t", "d_d-d"): _real(
        "owens_t", _owens_t, _owens_t_arb, (NAN, *OWENS_T_H_INF)
    ),
    ("expm1", "d-d"): _real(
        "expm1", specs.r_expm1, flint.arb.expm1, _limits(-1, math.inf)
    ),
    ("expm1", "cd-cd"): _complex(
        "expm1",
        _expm1,
        flint.acb.expm1,
        EXPM1_AT_INF,
        (
            Override(
                lambda a: LN_MAX < a[0] < math.inf,
                "xsf has inf * 0 = NaN in Im(e^z - 1) where e^x overflows, for a finite"
                f" Re z > {LN_MAX} (#1724)",
            ),
        ),
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
    ("log1p", "cd-cd"): _complex(
        "log1p",
        specs._c(mp.log1p),
        flint.acb.log1p,
        (
            *LOG1P_AT_INF,
            COMPLEX_NAN,
            Rule(
                lambda a: a[0] == -1 and a[1] == 0 and math.copysign(1, a[1]) > 0,
                complex(-math.inf, 0),
                "ln(1 + z) -> -inf + 0i for z -> -1 along Im z = +0, like C's clog",
            ),
        ),
    ),
    ("log1pmx", "d-d"): _real(
        "log1pmx", _log1pmx, lambda x: x.log1p() - x, (LOG1P_DOMAIN,)
    ),
    ("exp2", "d-d"): _real("exp2", specs.r_exp2, lambda x: flint.arb(2) ** x),
    ("exp10", "d-d"): _real("exp10", specs.r_exp10, lambda x: flint.arb(10) ** x),
    **_logistic("d-d", "", LN_MAX),
    **_logistic("f-f", "_f32", LN_MAX_F32),
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
    ("xlogy", "cd_cd-cd"): _complex(
        "xlogy", specs.c_xlogy, lambda x, y: x * y.log(), (COMPLEX_NAN, COMPLEX_X_ZERO)
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
    ("xlog1py", "cd_cd-cd"): _complex(
        "xlog1py", specs.c_xlog1py, lambda x, y: x * y.log1p()
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
    ("sinpi", "cd-cd"): _complex("sinpi", specs._c(mp.sinpi), flint.acb.sin_pi),
    ("cospi", "d-d"): _real(
        "cospi",
        specs.r_cospi,
        flint.arb.cos_pi,
        (Rule(lambda a: a[0] % 1 == 0.5, 0.0, "cos(pi (n + 1/2)) = +0, like in xsf"),),
    ),
    ("cospi", "cd-cd"): _complex("cospi", specs._c(mp.cospi), flint.acb.cos_pi),
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
    **_bessel_v(),
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
    ("expn", "d_d-d"): EXPN,
    ("expn", "p_d-d"): EXPN,
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
    ("expi", "cd-cd"): _complex(
        "expi",
        specs._c(mp.ei),
        flint.acb.ei,
        (COMPLEX_EI_AT_0,),
        (_e1_overflow(lambda a: -a[0]),),
    ),
    ("shichi", "cd-cd_cd"): Entry(
        (
            _complex_out(
                "shichi.0", _off_axis(mp.shi), flint.acb.shi, (), SHICHI_OVERFLOW
            ),
            _complex_out(
                "shichi.1",
                _off_axis(mp.chi),
                flint.acb.chi,
                (COMPLEX_LOG_AT_0,),
                SHICHI_OVERFLOW,
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
    ("ellipkinc", "d_d-d"): _incomplete(
        "ellipkinc",
        mp.ellipf,
        flint.acb.elliptic_f,
        "F",
        (
            Rule(
                lambda a: math.isinf(a[0]) and a[1] == -math.inf,
                math.nan,
                "no limit for phi -> +-inf and m -> -inf",
            ),
            Rule(
                lambda a: a[1] == -math.inf and 0 < a[0] < math.inf,
                0.0,
                "the limit for m -> -inf, with 0 < phi < inf",
            ),
            Rule(
                lambda a: a[1] == 1 and abs(a[0]) > math.pi / 2,
                math.inf,
                "F(phi|1) diverges for |phi| > pi/2, +inf like in xsf and mpmath",
            ),
        ),
        (
            Override(
                lambda a: a[1] == 1 and abs(a[0]) == math.pi / 2,
                "xsf returns +inf for m = 1 and |phi| = fl(pi/2) < pi/2, where F is finite",
            ),
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
    ("ellipj", "d_d-d_d_d_d"): Entry(tuple(_ellipj_out(i) for i in range(4))),
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
    ("struve_h", "d_d-d"): _struve(-1),
    ("struve_l", "d_d-d"): _struve(1),
    ("itstruve0", "d-d"): _real(
        "itstruve0",
        specs.r_itstruve0,
        _integral_arb(lambda t: _struve_arb(0.0, t, -1)),
        (INTEGRAL_AT_0,),
        (_mid_x_override(399),),
    ),
    ("itmodstruve0", "d-d"): _real(
        "itmodstruve0",
        specs.r_itmodstruve0,
        _integral_arb(lambda t: _struve_arb(0.0, t, 1)),
        (INTEGRAL_AT_0,),
    ),
    ("it1j0y0", "d-d_d"): Entry(
        (
            Out(
                "it1j0y0.0",
                (specs.r_itj0,),
                (_integral_arb(lambda t: t.bessel_j(0)),),
                (INTEGRAL_AT_0,),
                (_mid_x_override(446),),
            ),
            Out(
                "it1j0y0.1",
                (specs.r_ity0,),
                (_int_c0_arb(flint.arb.bessel_y, -1),),
                (INTEGRAL_AT_0, _cut("Y_0")),
                (_mid_x_override(446),),
            ),
        )
    ),
    ("it1i0k0", "d-d_d"): Entry(
        (
            Out(
                "it1i0k0.0",
                (specs.r_iti0,),
                (_integral_arb(lambda t: t.bessel_i(0)),),
                (INTEGRAL_AT_0,),
            ),
            Out(
                "it1i0k0.1",
                (specs.r_itk0,),
                (_int_c0_arb(flint.arb.bessel_k, 1),),
                (INTEGRAL_AT_0, _cut("K_0")),
                (_mid_x_override(379),),
            ),
        )
    ),
    ("it2i0k0", "d-d_d"): Entry(
        (
            Out(
                "it2i0k0.0",
                (specs.r_it2i0,),
                (_it2i0_arb,),
            ),
            Out(
                "it2i0k0.1",
                (specs.r_it2k0,),
                rules=(_cut("K_0"),),
                overrides=(
                    Override(
                        lambda a: 1 < a[0] < 40,
                        "xsf loses digits for 1 < x < 40, where xsf-rust differs "
                        "(scipy/xsf#379)",
                    ),
                ),
            ),
        ),
        Identity(lambda a: a[0] > 0, _it2k0_residual),
    ),
    ("itairy", "d-d_d_d_d"): Entry(tuple(_itairy_out(i) for i in range(4))),
    ("gammainc", "d_d-d"): Entry((_gammainc_out(False, 7e-16),)),
    ("gammaincc", "d_d-d"): Entry((_gammainc_out(True, 3e-15),)),
    ("gammaincinv", "d_d-d"): Entry((_gammaincinv_out(False, 2e-15),)),
    ("gammainccinv", "d_d-d"): Entry((_gammaincinv_out(True, 2e-14),)),
    ("chdtr", "d_d-d"): Entry(
        (_distribution("chdtr", False, _half, "(v / 2, x / 2)"),)
    ),
    ("chdtrc", "d_d-d"): Entry(
        (
            _distribution(
                "chdtrc",
                True,
                _half,
                "(v / 2, x / 2)",
                floors=(Floor(lambda a: _igam_fac_rgamma(*a), 3e-16, RGAMMA),),
            ),
        )
    ),
    ("chdtri", "d_d-d"): Entry((_chdtri(),)),
    ("gdtr", "d_d_d-d"): Entry(
        (_distribution("gdtr", False, _gamma_dist, "(b, a x)", (X_NEG,)),)
    ),
    ("gdtrc", "d_d_d-d"): Entry(
        (
            _distribution(
                "gdtrc",
                True,
                _gamma_dist,
                "(b, a x)",
                (X_NEG,),
                (Floor(lambda a: _igam_fac_rgamma(*a), 2e-15, RGAMMA),),
            ),
        )
    ),
    ("pdtr", "d_d-d"): Entry(
        (_distribution("pdtr", True, _poisson, "(floor(k) + 1, m)", (K_M_NEG,)),)
    ),
    ("pdtrc", "d_d-d"): Entry(
        (_distribution("pdtrc", False, _poisson, "(floor(k) + 1, m)", (K_M_NEG,)),)
    ),
    ("pdtri", "p_d-d"): Entry((_pdtri(),)),
    ("bdtr", "d_p_d-d"): Entry((_beta_dist("bdtr", _bdtr_args, BDTR_RULES),)),
    ("bdtrc", "d_p_d-d"): Entry((_beta_dist("bdtrc", _bdtrc_args, BDTRC_RULES),)),
    ("bdtri", "d_p_d-d"): Entry((_bdtri(),)),
    ("fdtr", "d_d_d-d"): Entry((_beta_dist("fdtr", _fdtr_args, (FDTR_DOMAIN,)),)),
    ("fdtrc", "d_d_d-d"): Entry((_beta_dist("fdtrc", _fdtrc_args, (FDTR_DOMAIN,)),)),
    ("fdtri", "d_d_d-d"): Entry((_fdtri(),)),
    ("nbdtr", "p_p_d-d"): Entry((_beta_dist("nbdtr", _nbdtr_args, NBDTR_RULES),)),
    ("nbdtrc", "p_p_d-d"): Entry((_beta_dist("nbdtrc", _nbdtrc_args, NBDTR_RULES),)),
    ("besselpoly", "d_d_d-d"): Entry(
        (
            Out(
                "besselpoly",
                (_besselpoly,),
                (_besselpoly_arb,),
                BESSELPOLY_RULES,
                BESSELPOLY_OVERRIDES,
            ),
        )
    ),
    ("hyperu", "d_d_d-d"): Entry(
        (
            Out(
                "hypu",
                (_hypu,),
                (_precise(_hypu_arb),),
                HYPU_RULES,
                HYPU_OVERRIDES,
                (
                    Floor(
                        lambda a: (
                            a[1] == 1
                            and 0 < a[0] < 1
                            and a[2] >= 1
                            and _chgul_digits(*a) < 12
                        ),
                        2e-15,
                        "the series of chgubi cancels, with digamma and gamma for a < 1 "
                        "(unlike in xsf, scipy/xsf#380), for b = 1, 0 < a < 1 and x >= 1 "
                        "where chgul has fewer than 12 digits",
                    ),
                ),
            ),
        )
    ),
    ("voigt_profile", "d_d_d-d"): Entry(
        (Out("voigt_profile", (_voigt,), (_voigt_arb,), VOIGT_RULES),)
    ),
    ("iv_ratio", "d_d-d"): Entry(
        (
            Out(
                "iv_ratio",
                (specs.r_iv_ratio,),
                (_precise(_sums_exact(_iv_ratio_ball)),),
                _iv_ratio_rules(False),
            ),
        )
    ),
    ("iv_ratio_c", "d_d-d"): Entry(
        (
            Out(
                "iv_ratio_c",
                (specs.r_iv_ratio_c,),
                (_precise(_sums_exact(_iv_ratio_c_ball)),),
                _iv_ratio_rules(True),
            ),
        )
    ),
    ("wright_bessel", "d_d_d-d"): Entry(
        (
            Out(
                "wright_bessel",
                (_wright_bessel,),
                (_wright_bessel_arb,),
                (WRIGHT_DOMAIN,),
                WRIGHT_OVERRIDES,
                WRIGHT_FLOORS,
            ),
        )
    ),
    ("log_wright_bessel", "d_d_d-d"): Entry(
        (
            Out(
                "log_wright_bessel",
                (_log_wright_bessel,),
                (_log_wright_bessel_arb,),
                (
                    WRIGHT_DOMAIN,
                    Rule(
                        lambda a: (
                            a[1] == 0 and (a[0] == 0 or a[2] == 0) and min(a) >= 0
                        ),
                        -math.inf,
                        "Phi(a, 0; x) = 0 for a = 0 or x = 0",
                    ),
                ),
                LOG_WRIGHT_OVERRIDES,
            ),
        )
    ),
    ("pbdv", "d_d-d_d"): Entry(
        tuple(
            Out(
                f"pbdv.{i}",
                (_pair_ref(_pbdv, i),),
                (_pbdv_arb(i),),
                (NAN,),
                floors=PBDV_FLOORS,
            )
            for i in range(2)
        )
    ),
    ("pmv", "d_d_d-d"): Entry((Out("pmv", (_pmv_scaled,), (_pmv_arb,), PMV_RULES),)),
    ("pbwa", "d_d-d_d"): Entry(
        tuple(
            Out(
                f"pbwa.{i}",
                (_pair_ref(specs.r_pbwa, i),),
                (_pbwa_arb(i),),
                (NAN, PBWA_NAN),
            )
            for i in range(2)
        )
    ),
    ("hyp2f1", "d_d_d_d-d"): Entry(
        (
            Out(
                "hyp2f1",
                (_hyp2f1,),
                (_hyp2f1_arb,),
                HYP2F1_RULES,
                HYP2F1_OVERRIDES,
                HYP2F1_FLOORS,
            ),
        )
    ),
    **_spheroidal("prolate_aswfa", (SW_MODE, SW_ANGLE)),
    **_spheroidal("oblate_aswfa", (SW_MODE, SW_ANGLE)),
    **_spheroidal("prolate_radial1", (SW_PROLATE_X,), (True,)),
    **_spheroidal("prolate_radial2", (SW_PROLATE_X,)),
    **_spheroidal("oblate_radial1", (), (False,)),
    **_spheroidal("oblate_radial2", (), (False,)),
}
