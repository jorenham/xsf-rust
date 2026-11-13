"""The references, conventions, and tolerance exceptions of the tables that `generate.py`
replaces: an `Entry` per table and signature. The other tables keep those of xsref.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from typing import Final, Protocol

import flint
import mpmath as mp
import specs

type Pred = Callable[[tuple[float, ...]], bool]  # of the f64 inputs
type Value = mp.mpf | tuple[mp.mpf, mp.mpf]  # a reference, or (reference, scale)
type Arb = Callable[..., flint.arb | flint.acb | None]  # None where Arb can't check
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


def _bessel_v() -> dict[tuple[str, str], Entry]:
    """The tables of J_v, Y_v, e^-|x| I_v, and e^x K_v of real order v (J_v and Y_v are not
    scaled for real x)."""

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
    ) -> tuple[tuple[str, str], Entry]:
        func = f"bessel_{f}"
        return (f"cyl_{func}", "d_d-d"), Entry((Out(func, (ref,), (arb,), rules),))

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

    return dict(
        [
            entry("je", j, specs.ARB["bessel_j"], not_real("J_v")),
            entry("ye", specs.r_yv, specs.ARB["bessel_y"], cut("Y_v")),
            entry("ie", specs.r_ie, ie, not_real("I_v")),
            entry(
                "ke",
                specs.r_ke,
                ke,
                cut("K_v"),
                Rule(
                    lambda a: a[1] == 0, math.inf, "the limit for x -> 0+, at x = +-0"
                ),
            ),
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


def _precise[*A](f: Callable[[*A], flint.arb]) -> Callable[[*A], flint.arb | None]:
    """`f`, but None where its ball is too wide to confirm a reference."""

    def arb(*args: *A) -> flint.arb | None:
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
    "for a < 1, xsf-rust computes 1 / Gamma(a) with rgamma(a) instead of exp(-lgam(a)) "
    "(scipy/xsf#456)"
)


def _igam_fac_rgamma(a: float, x: float) -> bool:
    """Whether `igam_fac(a, x)` uses rgamma(a)."""
    return a < 1 and abs(a - x) > 0.4 * a


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
        if not 0 < a < 1:
            return False
        with mp.workdps(20):
            x = float(ref(mp.mpf(a), mp.mpf(y)))
        return _igam_fac_rgamma(a, x) or _uses_rgamma(a, x, upper != (y > 0.9))

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

    def ball(*a: flint.arb) -> flint.arb | flint.acb | None:
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

    def ball(*a: flint.arb) -> flint.arb | None:
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


def _voigt_arb(x: flint.arb, sigma: flint.arb, gamma: flint.arb) -> flint.arb:
    """Re w(z) / (sigma sqrt(2 pi)), with z = (x + i gamma) / (sigma sqrt(2)); the limits for
    sigma = 0 and gamma = 0 like `_voigt`."""
    if sigma == 0:
        return gamma / (flint.arb.pi() * (x * x + gamma * gamma))
    if gamma == 0:
        return (-((x / sigma) ** 2) / 2).exp() / (
            abs(sigma) * (2 * flint.arb.pi()).sqrt()
        )
    z = flint.acb(x, gamma) / (sigma * flint.arb(2).sqrt())
    w = (-z * z).exp() * (-flint.acb(0, 1) * z).erfc()
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
    ("gammaincinv", "d_d-d"): Entry((_gammaincinv_out(False, 8e-16),)),
    ("gammainccinv", "d_d-d"): Entry((_gammaincinv_out(True, 2e-14),)),
    ("chdtr", "d_d-d"): Entry(
        (_distribution("chdtr", False, _half, "(v / 2, x / 2)"),)
    ),
    ("chdtrc", "d_d-d"): Entry(
        (_distribution("chdtrc", True, _half, "(v / 2, x / 2)"),)
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
}
