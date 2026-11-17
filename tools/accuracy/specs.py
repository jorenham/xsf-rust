"""Functions, domains, and high-precision mpmath references for the accuracy table.

Each `Spec` names an xsf-rust function (the same name is used by the Rust and C++ dump
programs), its argument names, a reference function, and one or more `Domain`s to sample.
The components of functions that return a tuple or an array are separate specs, e.g.
`fresnel.0` and `airy[2]`.

A reference function takes `mpf` arguments (the exact values of the f64 inputs) and returns
either the exact function value, or a `(value, scale)` tuple. With a scale, the error is
`min(|y - f| / |f|, |y - f| / scale)`; this is used for oscillatory functions, where the
relative error is meaningless near the zeros (see `NOTES`).

Where mpmath is known to be inaccurate or slow, the references use alternative formulas, e.g.
a direct Euler-Maclaurin sum for the Hurwitz zeta function, quadrature for `expn`, Arb for
`gammainc` and `gammaincc`, and a continued fraction for `betainc`. The references of
inverse functions find the root of the reference of the function itself (see `_invert`).
"""

import math
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from random import Random
from types import FunctionType
from typing import Any, ClassVar, Final, Protocol

import flint
import mpmath as mp
from mpmath.libmp import NoConvergence

DBL_MAX: Final = 1.7976931348623157e308
DPS: Final = 50  # the working digits of every reference

# ---------------------------------------------------------------- sampling

# Unicode math for the table, escaped to keep the source ASCII
IN, MINUS, PM, ELLIPSIS, SQRT = "\u2208", "\u2212", "\u00b1", "\u2026", "\u221a"
HALF, MU, SIGMA, GAMMA, EPSILON, PSI = (
    "\u00bd",
    "\u03bc",
    "\u03c3",
    "\u0393",
    "\u03b5",
    "\u03c8",
)
SUP = str.maketrans(
    "-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079"
)
GREEK = {
    "phi": "\u03c6",
    "lambda": "\u03bb",
    "delta": "\u03b4",
    "sigma": "\u03c3",
    "gamma": "\u03b3",
    "theta": "\u03b8",
    "nu": "\u03bd",
    "alpha": "\u03b1",
    "beta": "\u03b2",
    "re": "Re\u00a0z",
    "im": "Im\u00a0z",
    "a1": "a\u2081",
    "a2": "a\u2082",
    "a3": "a\u2083",
    "a4": "a\u2084",
    "a5": "a\u2085",
    "b1": "b\u2081",
    "b2": "b\u2082",
    "b3": "b\u2083",
    "b4": "b\u2084",
    "xr": "Re\u00a0x",
    "xi": "Im\u00a0x",
    "yr": "Re\u00a0y",
    "yi": "Im\u00a0y",
}


def nobreak(*parts: str) -> str:
    """Join the parts of one argument's domain with no-break spaces, and keep a leading sign
    attached to its interval, so that table cells only wrap between the arguments, or after
    a comma within the interval (the last part)."""
    *head, tail = (p.replace(" ", "\u00a0") for p in parts)
    s = "\u00a0".join([*head, tail.replace(",\u00a0", ", ")])
    return s.replace(f"{PM}[", f"{PM}\u2060[").replace(f"{MINUS}[", f"{MINUS}\u2060[")


def _var(argname: str) -> str:
    """The argument of a real or imaginary part, e.g. "z" for "re" (Re z), or "" if real."""
    name = GREEK.get(argname, argname)
    return name[3:] if name.startswith(("Re\u00a0", "Im\u00a0")) else ""


def wrap(s: str) -> str:
    """Put each argument's domain on its own line, i.e. break at the ", " between them (the ones
    outside brackets; within an argument's domain, `nobreak` leaves only those in intervals)."""
    parts: list[str] = []
    depth = start = 0
    for i, c in enumerate(s):
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == "," and depth == 0 and s[i + 1 : i + 2] == " ":
            parts.append(s[start:i])
            start = i + 2
    return "<br>".join([*parts, s[start:]])


class Dist(Protocol):
    """The distribution of one argument."""

    def sample(self, rng: Random) -> float: ...
    def desc(self) -> str: ...
    def pretty(self) -> str: ...


@dataclass(frozen=True)
class U:
    """Uniform on [lo, hi]."""

    lo: float
    hi: float

    def sample(self, rng: Random) -> float:
        return rng.uniform(self.lo, self.hi)

    def desc(self) -> str:
        return f"[{_g(self.lo)}, {_g(self.hi)}]"

    def pretty(self) -> str:
        return f"[{_u(self.lo)}, {_u(self.hi)}]"


@dataclass(frozen=True)
class LU:
    """Log-uniform on [lo, hi] (0 < lo < hi), optionally with a random sign."""

    lo: float
    hi: float
    sign: int = 1  # 1, -1, or 0 for a random sign

    def sample(self, rng: Random) -> float:
        x = 10 ** rng.uniform(math.log10(self.lo), math.log10(self.hi))
        s = self.sign or rng.choice((-1, 1))
        return s * min(x, self.hi)

    def desc(self) -> str:
        pre = {1: "", -1: "-", 0: "+-"}[self.sign]
        return f"{pre}[{_g(self.lo)}, {_g(self.hi)}] (log)"

    def pretty(self) -> str:
        pre = {1: "", -1: MINUS, 0: PM}[self.sign]
        return f"{pre}[{_u(self.lo)}, {_u(self.hi)}]"


@dataclass(frozen=True)
class Int:
    """Uniform integer in [lo, hi]."""

    lo: int
    hi: int

    def sample(self, rng: Random) -> float:
        return float(rng.randint(self.lo, self.hi))

    def desc(self) -> str:
        return f"{{{self.lo}, ..., {self.hi}}}"

    def pretty(self) -> str:
        return f"{{{self.lo}, {ELLIPSIS}, {self.hi}}}".replace("-", MINUS)


@dataclass(frozen=True)
class ILU:
    """Log-uniform integer in [lo, hi] (1 <= lo < hi)."""

    lo: int
    hi: int

    def sample(self, rng: Random) -> float:
        return float(min(int(LU(self.lo, self.hi + 1).sample(rng)), self.hi))

    def desc(self) -> str:
        return f"{{{_g(self.lo)}, ..., {_g(self.hi)}}} (log)"

    def pretty(self) -> str:
        return f"{{{_u(self.lo)}, {ELLIPSIS}, {_u(self.hi)}}}"


@dataclass(frozen=True)
class Each:
    """Each integer in [lo, hi] in turn, from lo again for each new `rng` (one per sampling
    run), so that `hi - lo + 1` samples have each one once."""

    lo: int
    hi: int
    # the `rng` of the current run, and the number of samples drawn with it
    _run: list[Any] = field(
        default_factory=lambda: [None, 0], compare=False, repr=False
    )

    def sample(self, rng: Random) -> float:
        if self._run[0] is not rng:
            self._run[:] = [rng, 0]
        i = self._run[1]
        self._run[1] += 1
        return float(self.lo + i % (self.hi - self.lo + 1))

    def desc(self) -> str:
        return f"each of {{{self.lo}, ..., {self.hi}}}"

    def pretty(self) -> str:
        return Int(self.lo, self.hi).pretty()


@dataclass(frozen=True)
class Offset:
    """c + sign * t (rounded to f64), where t is sampled from `dist`."""

    c: float
    dist: Dist
    sign: int = 1

    def sample(self, rng: Random) -> float:
        return self.c + self.sign * self.dist.sample(rng)

    def desc(self) -> str:
        return f"{_g(self.c)} {'+' if self.sign > 0 else '-'} {self.dist.desc()}"

    def pretty(self) -> str:
        return f"{_u(self.c)} {'+' if self.sign > 0 else MINUS} {self.dist.pretty()}"


def _g(x: float) -> str:
    """Short repr of a domain bound."""
    if x == int(x) and abs(x) < 1e6:
        return str(int(x))
    s = f"{x:.4g}"
    return s.replace("e+", "e").replace("e-0", "e-").replace("e0", "e")


def _u(x: float) -> str:
    """Unicode repr of a domain bound, e.g. 10^-3 with a superscript exponent."""
    sign = MINUS if x < 0 else ""
    a = abs(x)
    named = {
        math.pi: "\u03c0",
        2 * math.pi: "2\u03c0",
        math.pi / 2: "\u03c0/2",
        -INV_E: "1/e",
    }
    if a in named:
        return sign + named[a]
    e = round(math.log10(a)) if a else 0
    if (e <= -3 or e >= 4) and a == float(f"1e{e}"):
        return f"{sign}10{str(e).translate(SUP)}"
    return sign + _g(a)


def _is_log(d: Dist) -> bool:
    return isinstance(d, (LU, ILU)) or (isinstance(d, Offset) and _is_log(d.dist))


class Domain:
    """Independent per-argument distributions, or a custom `sampler(rng) -> args`.

    `call` overrides the function name passed to the dump programs (e.g. for an integer-degree
    variant), and `cpp=False` marks a variant that xsf does not have.
    """

    def __init__(
        self,
        *dists: Dist,
        n: int = 10000,
        sampler: Callable[[Random], Iterable[float]] | None = None,
        label: str = "",
        pretty: str = "",
        call: str = "",
        cpp: bool = True,
    ) -> None:
        self.dists, self.n, self.sampler = dists, n, sampler
        self.label, self.pretty_label = label, pretty
        self.call, self.cpp = call, cpp

    def desc(self, argnames: Iterable[str]) -> str:
        """Plain description; also seeds the sampling."""
        if self.label:
            return self.label
        return ", ".join(
            f"{a} in {d.desc()}" for a, d in zip(argnames, self.dists, strict=True)
        )

    def pretty(self, argnames: Iterable[str]) -> str:
        """Unicode description for the table; consecutive equal distributions are merged, but
        the real and imaginary parts of different arguments are not, to keep the lines short."""
        if self.pretty_label:
            return wrap(self.pretty_label)
        groups: list[tuple[list[str], Dist]] = []
        for a, d in zip(argnames, self.dists, strict=True):
            if groups and groups[-1][1] == d and _var(groups[-1][0][-1]) == _var(a):
                groups[-1][0].append(a)
            else:
                groups.append(([a], d))
        return wrap(
            ", ".join(
                nobreak(", ".join(GREEK.get(a, a) for a in names), IN, d.pretty())
                + (" (log)" if _is_log(d) else "")
                for names, d in groups
            )
        )

    def sample(self, rng: Random) -> tuple[float, ...]:
        if self.sampler is not None:
            return tuple(float(a) for a in self.sampler(rng))
        return tuple(d.sample(rng) for d in self.dists)


D = Domain


@dataclass
class Spec:
    name: str  # a suffix "#..." distinguishes specs of the same function, e.g. for complex input
    args: tuple[str, ...]
    ref: FunctionType  # its source is part of the cache key of the references
    domains: list[Domain]
    section: str
    std: bool = False  # has a Rust std equivalent
    cpp: bool = True  # has an xsf (C++) equivalent
    complex: bool = False  # complex-valued
    timeout: int = 30
    note: str = ""  # space-separated keys into NOTES
    rust: ClassVar[bool] = True

    @property
    def func(self) -> str:
        """The function name, without the "#..." suffix."""
        return self.name.split("#")[0]


# ---------------------------------------------------------------- references


def wdps(extra):
    return mp.workdps(mp.mp.dps + extra)


def xd(*a):
    """Extra digits needed for arguments of magnitude |a|."""
    m = max([abs(float(x)) for x in a] + [1.0])
    return int(math.log10(m)) + 5 if m > 1 else 5


def _envelope(v, x):
    """Bessel modulus sqrt(J_v^2 + Y_v^2), the local amplitude of J_v and Y_v."""
    return mp.sqrt(mp.besselj(v, x) ** 2 + mp.bessely(v, x) ** 2)


def _osc(f, v, x):
    """Use the Bessel modulus as error scale in the oscillatory region x > |v|."""
    return (f, _envelope(v, x)) if x > abs(v) else f


def r_gamma(x):
    with wdps(xd(x)):
        return mp.gamma(x)


def r_gammaln(x):
    with wdps(xd(x)):
        return mp.log(abs(mp.gamma(x))) if abs(x) < 1e5 else mp.re(mp.loggamma(x))


def r_rgamma(x):
    with wdps(xd(x)):
        return mp.rgamma(x)


def r_digamma(x):
    with wdps(xd(x)):
        return mp.digamma(x)


def r_trigamma(x):
    with wdps(xd(x)):
        return mp.psi(1, x)


def r_digammainv(y):
    # Newton on digamma(x) = y, from Minka's initial guess
    with wdps(20):
        x = mp.exp(y) + mp.mpf(1) / 2 if y >= -2.22 else -1 / (y + mp.euler)
        for _ in range(200):
            dx = (mp.digamma(x) - y) / mp.psi(1, x)
            x -= dx
            if abs(dx) < abs(x) * mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        return x


def r_beta(a, b):
    with wdps(xd(a, b)):
        return mp.beta(a, b)


def r_betaln(a, b):
    with wdps(xd(a, b)):
        return mp.log(abs(mp.beta(a, b)))


def r_binom(n, k):
    with wdps(xd(n, k)):
        return mp.binomial(n, k)


def r_poch(a, m):
    with wdps(xd(a, m)):
        return mp.rf(a, m)


def r_erf(x):
    return mp.erf(x)


def r_erfc(x):
    return mp.erfc(x)


def r_erfcx(x):
    with wdps(xd(x)):
        return mp.exp(x * x) * mp.erfc(x)


def _ncdf_inv(p):
    """x such that ncdf(x) = p, for 0 < p < 1; Newton in log space."""
    with wdps(20):
        x = mp.sqrt(2) * mp.erfinv(2 * p - 1) if p > 1e-20 else -mp.sqrt(-2 * mp.log(p))
        for _ in range(200):
            c = mp.ncdf(x)
            dx = (mp.log(c) - mp.log(p)) * c / mp.npdf(x)
            x -= dx
            if abs(dx) <= abs(x) * mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        return x


def r_erfinv(y):
    # erfinv(y) = ndtri((1 + y) / 2) / sqrt(2); 1 + y is exact at this precision
    with mp.workdps(400):
        if y < 0:
            return -_ncdf_inv((1 - y) / 2) / mp.sqrt(2)
        return _ncdf_inv((1 + y) / 2) / mp.sqrt(2)


def r_erfcinv(y):
    # erfcinv(y) = -ndtri(y / 2) / sqrt(2)
    with mp.workdps(400):
        if y > 1:
            return _ncdf_inv(1 - y / 2) / mp.sqrt(2)
        return -_ncdf_inv(y / 2) / mp.sqrt(2)


def r_ndtr(x):
    return mp.ncdf(x)


def r_ndtri(p):
    with mp.workdps(400):
        return _ncdf_inv(p)


def r_log_ndtr(x):
    with wdps(xd(x)):
        if x > 0:
            return mp.log1p(-mp.ncdf(-x))
        return mp.log(mp.ncdf(x))


def r_expm1(x):
    return mp.expm1(x)


def r_log1p(x):
    return mp.log1p(x)


def r_cbrt(x):
    return mp.cbrt(x) if x >= 0 else -mp.cbrt(-x)


def r_exp2(x):
    with wdps(5):
        return mp.power(2, x)


def r_exp10(x):
    with wdps(5):
        return mp.power(10, x)


def r_sinpi(x):
    return mp.sinpi(x)


def r_cospi(x):
    return mp.cospi(x)


def r_cosine_cdf(x):
    return (mp.pi + x + mp.sin(x)) / (2 * mp.pi)


def r_cosine_invcdf(p):
    # x = t - pi, where g(t) = t - sin(t) = 2 pi p; Newton on the convex g for t in [0, pi]
    if p > 0.5:
        return -r_cosine_invcdf(1 - p)
    with wdps(20):
        y = 2 * mp.pi * p
        eps = mp.mpf(10) ** (-mp.mp.dps + 5)

        def g(t):
            if t > 0.5:
                return t - mp.sin(t)
            # series, to avoid the cancellation for small t
            s, term, k = mp.mpf(0), t**3 / 6, 1
            while abs(term) > eps * s:
                s += term
                term *= -t * t / ((2 * k + 2) * (2 * k + 3))
                k += 1
            return s

        t = min(mp.cbrt(6 * y), mp.pi)
        for _ in range(1000):
            dt = (g(t) - y) / (2 * mp.sin(t / 2) ** 2)
            t -= dt
            if abs(dt) <= eps * t:
                break
        return t - mp.pi


def r_j0(x):
    return _osc(mp.besselj(0, x), 0, x)


def r_j1(x):
    return _osc(mp.besselj(1, x), 1, x)


def r_y0(x):
    return _osc(mp.bessely(0, x), 0, x)


def r_y1(x):
    return _osc(mp.bessely(1, x), 1, x)


def r_jv(v, x):
    with wdps(xd(v, x)):
        return _osc(mp.besselj(v, x), v, x)


def r_yv(v, x):
    with wdps(xd(v, x)):
        return _osc(mp.bessely(v, x), v, x)


def r_iv(v, x):
    with wdps(xd(v, x)):
        return mp.besseli(v, x)


def r_kv(v, x):
    with wdps(xd(v, x)):
        return mp.besselk(v, x)


def r_i0(x):
    return mp.besseli(0, x)


def r_i1(x):
    return mp.besseli(1, x)


def r_k0(x):
    return mp.besselk(0, x)


def r_k1(x):
    return mp.besselk(1, x)


def r_struve_h(v, x):
    with wdps(xd(v, x) + 20):
        return _osc(mp.struveh(v, x), v, x)


def r_struve_l(v, x):
    with wdps(xd(v, x) + 20):
        return mp.struvel(v, x)


def r_hurwitz_zeta(x, q):
    # direct Euler-Maclaurin; mpmath's zeta(s, a) is inaccurate for large a at dps < ~80
    with wdps(xd(x, q) + 30):
        N = 60
        s = mp.fsum((q + k) ** -x for k in range(N))
        w = q + N
        s += w ** (1 - x) / (x - 1) + w**-x / 2
        fac = x
        t = w ** (-x - 1)
        for k in range(1, 25):
            s += mp.bernoulli(2 * k) / mp.factorial(2 * k) * fac * t
            fac *= (x + 2 * k - 1) * (x + 2 * k)
            t /= w * w
        return s


def r_riemann_zeta(x):
    with wdps(xd(x)):
        return mp.zeta(x)


def r_zetac(x):
    # zeta(x) - 1 ~ 2^-x for large x: add enough digits to resolve it
    with wdps(xd(x) + (int(x * 0.302) if x > 0 else 0)):
        return mp.zeta(x) - 1


def r_expi(x):
    with wdps(xd(x)):
        return mp.ei(x)


def r_expn(n, x):
    # quadrature: mpmath's expint(n, x) is badly wrong for n ~ x >~ 100
    with wdps(10):
        if n <= 1 or x < 1:
            return mp.expint(int(n), x)
        s = x + n
        pts = [1, 1 + 1 / s, 1 + 10 / s, 1 + 100 / s, 1 + 1000 / s, mp.inf]
        return mp.exp(-x) * mp.quad(lambda t: mp.exp(-x * (t - 1)) * t**-n, pts)


def _gammainc_series(a, x):
    """P(a, x) via the power series (DLMF 8.11.4)."""
    with wdps(xd(a, x) * 2):
        t = s = mp.mpf(1)
        k = 0
        while t > mp.eps * s:
            k += 1
            t *= x / (a + k)
            s += t
        return mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1)) * s


def _gammainc_pq(a, x, upper):
    """P(a, x) or Q(a, x), to the mpmath precision: Q with Arb, and P as 1 - Q, or where P is
    too small for that, with the series (Arb's P is wrong for large a and x < a)."""
    # 1 - Q for P down to 2^(-15 mp.prec), i.e. about 1e-760 at 50 digits
    for k in (2, 4, 8, 16):
        with flint.ctx.workprec(k * mp.mp.prec):
            q = flint.arb(x).gamma_upper(flint.arb(a), regularized=1)
            y = q if upper else 1 - q
            if y.rel_accuracy_bits() >= mp.mp.prec:
                return mp.mpf(y.mid())
    assert not upper, (a, x)
    return _gammainc_series(a, x)


def r_gammainc(a, x):
    return _gammainc_pq(a, x, False)


def r_gammaincc(a, x):
    return _gammainc_pq(a, x, True)


def _log_gammainc_pq(a, x, upper):
    """ln P(a, x) or ln Q(a, x), as the log1p of minus the complement where that is smaller."""
    y = _gammainc_pq(a, x, upper)
    return mp.log1p(-_gammainc_pq(a, x, not upper)) if y > 0.5 else mp.log(y)


def r_log_gammainc(a, x):
    return _log_gammainc_pq(a, x, False)


def r_log_gammaincc(a, x):
    return _log_gammainc_pq(a, x, True)


def _cf_betainc(a, b, x, y):
    """I_x(a, b), with y = 1 - x, via the modified Lentz continued fraction (DLMF 8.17.22)."""
    if x == 0:
        return mp.mpf(0)
    lbeta = mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
    lpre = a * mp.log(x) + b * mp.log(y) - mp.log(a) - lbeta
    tiny = mp.mpf(10) ** (-2 * mp.mp.dps)
    eps = mp.mpf(10) ** (-mp.mp.dps + 3)

    def clamp(t):
        return tiny if abs(t) < tiny else t

    c = mp.mpf(1)
    d = 1 / clamp(1 - (a + b) * x / (a + 1))
    h = d
    m = 0
    while True:
        m += 1
        m2 = 2 * m
        aa = m * (b - m) * x / ((a + m2 - 1) * (a + m2))
        d = 1 / clamp(1 + aa * d)
        c = clamp(1 + aa / c)
        h *= d * c
        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1))
        d = 1 / clamp(1 + aa * d)
        c = clamp(1 + aa / c)
        de = d * c
        h *= de
        if abs(de - 1) < eps or m > 10**7:
            break
    return mp.exp(lpre) * h


def _betainc(a, b, x, y):
    """I_x(a, b), given y = 1 - x, so that x and y can both be computed without cancellation."""
    with wdps(xd(a, b) * 2 + 10):
        if x < (a + 1) / (a + b + 2):
            return _cf_betainc(a, b, x, y)
        return 1 - _cf_betainc(b, a, y, x)


def r_betainc(a, b, x):
    return _betainc(a, b, x, 1 - x)


def r_betaincc(a, b, x):
    # 1 - I_x(a, b) = I_{1-x}(b, a) (DLMF 8.17.4), so the CF is of the smaller of the two
    return _betainc(b, a, 1 - x, x)


def r_fdtr(a, b, x):
    return _betainc(a / 2, b / 2, a * x / (b + a * x), b / (b + a * x))


def r_fdtrc(a, b, x):
    return _betainc(b / 2, a / 2, b / (b + a * x), a * x / (b + a * x))


def r_ellipk(m):
    with wdps(xd(m)):
        return mp.ellipk(m)


def r_ellipe(m):
    with wdps(xd(m)):
        return mp.ellipe(m)


def r_ellipkinc(phi, m):
    with wdps(xd(phi, m)):
        return mp.ellipf(phi, m)


def r_ellipeinc(phi, m):
    with wdps(xd(phi, m)):
        return mp.ellipe(phi, m)


def _lame(h2, k2, n, p):
    """E^p_n as (e0, e1, e2, R): E(s) = s^e0 |s^2 - h2|^(e1/2) |s^2 - k2|^(e2/2) R(s^2), the
    type (K, L, M, N, in SciPy's order of p, see DLMF 29.12.1) and the monic polynomial R.

    Not SciPy's method: in y = s^2, with A = y - h2 and B = y - k2, the Lame equation for E is
    T R = -a R, with T = 4 y A B R'' + Q1 R' + Q0 R, where
    Q1 = 2 (1 + 2 e0) A B + 2 y ((1 + 2 e1) B + (1 + 2 e2) A) and
    Q0 = (e0 + e1 + 2 e0 e1) B + (e0 + e2 + 2 e0 e2) A + (e1 + e2 + 2 e1 e2 - n (n + 1)) y.
    T maps the polynomials of degree m = (n - e0 - e1 - e2) / 2 to themselves, so R is an
    eigenvector of its matrix; p counts up the eigenvalues a within each type."""
    n, p = int(n), int(p)
    r = n // 2
    for e1, e2, size in [(0, 0, r + 1), (1, 0, n - r), (0, 1, n - r), (1, 1, r)]:
        if p <= size:
            break
        p -= size
    e0 = (n - e1 - e2) % 2
    m = (n - e0 - e1 - e2) // 2
    prod, total = h2 * k2, h2 + k2
    # the coefficients (ascending in y) of R, R', and R'' in T R
    q0 = [
        -(e0 + e1 + 2 * e0 * e1) * k2 - (e0 + e2 + 2 * e0 * e2) * h2,
        2 * (e0 + e1 + e2 + e0 * e1 + e0 * e2 + e1 * e2) - n * (n + 1),
    ]
    q1 = [
        2 * (1 + 2 * e0) * prod,
        -2 * (1 + 2 * e0) * total - 2 * (1 + 2 * e1) * k2 - 2 * (1 + 2 * e2) * h2,
        2 * (1 + 2 * e0) + 4 * (1 + e1 + e2),
    ]
    q2 = [0, 4 * prod, -4 * total, 4]
    t = mp.zeros(m + 1, m + 1)
    for j in range(m + 1):
        # T y^j = Q0 y^j + j Q1 y^(j - 1) + j (j - 1) 4 y A B y^(j - 2), without the y^(m + 1)
        # term, which vanishes for j = m
        for d, q in enumerate([q0, q1, q2]):
            for i, c in enumerate(q):
                if 0 <= i + j - d <= m:
                    t[i + j - d, j] += mp.ff(j, d) * c
    ev, er = mp.eig(t)
    # the p-th smallest a = -ev
    j = sorted(range(m + 1), key=lambda i: -mp.re(ev[i]))[p - 1]
    return e0, e1, e2, [mp.re(er[i, j] / er[m, j]) for i in range(m + 1)]


def _lame_eval(lame, s, a, b):
    """E(s) of `_lame` with the prefactor signs +1, given a = s^2 - h2 and b = s^2 - k2."""
    e0, e1, e2, r = lame
    psi = s**e0 * mp.sqrt(abs(a)) ** e1 * mp.sqrt(abs(b)) ** e2
    return psi * mp.polyval(r[::-1], s * s)


def r_ellip_harm(h2, k2, n, p, s, signm, signn):
    """E^p_n(s) of `_lame`, with the scale max |E| on [0, h] or [h, k] (whichever has s) if s < k;
    its zeros are there."""
    lame = _lame(h2, k2, n, p)
    sign = signm ** lame[1] * signn ** lame[2]
    e = sign * _lame_eval(lame, s, s * s - h2, s * s - k2)
    if s * s >= k2:
        return e
    # a grid of t = h sin(phi) or t^2 = h2 cos^2(theta) + k2 sin^2(theta), denser at the ends
    xs = [mp.pi / 2 * i / 128 for i in range(129)]
    if s * s < h2:
        ts = [(mp.sqrt(h2) * mp.sin(x), -h2 * mp.cos(x) ** 2) for x in xs]
        pts = [(t, a, t * t - k2) for t, a in ts]
    else:
        d = k2 - h2
        pts = [
            (mp.sqrt(h2 + d * mp.sin(x) ** 2), d * mp.sin(x) ** 2, -d * mp.cos(x) ** 2)
            for x in xs
        ]
    return e, max(abs(_lame_eval(lame, *pt)) for pt in pts)


def _quad_rel(f, pts):
    """`mp.quad` to a relative tolerance: its tolerance is absolute, so f is first scaled by a
    rough integral."""
    with mp.workdps(15):
        c = mp.fabs(mp.quad(f, pts))
    return c * mp.quad(lambda x: f(x) / c, pts)


def r_ellip_harm_2(h2, k2, n, p, s):
    """(2n + 1) E(s) times the integral over t in [s, inf) of 1 / (E(t)^2 sqrt((t^2 - h2)
    (t^2 - k2))), with t^2 = k2 + v^2 and v = v_s e^x, which leaves v / (t E(t)^2 sqrt(t^2 - h2))
    dx over x in [0, inf): smooth, also for s ~ k (v_s ~ 0)."""
    lame = _lame(h2, k2, n, p)
    vs = mp.sqrt(s * s - k2)

    def f(x):
        v = vs * mp.exp(x)
        t = mp.sqrt(k2 + v * v)
        a = k2 - h2 + v * v
        return v / (t * _lame_eval(lame, t, a, v * v) ** 2 * mp.sqrt(a))

    # v ~ sqrt(k2) is where the integrand changes from ~1 to ~v^(-2n - 1)
    xk = mp.log(mp.sqrt(k2) / vs)
    pts = [0, xk, xk + 3, mp.inf] if xk > 0 else [0, 3, mp.inf]
    e = _lame_eval(lame, s, s * s - h2, vs * vs)
    return (2 * n + 1) * e * _quad_rel(f, pts)


def r_ellip_normal(h2, k2, n, p):
    """8 (J2 I0 - J0 I2), the double integral split by y^2 - x^2, with the integrals Ji (of
    y^i E(y)^2 / sqrt((k2 - y^2)(y^2 - h2)) over [h, k]) and Ii (of x^i E(x)^2 /
    sqrt((h2 - x^2)(k2 - x^2)) over [0, h]); with y^2 = h2 cos^2(theta) + k2 sin^2(theta) and
    x = h sin(phi), so that their integrands are smooth."""
    lame = _lame(h2, k2, n, p)
    d = k2 - h2

    def j(i):
        def f(theta):
            a, b = d * mp.sin(theta) ** 2, -d * mp.cos(theta) ** 2
            y2 = h2 + a
            return y2 ** (i / 2 - 1 / 2) * _lame_eval(lame, mp.sqrt(y2), a, b) ** 2

        return _quad_rel(f, [0, mp.pi / 4, mp.pi / 2])

    def i_(i):
        def f(phi):
            x = mp.sqrt(h2) * mp.sin(phi)
            a, b = -h2 * mp.cos(phi) ** 2, x * x - k2
            return x**i * _lame_eval(lame, x, a, b) ** 2 / mp.sqrt(-b)

        return _quad_rel(f, [0, mp.pi / 4, mp.pi / 2])

    return 8 * (j(2) * i_(0) - j(0) * i_(2))


def _arb(name, *args):
    """`ARB[name]` with increasing precision, until it is accurate to the working precision."""
    prec = mp.mp.prec
    while prec <= 2**14:
        with flint.ctx.workprec(prec):
            y = ARB[name](*map(flint.arb, args))
        if y.rel_accuracy_bits() >= mp.mp.prec:
            return mp.mpf(y)
        prec *= 2
    raise NoConvergence


def _hyp0f1_series(b, z):
    """The power series of hyp0f1 (in the working precision), and its largest term.

    For k > -b, the term ratios z / ((b + k) (k + 1)) decrease monotonically in absolute value,
    so they are at most rho (below), which bounds the rest geometrically."""
    s, t, top, k = mp.mpf(0), mp.mpf(1), mp.mpf(1), 0
    while True:
        s += t
        top = max(top, abs(t))
        t *= z / ((b + k) * (k + 1))
        k += 1
        if k > -b:
            rho = abs(z) / ((b + k) * (k + 1))
            if rho < 1 and abs(t) * rho / (1 - rho) <= mp.eps * abs(s):
                return s + t, top


def r_hyp0f1(b, z):
    # mpmath's hyp0f1 can be far off for b < 0 < z (mpmath/mpmath#1216), so for -10^4 < b < 0 < z
    # sum the series directly, with 20 more digits than it loses to cancellation, up to 300
    # digits (else use Arb). For b <= -10^4 (where mpmath agreed with the series on all samples)
    # that takes too many terms.
    if -(10**4) < b < 0 < z:
        extra = 20
        while True:
            with wdps(extra):
                s, top = _hyp0f1_series(b, z)
                if s != 0 and top / abs(s) < mp.mpf(10) ** (extra - 15):
                    return s
            extra = int(mp.log10(top / abs(s))) + 20 if s != 0 else 2 * extra
            if mp.mp.dps + extra > 300:
                return _arb("hyp0f1", b, z)
    # for z < 0, 0F1(; b; z) = Gamma(b) (-z)^(-v/2) J_v(x), with v = b - 1 and x = 2 sqrt(-z)
    with wdps(xd(b, z)):
        f = mp.hyp0f1(b, z)
        if z < 0 and -4 * z > (b - 1) ** 2:
            v, x = b - 1, 2 * mp.sqrt(-z)
            return f, abs(mp.gamma(b)) * (-z) ** (-v / 2) * _envelope(v, x)
        return f


def r_lambertw(x):
    with wdps(xd(x)):
        return mp.lambertw(x)


def r_spence(x):
    # spence(x) = Li2(1 - x); for x > 1 use Li2(1 - x) = -Li2(1 - 1/x) - ln(x)^2 / 2
    with mp.workdps(mp.mp.dps + 10 + max(0, int(mp.log10(x)))):
        if x <= 1:
            return mp.polylog(2, 1 - x)
        return -mp.polylog(2, 1 - 1 / x) - mp.log(x) ** 2 / 2


# ---------------------------------------------------------------- root finding


def _solve(G, a, b, ga, gb):
    """Root of G in [a, b], where G(a) = ga and G(b) = gb have opposite signs (the Illinois
    method, with a bisection step every fourth iteration to guarantee progress)."""
    eps = mp.mpf(10) ** (-mp.mp.dps + 5)
    tiny = mp.mpf(10) ** (-3 * mp.mp.dps)
    side = 0
    for i in range(10000):
        if ga == 0:
            return a
        if gb == 0:
            return b
        if abs(b - a) <= eps * max(abs(a), abs(b)) + tiny:
            break
        c = (a + b) / 2 if i % 4 == 3 else (a * gb - b * ga) / (gb - ga)
        gc = G(c)
        if (gc < 0) == (ga < 0):
            a, ga = c, gc
            if side == -1:
                gb /= 2
            side = -1
        else:
            b, gb = c, gc
            if side == 1:
                ga /= 2
            side = 1
    return (a * gb - b * ga) / (gb - ga)


def _invert(G, t0, step=1):
    """Root of the increasing function G, starting the search for a bracket at t0."""
    with wdps(10):
        g0 = G(t0)
        if g0 == 0:
            return t0
        d = step if g0 < 0 else -step
        for _ in range(2000):
            t1 = t0 + d
            g1 = G(t1)
            if (g1 < 0) != (g0 < 0):
                return _solve(G, t0, t1, g0, g1) if d > 0 else _solve(G, t1, t0, g1, g0)
            t0, g0, d = t1, g1, 2 * d
        raise ArithmeticError("no bracket")


def _invert_pos(F, target, x0, decreasing=False):
    """x > 0 such that F(x) = target, for a monotonic F; the root is found in log space."""
    sign = -1 if decreasing else 1
    return mp.exp(_invert(lambda t: sign * (F(mp.exp(t)) - target), mp.log(x0)))


def _invert_unit(F, target, x0, decreasing=False):
    """x in (0, 1) such that F(x, 1 - x) = target, for a monotonic F; the root is found in logit
    space, and F gets both x and 1 - x, so that neither loses precision near 0 or 1."""
    sign = -1 if decreasing else 1

    def G(s):
        x, y = 1 / (1 + mp.exp(-s)), 1 / (1 + mp.exp(s))
        return sign * (F(x, y) - target)

    s = _invert(G, mp.log(x0) - mp.log1p(-x0))
    return 1 / (1 + mp.exp(-s))


def _only_root(G, lo, hi, n):
    """Root of G in [lo, hi], which must be the only sign change of G (at 30 fewer digits) on a
    grid of n + 1 points, so that the root is bracketed and, up to the grid, unique."""
    ts = [lo + (hi - lo) * k / n for k in range(n + 1)]
    with wdps(-30):
        neg = [G(t) < 0 for t in ts]
    ks = [k for k in range(n) if neg[k] != neg[k + 1]]
    if len(ks) != 1:
        raise ArithmeticError(f"{len(ks)} sign changes")
    a, b = ts[ks[0]], ts[ks[0] + 1]
    with wdps(10):
        ga, gb = G(a), G(b)
        if (ga < 0) == (gb < 0):
            raise ArithmeticError("no bracket")
        return _solve(G, a, b, ga, gb)


def _only_root_pos(F, target):
    """The only x in [1e-100, 1e12] with F(x) = target (see `_only_root`), on a grid in log space
    with two points per decade."""
    lo, hi = -100 * mp.ln10, 12 * mp.ln10
    return mp.exp(_only_root(lambda t: F(mp.exp(t)) - target, lo, hi, 224))


# ---------------------------------------------------------------- gamma and related


def r_loggamma(x):
    with wdps(xd(x)):
        return mp.loggamma(x)


def r_gammasgn(x):
    if x > 0:
        return mp.mpf(1)
    return mp.mpf(-1) if int(mp.floor(-x)) % 2 == 0 else mp.mpf(1)


def r_lgam1p(x):
    # log(Gamma(1 + x)), with the Taylor series around 0, since 1 + x is inexact for tiny x
    if abs(x) < 0.5:
        with wdps(10):
            s, k, xk = -mp.euler * x, 2, x * x
            while True:
                t = (-1) ** k * mp.zeta(k) * xk / k
                s += t
                if abs(t) < mp.eps * abs(s):
                    return s
                k, xk = k + 1, xk * x
    with wdps(xd(x) + 30):
        return mp.loggamma(1 + x)


def r_polygamma(n, x):
    with wdps(xd(x) + 10):
        return mp.polygamma(int(n), x)


def r_multigammaln(a, d):
    d = int(d)
    with wdps(xd(a) + 10):
        return d * (d - 1) / 4 * mp.log(mp.pi) + mp.fsum(
            mp.loggamma(a - mp.mpf(j) / 2) for j in range(d)
        )


def r_ff(x, m):
    with wdps(xd(x, m)):
        return mp.ff(x, m)


def r_stirling2(n, k):
    # exact, with integer arithmetic
    n, k = int(n), int(k)
    s = sum((-1) ** i * math.comb(k, i) * (k - i) ** n for i in range(k + 1))
    return mp.mpf(s // math.factorial(k))


# ---------------------------------------------------------------- error function and normal


def r_erfi(x):
    with wdps(xd(x)):
        return mp.erfi(x)


def r_dawsn(x):
    if abs(x) > 100:
        # asymptotic series 1/(2x) sum (2k - 1)!! / (2x^2)^k; the terms decrease until k ~ x^2
        s, t, k = mp.mpf(0), 1 / (2 * x), 0
        while abs(t) > mp.eps * abs(s):
            s += t
            k += 1
            t *= (2 * k - 1) / (2 * x * x)
        return s
    with wdps(10):
        return mp.sqrt(mp.pi) / 2 * mp.exp(-x * x) * mp.erfi(x)


def r_voigt_profile(x, sigma, gamma):
    # Re w(z) / (sigma sqrt(2 pi)), with z = (x + i gamma) / (sigma sqrt(2)) and the Faddeeva
    # function w(z) = exp(-z^2) erfc(-iz)
    with wdps(20):
        z = mp.mpc(x, gamma) / (sigma * mp.sqrt(2))
        return mp.re(mp.exp(-z * z) * mp.erfc(-1j * z)) / (sigma * mp.sqrt(2 * mp.pi))


def _log_ncdf_asym(x):
    """log(ncdf(x)) for x < -1e4, with the asymptotic series of the Mills ratio."""
    s, t, k = mp.mpf(0), mp.mpf(1), 0
    while abs(t) > mp.eps * abs(s):
        s += t
        k += 1
        t *= -(2 * k - 1) / (x * x)
    return -x * x / 2 - mp.log(-x) - mp.log(2 * mp.pi) / 2 + mp.log(s)


def r_ndtri_exp(y):
    with wdps(20):
        if y < -1e8:
            # exp(y) is too small for mpmath, so solve log(ncdf(x)) = y directly
            x0 = -mp.sqrt(-2 * y)
            return _invert(lambda x: _log_ncdf_asym(x) - y, x0, abs(x0) * mp.mpf(1e-3))
        if y < -mp.log(2):
            return _ncdf_inv(mp.exp(y))
        return -_ncdf_inv(-mp.expm1(y))


def _owens_series(h, a):
    """T(h, a) for 0 <= a <= 1 via Owen's series, at a precision that absorbs the cancellation."""
    H = h * h / 2
    if H > 2000:
        return mp.mpf(0)
    with mp.workdps(mp.mp.dps + int(H / 2.3) + 20):
        eH = mp.exp(-H)
        part, fac, tot, j, ap = mp.mpf(1), mp.mpf(1), mp.mpf(0), 0, a
        tol = mp.mpf(10) ** (-mp.mp.dps)
        while True:
            c = (1 - eH * part) / (2 * j + 1) * ap
            tot += c if j % 2 == 0 else -c
            if abs(c) < tol * abs(mp.atan(a) - tot) * mp.mpf(1e-5) or j > 200000:
                break
            j += 1
            fac *= H / j
            part += fac
            ap *= a * a
        return (mp.atan(a) - tot) / (2 * mp.pi)


def r_owens_t(h, a):
    with wdps(20):
        h, sg, a = abs(h), mp.sign(a), abs(a)
        if a <= 1:
            return sg * _owens_series(h, a)
        # T(h, a) = (P(-h) + P(-ah)) / 2 - P(-h) P(-ah) - T(ah, 1/a), for a > 1
        rest = _owens_series(a * h, 1 / a)
        with wdps(40):
            q1, q2 = mp.ncdf(-h), mp.ncdf(-a * h)
            return sg * ((q1 + q2) / 2 - q1 * q2 - rest)


def _bvn(h, k, r):
    """P(X <= h, Y <= k) for a standard bivariate normal with correlation r, by quadrature of
    Sheppard's formula Phi(h) Phi(k) + 1/(2 pi) int_0^asin(r) exp(...) dt."""

    def f(t):
        c = mp.cos(t)
        return mp.exp(-(h * h + k * k - 2 * h * k * mp.sin(t)) / (2 * c * c))

    return mp.ncdf(h) * mp.ncdf(k) + mp.quad(f, [0, mp.asin(r)]) / (2 * mp.pi)


def _bvn_neg(h, k, r):
    """P(X <= h, Y <= k) for -1 < r < 0 and h + k < 0, without cancellation: P = 0 at r = -1,
    so by Plackett's identity P is the integral of the density phi2(h, k; s) over s in [-1, r].
    With a = (h + k)^2 / 4, b = (h - k)^2 / 4, and 1 / (1 + s) = w0 + t / a, w0 = 1 / (1 + r),
    this is int_0^inf exp(-t - a w0 - b w / (2 w - 1)) / (w sqrt(2 w - 1)) dt / (2 pi a).
    The integrand varies on the scale t ~ a w0, so this needs a w0 not too small."""
    a, b, w0 = (h + k) ** 2 / 4, (h - k) ** 2 / 4, 1 / (1 + r)
    e0 = b * w0 / (2 * w0 - 1)

    def f(t):
        w = w0 + t / a
        return mp.exp(e0 - t - b * w / (2 * w - 1)) / (w * mp.sqrt(2 * w - 1))

    s = mp.quad(f, [0, 0.25, 1, 4, 16, 64, 256, mp.inf])
    return s * mp.exp(-a * w0 - e0) / (2 * mp.pi * a)


def r_bivariate_normal_cdf(h, k, r):
    # in the tail, the quadrature of Sheppard's formula below misses the peak at t = asin(r),
    # e.g. by 1e-7 for (-9.73, -4.44, -0.913)
    if -1 < r < 0 and h + k < 0 and (h + k) ** 2 > 4 * (1 + r):
        with wdps(10):
            return _bvn_neg(h, k, r)
    # the terms of Sheppard's formula cancel where the probability is tiny, so add digits
    # until the result is well above the cancellation error
    for extra in (30, 100, 300, 1000):
        with wdps(extra):
            p = _bvn(h, k, r)
            if p != 0 and abs(p) > mp.ncdf(h) * mp.ncdf(k) * mp.mpf(10) ** (20 - extra):
                return p
    return p


def r_trivariate_normal_cdf(h1, h2, h3, r12, r13, r23):
    # Plackett's identity: scale r12 and r13 by t in [0, 1]; then Phi3 is Phi(h1) Phi2(h2, h3)
    # plus the integral over t of r12 dPhi3/dr12 + r13 dPhi3/dr13, where dPhi3/dr_ij is the
    # bivariate normal density at (h_i, h_j) times a conditional univariate normal CDF.
    # The terms cancel where the probability is tiny, so add digits until the result converges.
    prev = None
    for extra in (20, 60, 150, 400):
        with wdps(extra):
            p = _tvn(h1, h2, h3, r12, r13, r23)
        if prev is not None and abs(p - prev) <= abs(p) * mp.mpf(10) ** (
            -mp.mp.dps - 5
        ):
            return p
        prev = p
    return p


def _tvn(h1, h2, h3, r12, r13, r23):
    """Plackett's identity for r_trivariate_normal_cdf, at the working precision."""

    def phi2(a, b, r):
        return mp.exp(-(a * a - 2 * r * a * b + b * b) / (2 * (1 - r * r))) / (
            2 * mp.pi * mp.sqrt(1 - r * r)
        )

    def cond(a, b, c, rab, rac, rbc):
        # P(Z <= c | X = a, Y = b), with the correlations of (X, Y), (X, Z), (Y, Z)
        b1 = (rac - rab * rbc) / (1 - rab * rab)
        b2 = (rbc - rab * rac) / (1 - rab * rab)
        s2 = 1 - rac * b1 - rbc * b2
        return mp.ncdf((c - b1 * a - b2 * b) / mp.sqrt(s2))

    def f(t):
        p, q = t * r12, t * r13
        d12 = r12 * phi2(h1, h2, p) * cond(h1, h2, h3, p, q, r23)
        d13 = r13 * phi2(h1, h3, q) * cond(h1, h3, h2, q, p, r23)
        return d12 + d13

    return mp.ncdf(h1) * _bvn(h2, h3, r23) + mp.quad(f, [0, 1])


def r_nrdtrimn(p, std, x):
    with mp.workdps(400):
        return x - std * _ncdf_inv(p)


def r_nrdtrisd(mean, p, x):
    with mp.workdps(400):
        return (x - mean) / _ncdf_inv(p)


# ---------------------------------------------------------------- elementary


def r_cosm1(x):
    with wdps(xd(x)):
        return -2 * mp.sin(x / 2) ** 2


def r_powm1(x, y):
    return mp.powm1(x, y)


def _deg(x):
    """x mod 360 (exact, with the sign of x), divided by 180."""
    with wdps(xd(x)):
        return mp.sign(x) * mp.fmod(abs(x), 360) / 180


def r_sindg(x):
    return mp.sinpi(_deg(x))


def r_cosdg(x):
    return mp.cospi(_deg(x))


def r_tandg(x):
    t = _deg(x)
    return mp.sinpi(t) / mp.cospi(t)


def r_cotdg(x):
    t = _deg(x)
    return mp.cospi(t) / mp.sinpi(t)


def r_radian(d, m, s):
    return (d + m / 60 + s / 3600) * mp.pi / 180


def r_sinc(x):
    with wdps(xd(x)):
        return mp.sinpi(x) / (mp.pi * x)


def r_diric(n, x):
    with wdps(xd(n * x)):
        return mp.sin(n * x / 2) / (n * mp.sin(x / 2))


def r_agm(a, b):
    return mp.agm(a, b)


# ---------------------------------------------------------------- exponential and logarithmic


def r_expit(x):
    return 1 / (1 + mp.exp(-x))


def r_logit(p):
    return mp.log(p) - mp.log1p(-p)


def r_log_expit(x):
    return -mp.log1p(mp.exp(-x))


def r_exprel(x):
    return mp.expm1(x) / x


def r_log1mexp(x):
    # log1p(-exp(x)) for x < -log(2), where 1 - exp(x) can round to 1
    return mp.log1p(-mp.exp(x)) if x < -mp.log(2) else mp.log(-mp.expm1(x))


def r_log1pmx(x):
    if abs(x) < 1e-3:
        # the series -x^2/2 + x^3/3 - ..., since log1p(x) - x cancels
        s, k, xk = mp.mpf(0), 2, x * x
        while True:
            t = (-1) ** (k + 1) * xk / k
            s += t
            if abs(t) < mp.eps * abs(s):
                return s
            k, xk = k + 1, xk * x
    with wdps(10):
        return mp.log1p(x) - x


def r_xlogy(x, y):
    return x * mp.log(y)


def r_xlog1py(x, y):
    return x * mp.log1p(y)


def c_xlogy(xr, xi, yr, yi):
    """x ln(y), but 0 if x = 0, like SciPy (also for y = 0)."""
    x = mp.mpc(xr, xi)
    return x * mp.log(mp.mpc(yr, yi)) if x else mp.mpc(0)


def c_xlog1py(xr, xi, yr, yi):
    """x ln(1 + y), but 0 if x = 0, like SciPy (also for y = -1)."""
    x = mp.mpc(xr, xi)
    return x * mp.log1p(mp.mpc(yr, yi)) if x else mp.mpc(0)


def r_softplus(x):
    return mp.log1p(mp.exp(x))


def r_logaddexp(x, y):
    # max + log1p(exp(min - max)), since 1 + exp(min - max) can round to 1
    with wdps(20):
        return max(x, y) + mp.log1p(mp.exp(-abs(x - y)))


def r_logaddexp2(x, y):
    with wdps(20):
        return max(x, y) + mp.log1p(mp.power(2, -abs(x - y))) / mp.log(2)


def r_boxcox(x, lmbda):
    with wdps(20):
        return mp.expm1(lmbda * mp.log(x)) / lmbda


def r_boxcox1p(x, lmbda):
    with wdps(20):
        return mp.expm1(lmbda * mp.log1p(x)) / lmbda


def r_inv_boxcox(y, lmbda):
    with wdps(20):
        return mp.exp(mp.log1p(lmbda * y) / lmbda)


def r_inv_boxcox1p(y, lmbda):
    with wdps(20):
        return mp.expm1(mp.log1p(lmbda * y) / lmbda)


def r_entr(x):
    return -x * mp.log(x)


def r_rel_entr(x, y):
    with wdps(20):
        return x * (mp.log(x) - mp.log(y))


def r_kl_div(x, y):
    with wdps(40):
        return x * (mp.log(x) - mp.log(y)) - x + y


def r_huber(delta, r):
    return r * r / 2 if abs(r) <= delta else delta * (abs(r) - delta / 2)


def r_pseudo_huber(delta, r):
    return r * r / (mp.sqrt(1 + (r / delta) ** 2) + 1)


# ---------------------------------------------------------------- Bessel


def _sc(f, v, x):
    """I_v(x) e^-|x| or K_v(x) e^x, for f = mp.besseli or mp.besselk."""
    with wdps(xd(v, x)):
        return f(v, x) * mp.exp(-abs(x) if f is mp.besseli else x)


def r_i0e(x):
    return _sc(mp.besseli, 0, x)


def r_i1e(x):
    return _sc(mp.besseli, 1, x)


def r_k0e(x):
    return _sc(mp.besselk, 0, x)


def r_k1e(x):
    return _sc(mp.besselk, 1, x)


def r_ie(v, x):
    return _sc(mp.besseli, v, x)


def r_ke(v, x):
    return _sc(mp.besselk, v, x)


def _bessel_diff(f, v, x, n, c):
    """The n-th derivative c^n 2^-n sum_k c^k binom(n, k) f(v - n + 2k, x), DLMF 10.6.7 and
    10.29.5, with c = -1 for J, Y, and K, and c = 1 for I."""
    n = int(n)
    with wdps(10):
        s = mp.fsum(
            c**k * mp.binomial(n, k) * f(v - n + 2 * k, x) for k in range(n + 1)
        )
        return (s if f is not mp.besselk else (-1) ** n * s) / 2**n


def r_jvp(v, x, n):
    with wdps(xd(v, x)):
        return _osc(_bessel_diff(mp.besselj, v, x, n, -1), v, x)


def r_yvp(v, x, n):
    with wdps(xd(v, x)):
        return _osc(_bessel_diff(mp.bessely, v, x, n, -1), v, x)


def r_ivp(v, x, n):
    with wdps(xd(v, x)):
        return _bessel_diff(mp.besseli, v, x, n, 1)


def r_kvp(v, x, n):
    with wdps(xd(v, x)):
        return _bessel_diff(mp.besselk, v, x, n, 1)


def _iv_ratio(v, x):
    return mp.besseli(v, x) / mp.besseli(v - 1, x)


def r_iv_ratio(v, x):
    with wdps(xd(v, x)):
        return _iv_ratio(v, x)


def r_iv_ratio_c(v, x):
    with wdps(2 * xd(v, x)):
        return 1 - _iv_ratio(v, x)


def r_iv_ratioinv(v, r):
    if r == 0:
        return mp.mpf(0)
    return _invert_pos(lambda x: _iv_ratio(v, x), r, v)


def _sph(f, n, x, derivative=False):
    """sqrt(pi / (2x)) f(n + 1/2, x), or its derivative, as n/x f_n(x) - f_(n+1)(x) (or with
    +f_(n+1)(x) for i_n), DLMF 10.51.2 and 10.51.5, without cancellation for small x."""
    v = n + mp.mpf(1) / 2
    c = mp.sqrt(mp.pi / (2 * x))
    if derivative:
        return c * (n / x * f(v, x) + (1 if f is mp.besseli else -1) * f(v + 1, x))
    return c * f(v, x)


def _sph_osc(f, n, x, derivative=False):
    with wdps(xd(n, x)):
        y = _sph(f, n, x, derivative)
        if x > n + 0.5:
            return y, mp.sqrt(mp.pi / (2 * x)) * _envelope(n + mp.mpf(1) / 2, x)
        return y


def r_sph_j(n, x):
    return _sph_osc(mp.besselj, n, x)


def r_sph_y(n, x):
    return _sph_osc(mp.bessely, n, x)


def r_sph_i(n, x):
    with wdps(xd(n, x)):
        return _sph(mp.besseli, n, x)


def r_sph_k(n, x):
    with wdps(xd(n, x)):
        return _sph(mp.besselk, n, x)


def r_sph_jp(n, x):
    return _sph_osc(mp.besselj, n, x, True)


def r_sph_yp(n, x):
    return _sph_osc(mp.bessely, n, x, True)


def r_sph_ip(n, x):
    with wdps(xd(n, x)):
        return _sph(mp.besseli, n, x, True)


def r_sph_kp(n, x):
    with wdps(xd(n, x)):
        return _sph(mp.besselk, n, x, True)


def _riccati(f, n, x, derivative=False):
    """x f_n(x) for the spherical Bessel function f_n of `_sph`, or its derivative
    x f_(n-1)(x) - n f_n(x) (DLMF 10.51.2), with the scale sqrt(pi x / 2) times the Bessel
    modulus for x > n + 1/2. Unlike (n + 1) f_n(x) - x f_(n+1)(x), this doesn't cancel for y_0 at
    small x."""
    with wdps(xd(n, x)):
        if derivative:
            y = x * _sph(f, n - 1, x) - n * _sph(f, n, x)
        else:
            y = x * _sph(f, n, x)
        v = n + mp.mpf(1) / 2
        return (y, mp.sqrt(mp.pi * x / 2) * _envelope(v, x)) if x > v else y


def r_riccati_j(n, x):
    return _riccati(mp.besselj, n, x)


def r_riccati_y(n, x):
    return _riccati(mp.bessely, n, x)


def r_riccati_jp(n, x):
    return _riccati(mp.besselj, n, x, True)


def r_riccati_yp(n, x):
    return _riccati(mp.bessely, n, x, True)


def _lmbda(v, x, derivative=False):
    """The Jahnke-Emden function c J_v(x) with c = Gamma(v + 1) (2/x)^v, or its derivative
    -c J_(v+1)(x) (DLMF 10.6.6), with the scale c times the Bessel modulus where it oscillates."""
    with wdps(xd(v, x)):
        c = mp.gamma(v + 1) * (2 / x) ** v
        u = v + 1 if derivative else v
        y = (-c if derivative else c) * mp.besselj(u, x)
        return (y, c * _envelope(u, x)) if x > u else y


def r_lmbda(v, x):
    return _lmbda(v, x)


def r_lmbdap(v, x):
    return _lmbda(v, x, True)


# ---------------------------------------------------------------- integrals of Bessel and Struve


def _int_j0(f, x):
    # int_0^x C_0(t) dt = x C_0(x) + pi x / 2 (C_1(x) H_0(x) - C_0(x) H_1(x)), DLMF 10.22.2
    with wdps(xd(x) + 10):
        c0, c1 = f(0, x), f(1, x)
        return x * c0 + mp.pi * x / 2 * (c1 * mp.struveh(0, x) - c0 * mp.struveh(1, x))


def r_itj0(x):
    return _int_j0(mp.besselj, x)


def r_ity0(x):
    return _int_j0(mp.bessely, x)


def r_it2j0(x):
    # int_0^x (1 - J_0(t)) / t dt = x^2/8 2F3(1, 1; 2, 2, 2; -x^2/4)
    with wdps(xd(x) + 10):
        return x * x / 8 * mp.hyper([1, 1], [2, 2, 2], -x * x / 4)


def r_it2y0(x):
    # int_x^inf Y_0(t) / t dt = pi/6 - G(x), where G is the antiderivative of Y_0(t) / t that
    # follows from the power series of Y_0 (DLMF 10.8.2), with L = log(x/2) + euler:
    #   G(x) = L^2/pi + 2/pi sum_k>=1 (-1)^k (x/2)^2k / k!^2 (L/(2k) - 1/(4k^2) - H_k/(2k)).
    # The terms grow to about e^x, so this needs about x log10(e) extra digits.
    with wdps(int(0.4343 * x) + 20):
        L = mp.log(x / 2) + mp.euler
        s, h, k, p = mp.mpf(0), mp.mpf(0), 1, mp.mpf(1)
        tol = mp.mpf(10) ** (-mp.mp.dps)
        while True:
            p *= (x / 2) ** 2 / (k * k)
            h += mp.mpf(1) / k
            t = (-1) ** k * p * (L / (2 * k) - mp.mpf(1) / (4 * k * k) - h / (2 * k))
            s += t
            if k > x and abs(t) < tol:
                return mp.pi / 6 - L * L / mp.pi - 2 / mp.pi * s
            k += 1


def r_iti0(x):
    # int_0^x I_0(t) dt = x 1F2(1/2; 1, 3/2; x^2/4)
    with wdps(xd(x) + 10):
        return x * mp.hyp1f2(0.5, 1, 1.5, x * x / 4)


def r_itk0(x):
    # int_0^x K_0(t) dt = x K_0(x) + pi x / 2 (K_1(x) L_0(x) + K_0(x) L_1(x)), DLMF 10.43.2
    with wdps(xd(x) + 10):
        k0, k1 = mp.besselk(0, x), mp.besselk(1, x)
        return x * k0 + mp.pi * x / 2 * (k1 * mp.struvel(0, x) + k0 * mp.struvel(1, x))


def r_it2i0(x):
    # int_0^x (I_0(t) - 1) / t dt = x^2/8 2F3(1, 1; 2, 2, 2; x^2/4)
    with wdps(xd(x) + 10):
        return x * x / 8 * mp.hyper([1, 1], [2, 2, 2], x * x / 4)


def r_it2k0(x):
    # int_x^inf K_0(t) / t dt = pi^2/24 - G(x), where G is the antiderivative of K_0(t) / t that
    # follows from the power series of K_0 (DLMF 10.31.2), with L = log(x/2) + euler:
    #   G(x) = -L^2/2 - sum_k>=1 (x/2)^2k / k!^2 (L/(2k) - 1/(4k^2) - H_k/(2k)).
    # The terms grow to about e^x, and the result is about e^-x, so this needs about
    # 2x log10(e) extra digits.
    with wdps(int(0.8686 * x) + 20):
        L = mp.log(x / 2) + mp.euler
        s, h, k, p = mp.mpf(0), mp.mpf(0), 1, mp.mpf(1)
        tol = mp.mpf(10) ** (-mp.mp.dps)
        while True:
            p *= (x / 2) ** 2 / (k * k)
            h += mp.mpf(1) / k
            t = p * (L / (2 * k) - mp.mpf(1) / (4 * k * k) - h / (2 * k))
            s += t
            if k > x and abs(t) < tol:
                return mp.pi**2 / 24 + L * L / 2 + s
            k += 1


def r_besselpoly(a, lmbda, v):
    # int_0^1 t^lambda J_v(2 a t) dt
    #   = a^v / ((lambda + v + 1) Gamma(v + 1)) 1F2(c; v + 1, c + 1; -a^2), c = (lambda + v + 1)/2
    with wdps(xd(a) * 2 + 10):
        c = (lmbda + v + 1) / 2
        return (
            a**v
            / ((lmbda + v + 1) * mp.gamma(v + 1))
            * mp.hyp1f2(c, v + 1, c + 1, -a * a)
        )


def r_itstruve0(x):
    # int_0^x H_0(t) dt = x^2/pi 2F3(1, 1; 3/2, 3/2, 2; -x^2/4)
    with wdps(xd(x) + 10):
        return x * x / mp.pi * mp.hyper([1, 1], [1.5, 1.5, 2], -x * x / 4)


def r_it2struve0(x):
    # int_x^inf H_0(t) / t dt = pi/2 - 2x/pi 2F3(1/2, 1; 3/2, 3/2, 3/2; -x^2/4)
    with wdps(xd(x) + 10):
        return mp.pi / 2 - 2 * x / mp.pi * mp.hyper(
            [0.5, 1], [1.5, 1.5, 1.5], -x * x / 4
        )


def r_itmodstruve0(x):
    # int_0^x L_0(t) dt = x^2/pi 2F3(1, 1; 3/2, 3/2, 2; x^2/4)
    with wdps(xd(x) + 10):
        return x * x / mp.pi * mp.hyper([1, 1], [1.5, 1.5, 2], x * x / 4)


# ---------------------------------------------------------------- Kelvin and Airy functions


def _kelvin(kind, x):
    """ber, bei, ker, kei (kind 0 to 3), or their derivatives (kind 4 to 7), with the modulus of
    the complex function (e.g. |ber + i bei|) as scale. ber and bei are even.

    ker_v + i kei_v = (-i)^v K_v(x e^(i pi/4)) (DLMF 10.61.2), since mpmath's ker and kei are slow
    for large x: they work at up to 2600 digits for x = 1000."""
    s2 = mp.sqrt(2)
    # the sums of the order-1 functions cancel for small x, by up to a factor x^2
    with wdps(xd(x) + 2 * max(0, int(-mp.log10(abs(x))))):
        a = abs(x)
        if kind < 2:
            f = mp.ber(0, a), mp.bei(0, a)
        elif kind < 4:
            k = mp.besselk(0, a * mp.expjpi(0.25))
            f = k.real, k.imag
        elif kind < 6:  # ber'(x) = (ber_1(x) + bei_1(x)) / sqrt(2), DLMF 10.63.1
            b1, c1 = mp.ber(1, a), mp.bei(1, a)
            f = mp.sign(x) * (b1 + c1) / s2, mp.sign(x) * (c1 - b1) / s2
        else:
            k = mp.besselk(1, a * mp.expjpi(0.25))
            b1, c1 = k.imag, -k.real
            f = (b1 + c1) / s2, (c1 - b1) / s2
        return f[kind % 2], mp.hypot(*f)


def _airy(x, scaled=False):
    """Ai, Ai', Bi, Bi', with the scale sqrt(Ai^2 + Bi^2) (or of the derivatives) for x < 0."""
    with wdps(xd(x) + 10):
        ai, aip = mp.airyai(x), mp.airyai(x, derivative=1)
        bi, bip = mp.airybi(x), mp.airybi(x, derivative=1)
        if scaled and x > 0:
            z = 2 * x * mp.sqrt(x) / 3
            ai, aip, bi, bip = (
                ai * mp.exp(z),
                aip * mp.exp(z),
                bi * mp.exp(-z),
                bip * mp.exp(-z),
            )
        if x >= 0:
            return [ai, aip, bi, bip]
        s, sp = mp.hypot(ai, bi), mp.hypot(aip, bip)
        return [(ai, s), (aip, sp), (bi, s), (bip, sp)]


def _airy_int(x):
    """int_0^x Ai(t), Bi(t), Ai(-t), and Bi(-t) dt. With Ai = c1 f - c2 g and
    Bi = sqrt(3) (c1 f + c2 g), DLMF 9.4.1-9.4.3, these follow from
    int_0^t f = t 1F2(1/3; 2/3, 4/3; t^3/9) and int_0^t g = t^2/2 1F2(2/3; 4/3, 5/3; t^3/9);
    mpmath's airyai(x, derivative=-1) does not converge for |x| >~ 300. For int_0^x Ai(t) dt,
    where the two terms cancel for large x, use mpmath for x <= 100, and otherwise
    1/3 - int_x^inf Ai(t) dt with int_x^inf Ai(t) dt ~ Ai(x) / sqrt(x) (< 1e-280)."""
    third = mp.mpf(1) / 3
    with wdps(xd(x) + 20):
        c1 = 1 / (mp.cbrt(9) * mp.gamma(2 * third))
        c2 = 1 / (mp.cbrt(3) * mp.gamma(third))

        def fg(t):
            z = t**3 / 9
            return c1 * t * mp.hyp1f2(
                third, 2 * third, 4 * third, z
            ), c2 * t * t / 2 * mp.hyp1f2(2 * third, 4 * third, 5 * third, z)

        f, g = fg(x)
        fn, gn = fg(-x)
        ai = (
            mp.airyai(x, derivative=-1)
            if x <= 100
            else third - mp.airyai(x) / mp.sqrt(x)
        )
        return [ai, mp.sqrt(3) * (f + g), gn - fn, -mp.sqrt(3) * (fn + gn)]


# ---------------------------------------------------------------- parabolic cylinder functions


def r_pbdv(v, x):
    # D'_v(x) = x/2 D_v(x) - D_{v+1}(x), DLMF 12.8.3
    with wdps(20):
        d = mp.pcfd(v, x)
        return d, x / 2 * d - mp.pcfd(v + 1, x)


def r_pbvv(v, x):
    # V_v(x) is V(-v - 1/2, x) in DLMF notation
    with wdps(20):
        a = -v - mp.mpf(1) / 2
        return mp.pcfv(a, x), mp.diff(lambda t: mp.pcfv(a, t), x)


def c_pbdn(n, re, im):
    # D_n(z) and D'_n(z) = z/2 D_n(z) - D_{n+1}(z) (DLMF 12.8.3)
    with wdps(20 + xd(re, im)):
        z = mp.mpc(re, im)
        d = mp.pcfd(n, z)
        return d, z / 2 * d - mp.pcfd(n + 1, z)


def r_pbwa(a, x):
    with wdps(20):
        return mp.pcfw(a, x), mp.diff(lambda t: mp.pcfw(a, t), x)


# ---------------------------------------------------------------- Mathieu functions


def _mathieu_matrix(odd, m, q, n):
    """The first `n` rows of the symmetric tridiagonal matrix of the recurrence relations of
    the Fourier coefficients (DLMF 28.4.5-8) of ce_m (or se_m, if `odd`), as its diagonal, its
    off-diagonal, and the harmonics. For ce_2j it is symmetrized with sqrt(2) A_0, so that the
    normalization (DLMF 28.4.9-12) is the 2-norm in each case."""
    h0 = m % 2 or 2 * odd
    harm = [h0 + 2 * k for k in range(n)]
    d = [h * h for h in harm]
    e = [q] * n
    if h0 == 0:
        e[0] = math.sqrt(2) * q if isinstance(q, float) else mp.sqrt(2) * q
    elif h0 == 1:
        d[0] = 1 - q if odd else 1 + q
    return d, e, harm


def _mathieu_size(h0, lam, q, tol, k0=0):
    """A number of rows after which the Fourier coefficients are smaller than `tol` times the
    one of index `k0` (or of the turning point, if later): beyond it, each is at most
    |q| / (h^2 - lam - |q|) times the previous one, for the harmonic h."""
    q, k, tail = abs(q), 0, 1.0
    while tail > tol:
        gap = (h0 + 2 * k) ** 2 - lam - q
        if k > k0 and gap > 2 * q:
            tail *= q / gap
        k += 1
    return k + 2


def _mathieu_bisect(d, e, i):
    """Eigenvalue `i` (from 0, ascending) of the symmetric tridiagonal matrix, in f64, by
    bisection on the number of negative pivots (a Sturm sequence)."""
    r = 2 * max(map(abs, e))
    lo, hi = min(d) - r, sorted(d)[i] + r
    while hi - lo > 1e-13 * max(abs(lo), abs(hi), 1.0):
        x = (lo + hi) / 2
        count, p = 0, 1.0
        for k, dk in enumerate(d):
            p = dk - x - (e[k - 1] ** 2 / p if k else 0.0)
            p = p or 1e-300
            count += p < 0
        if count > i:
            hi = x
        else:
            lo = x
    return (lo + hi) / 2


def _mathieu_ratios(lam, d, e, p):
    """The solution x of the recurrence relations relative to x_p, by forward recurrence of
    x_(k-1) / x_k up to k = p and backward recurrence of x_(k+1) / x_k down to it (both stable
    if x_p is the largest), and the residual of row p, which is zero iff `lam` is an
    eigenvalue."""
    n = len(d)
    x = [mp.mpf(1)] * n
    s = r = mp.mpf(0)
    for k in range(p):
        s = e[k] / (lam - d[k] - (e[k - 1] * s if k else 0))
        x[k] = s
    for k in range(n - 1, p, -1):
        r = e[k - 1] / (lam - d[k] - e[k] * r)
        x[k] = r
    for k in range(p - 1, -1, -1):
        x[k] *= x[k + 1]
    for k in range(p + 1, n):
        x[k] *= x[k - 1]
    return x, lam - d[p] - (e[p - 1] * s if p else 0) - e[p] * r


def _mathieu_eig(odd, m, q, k0=0):
    """The characteristic value a_m(q) (or b_m(q), if `odd`), the Fourier coefficients of ce_m
    (se_m) as (harmonic, coefficient) pairs, and the index of the largest one. Each coefficient
    is accurate to the working precision relative to itself, up to index `k0` at least.

    The f64 bisection only selects the eigenvalue; it is found with Newton's method on the
    residual g of `_mathieu_ratios`, whose derivative is the sum of x^2."""
    i = (m - odd) // 2
    qf = float(q)
    h0 = m % 2 or 2 * odd
    # the size of the f64 matrix, from an upper bound of the eigenvalue
    n = _mathieu_size(h0, (h0 + 2 * i) ** 2 + 4 * abs(qf), qf, 1e-20)
    lam0 = _mathieu_bisect(*_mathieu_matrix(odd, m, qf, n)[:2], i)
    tol = mp.mpf(10) ** -(mp.mp.dps - 5)
    n = _mathieu_size(h0, lam0, qf, float(tol) ** 2, k0)
    d, e, harm = _mathieu_matrix(odd, m, q, n)
    lam, p = mp.mpf(lam0), min(range(n), key=lambda k: abs(lam0 - d[k]))
    for _ in range(50):
        x, g = _mathieu_ratios(lam, d, e, p)
        pmax = max(range(n), key=lambda k: abs(x[k]))
        if pmax != p:
            p = pmax
            continue
        norm2 = mp.fsum(t * t for t in x)
        lam -= g / norm2
        if abs(g) <= tol * norm2 * abs(lam):
            break
    else:
        raise NoConvergence
    if abs(lam - lam0) > 1e-6 * max(abs(lam0), 1.0):
        raise NoConvergence
    # the sign of the first coefficient is that of q^i, by continuity in q from q = 0
    # (DLMF 28.4.13-16), as it is never zero; for q > 0 it is positive, as in specfun
    s = mp.sign(x[0]) * (-1 if q < 0 and i % 2 else 1) / mp.sqrt(norm2)
    c = [s * t for t in x]
    if h0 == 0:
        c[0] /= mp.sqrt(2)
    return lam, list(zip(harm, c, strict=True)), p


def _mathieu_cv(odd, m, q):
    with wdps(15):
        return _mathieu_eig(odd, int(m), q)[0]


def r_mathieu_a(m, q):
    return _mathieu_cv(0, m, q)


def r_mathieu_b(m, q):
    return _mathieu_cv(1, m, q)


def _mathieu_coef(odd, m, q, k):
    with wdps(15):
        return _mathieu_eig(odd, int(m), q, int(k))[1][int(k)][1]


def r_mathieu_even_coef(m, q, k):
    return _mathieu_coef(0, m, q, k)


def r_mathieu_odd_coef(m, q, k):
    return _mathieu_coef(1, m, q, k)


def _mathieu_sum(f):
    """The value of f() -> (value, terms, scale), where the value is the sum of the terms (each
    accurate to the working precision), with as many extra digits as the sum cancels relative
    to the scale, and so that the last terms are negligible."""
    extra = 15
    while True:
        with wdps(extra):
            y, terms, scale = f()
            lost = mp.log10(mp.fsum(map(abs, terms)) / scale)
            last = max(map(abs, terms[-3:])) / scale
            if lost + 15 <= extra and last < mp.mpf(10) ** -(mp.mp.dps - 10):
                return y
        extra = max(int(lost) + 25, extra + 20)
        if mp.mp.dps + extra > 300:
            raise NoConvergence


def _mathieu_angular(odd, m, q, z, derivative=False):
    """ce_m(z, q) (or se_m(z, q), if `odd`) or its derivative, as the Fourier series
    DLMF 28.4.1-4, with the scale sqrt(y^2 + y'^2 / w), or sqrt(w y^2 + y'^2) for the
    derivative, where w = max(|a - 2 q cos(2 z)|, 1) is the squared local frequency of
    y'' + (a - 2 q cos(2 z)) y = 0 (DLMF 28.2.1)."""

    def f():
        lam, coef, _ = _mathieu_eig(odd, int(m), q)
        if odd:
            y = [c * mp.sin(h * z) for h, c in coef]
            dy = [h * c * mp.cos(h * z) for h, c in coef]
        else:
            y = [c * mp.cos(h * z) for h, c in coef]
            dy = [-h * c * mp.sin(h * z) for h, c in coef]
        w = max(abs(lam - 2 * q * mp.cos(2 * z)), 1)
        y1, dy1 = mp.fsum(y), mp.fsum(dy)
        scale = mp.sqrt(y1 * y1 + dy1 * dy1 / w)
        if derivative:
            return (dy1, mp.sqrt(w) * scale), dy, mp.sqrt(w) * scale
        return (y1, scale), y, scale

    return _mathieu_sum(f)


def r_mathieu_cem(m, q, x):
    return _mathieu_angular(0, m, q, mp.radians(x))


def r_mathieu_cemp(m, q, x):
    return _mathieu_angular(0, m, q, mp.radians(x), True)


def r_mathieu_sem(m, q, x):
    return _mathieu_angular(1, m, q, mp.radians(x))


def r_mathieu_semp(m, q, x):
    return _mathieu_angular(1, m, q, mp.radians(x), True)


def r_mathieu_ce(m, q, x):
    return _mathieu_angular(0, m, q, x)


def r_mathieu_cep(m, q, x):
    return _mathieu_angular(0, m, q, x, True)


def r_mathieu_se(m, q, x):
    return _mathieu_angular(1, m, q, x)


def r_mathieu_sep(m, q, x):
    return _mathieu_angular(1, m, q, x, True)


def _mathieu_bessel(f, n, x):
    """f_k(x) for k = 0, ..., n and their derivatives, for f = J or Y, by recurrence in the
    stable direction: backward for J, forward for Y."""
    if f is mp.besselj:
        b = [f(n + 1, x), f(n, x)]
        for k in range(n, 0, -1):
            b.append(2 * k / x * b[-1] - b[-2])
        b.reverse()
    else:
        b = [f(0, x), f(1, x)]
        for k in range(1, n + 1):
            b.append(2 * k / x * b[-1] - b[-2])
    return b, [-b[1]] + [(b[k - 1] - b[k + 1]) / 2 for k in range(1, n + 1)]


def _mathieu_radial(odd, kind, m, q, z, derivative=False):
    """Mc_m^(kind)(z, h) (or Ms_m^(kind)(z, h), if `odd`) with h = sqrt(q), or its derivative,
    as the series of Bessel function products DLMF 28.24.2-5, with s the index of the largest
    Fourier coefficient, so that the series doesn't cancel (unlike for s = 0, as in specfun).
    The scale, where the function oscillates (2 q cosh(2 z) > a), is the modulus of the first
    and second kind."""

    def f():
        lam, coef, s = _mathieu_eig(odd, int(m), q)
        n, h0 = len(coef), coef[0][0]
        u1, u2 = mp.sqrt(q) * mp.exp(-z), mp.sqrt(q) * mp.exp(z)
        j1, dj1 = _mathieu_bessel(mp.besselj, n + s + h0, u1)

        def at(b, k):
            # J_(-k) = (-1)^k J_k, and so for Y and the derivatives
            return b[k] if k >= 0 else (-1) ** k * b[-k]

        def series(c2, dc2):
            terms = []
            for k in range(n):
                t = 0
                for sign, v1, v2 in [
                    (1, k - s, k + s + h0),
                    (-1 if odd else 1, k + s + h0, k - s),
                ]:
                    if derivative:
                        t += sign * u2 * at(j1, v1) * at(dc2, v2)
                        t -= sign * u1 * at(dj1, v1) * at(c2, v2)
                    else:
                        t += sign * at(j1, v1) * at(c2, v2)
                t *= (-1) ** (k + (int(m) - odd) // 2) * coef[k][1] / coef[s][1]
                terms.append(t / 2 if h0 == s == 0 else t)
            return terms

        t1 = series(*_mathieu_bessel(mp.besselj, n + s + h0, u2))
        t2 = series(*_mathieu_bessel(mp.bessely, n + s + h0, u2))
        y1, y2 = mp.fsum(t1), mp.fsum(t2)
        y, terms = (y1, t1) if kind == 1 else (y2, t2)
        if 2 * q * mp.cosh(2 * z) > lam:
            scale = mp.sqrt(y1 * y1 + y2 * y2)
            return (y, scale), terms, scale
        return y, terms, abs(y)

    return _mathieu_sum(f)


# ---------------------------------------------------------------- spheroidal wave functions

# mpmath and Arb don't have these. With kd = 1 for the prolate and -1 for the oblate ones, they are
# the functions of DLMF 30 with gamma^2 = kd c^2, in the normalization of Flammer, "Spheroidal
# Wave Functions" (1957), which SciPy uses.


def _sw_rows(m, n, g2, L):
    """The rows (a, b, c) of the recurrence a e_(l-2) + (b - cv) e_l + c e_(l+2) = 0 of the
    coefficients e_l of P_l^m in the angular function (DLMF 30.8.1) with gamma^2 = g2, for
    l = n - 2 floor((n - m) / 2), ..., L: -A_k, B_k + gamma^2, and -C_k of DLMF 30.8.3 with
    l = n + 2k, since e_l = (-1)^k a_(n,k) and cv = lambda + gamma^2."""
    rows = []
    for l in range(m + (n - m) % 2, L + 1, 2):
        ll = l * (l + 1)
        a = g2 * ((l - m - 1) * (l - m)) / ((2 * l - 3) * (2 * l - 1))
        b = ll + g2 * (1 - mp.mpf(2 * (ll - 1 + m * m)) / ((2 * l - 1) * (2 * l + 3)))
        c = g2 * ((l + m + 1) * (l + m + 2)) / ((2 * l + 3) * (2 * l + 5))
        rows.append((a, b, c))
    return rows


def _sw_count(rows, x):
    """The number of eigenvalues below x of the tridiagonal matrix `rows`: the negative terms of
    its Sturm sequence."""
    count, q, c0 = 0, 1, 0
    for a, b, c in rows:
        q = b - x - a * c0 / q or 1e-300
        count += q < 0
        c0 = c
    return count


def _sw_guess(rows, i):
    """The eigenvalue of index i of the tridiagonal matrix `rows`, by bisection in f64."""
    rows = [tuple(map(float, row)) for row in rows]
    lo = min(b - abs(a) - abs(c) for a, b, c in rows)
    hi = max(b + abs(a) + abs(c) for a, b, c in rows)
    while True:
        x = (lo + hi) / 2
        if not lo < x < hi:
            return x
        if _sw_count(rows, x) > i:
            hi = x
        else:
            lo = x


def _sw_coef(rows, i, cv):
    """The coefficients e_l of `_sw_rows` with e_n = 1 (row i), from the ratios e_l / e_(l+2)
    below n (upward) and e_l / e_(l-2) above n (downward), which are the stable directions; and
    the Newton step f / f' of the residual f(cv) of row i, which is zero at an eigenvalue."""
    e = [mp.mpf(1)] * len(rows)
    r = dr = 0  # the ratio, and its derivative w.r.t. cv
    for j in range(i):
        a, b, c = rows[j]
        d = b - cv + a * r
        r, dr = -c / d, c * (a * dr - 1) / d**2
        e[j] = r
    a = rows[i][0]
    f, df = rows[i][1] - cv + a * r, a * dr - 1
    r = dr = 0
    for j in range(len(rows) - 1, i, -1):
        a, b, c = rows[j]
        d = b - cv + c * r
        r, dr = -a / d, a * (c * dr - 1) / d**2
        e[j] = r
    c = rows[i][2]
    f, df = f + c * r, df + c * dr
    for j in range(i - 1, -1, -1):
        e[j] *= e[j + 1]
    for j in range(i + 1, len(rows)):
        e[j] *= e[j - 1]
    return e, f / df


def _sw_cv(m, n, g2):
    """The characteristic value of SciPy, lambda_n^m(gamma^2) + gamma^2 of DLMF 30.3 with
    gamma^2 = g2, the coefficients e_l of `_sw_coef`, and their last degree L.

    The matrix of `_sw_rows` is similar to a symmetric one with a nonzero off-diagonal, so its
    eigenvalues are simple, and ordered like at gamma = 0, where they are l (l + 1). The
    truncation at L is doubled until e_L is negligible."""
    i = (n - m) // 2
    tol = mp.mpf(10) ** -(mp.mp.dps + 10)
    # e_(l+2) / e_l ~ -g2 / (4 l^2) (DLMF 30.8.7), so e_L is of the order of (c e / (2 L))^L
    k, c, cv = 4, float(mp.sqrt(abs(g2))), None
    while (n + 2 * k) * math.log10((n + 2 * k) / (1.4 * (1 + c))) < (
        mp.mp.dps + 12 + 2 * m * math.log10(n + 2 * k)
    ):
        k += 1
    while True:
        L = n + 2 * k
        rows = _sw_rows(m, n, g2, L)
        if cv is None:
            cv = mp.mpf(_sw_guess(rows, i))
        for _ in range(30):
            e, step = _sw_coef(rows, i, cv)
            cv -= step
            if abs(step) <= 1000 * mp.eps * max(1, abs(cv)):
                break
        else:
            raise NoConvergence
        if abs(e[-1]) * math.perm(L + m, 2 * m) <= tol * max(map(abs, e)):
            break
        k *= 2
    d = mp.sqrt(tol) * max(1, abs(cv))
    if (_sw_count(rows, cv - d), _sw_count(rows, cv + d)) != (i, i + 1):
        raise NoConvergence
    return cv, e, L


def _sw_legendre(m, L, x):
    """P_l^m(x) for l = 0, ..., L without the Condon-Shortley phase, by upward recurrence
    (DLMF 14.10.3)."""
    p = [mp.mpf(0)] * (L + 1)
    p[m] = mp.fac2(2 * m - 1) * mp.sqrt(1 - x * x) ** m
    for l in range(m, L):
        p[l + 1] = ((2 * l + 1) * x * p[l] - (l + m) * p[l - 1]) / (l - m + 1)
    return p


def _sw_ang(m, n, c, x, kd):
    """The angular function of the first kind S_mn(c, x), its derivative, and the characteristic
    value. Normalized like in Flammer: S and P_n^m, or their derivatives if n - m is odd, are
    equal at x = 0. The extra digits are for the cancellation, which grows like e^c."""
    m, n = int(m), int(n)
    with wdps(10 + int(c) // 2):
        cv, e, L = _sw_cv(m, n, kd * c * c)
        ls = range(m + (n - m) % 2, L + 1, 2)
        p0, p = _sw_legendre(m, L, mp.mpf(0)), _sw_legendre(m, L, x)
        # the values at 0, or the derivatives (DLMF 14.10.5)
        at0 = [(l + m) * p0[l - 1] if (n - m) % 2 else p0[l] for l in ls]
        norm = at0[(n - m) // 2] / mp.fdot(e, at0)
        s = mp.fdot(e, [p[l] for l in ls])
        ds = mp.fdot(e, [(l + m) * p[l - 1] - l * x * p[l] for l in ls]) / (1 - x * x)
        return norm * s, norm * ds, cv


def _sw_angular(m, n, c, x, kd):
    """The angular function of the first kind S and its derivative. Where S oscillates, which
    is where q = cv - kd c^2 x^2 - m^2 / p > 0 in (p S')' + q S = 0 with p = 1 - x^2, their
    scales are A = sqrt(S^2 + p S'^2 / q), the envelope of S through its local extrema, and
    sqrt(q / p) A. A diverges at a turning point (q = 0), so q is at least (p q'^2)^(1/3),
    which it is at one Airy width from it."""
    s, ds, cv = _sw_ang(m, n, c, x, kd)
    p = 1 - x * x
    q = cv - kd * (c * x) ** 2 - m * m / p
    if q <= 0:
        return s, ds
    dq = -2 * x * (kd * c * c + (m / p) ** 2)
    q = max(q, mp.cbrt(p * dq * dq))
    a = mp.sqrt(s * s + p * ds * ds / q)
    return (s, a), (ds, mp.sqrt(q / p) * a)


def _sw_sph(kf, L, z):
    """The spherical Bessel functions j_l(z) (kf = 1) or y_l(z) (kf = 2) for l = 0, ..., L + 1,
    by downward resp. upward recurrence (DLMF 10.51.1), the stable directions."""
    f = [mp.mpf(0)] * (L + 2)
    if kf == 1:
        s = mp.sqrt(mp.pi / (2 * z))
        f[L + 1], f[L] = s * mp.besselj(L + 1.5, z), s * mp.besselj(L + 0.5, z)
        for l in range(L, 0, -1):
            f[l - 1] = (2 * l + 1) / z * f[l] - f[l + 1]
    else:
        f[0] = -mp.cos(z) / z
        f[1] = (f[0] - mp.sin(z)) / z
        for l in range(1, L + 1):
            f[l + 1] = (2 * l + 1) / z * f[l] - f[l - 1]
    return f


def _sw_series(kf, m, n, c, x, kd, e, L):
    """The radial function of the first or second kind (kf) and its derivative, by the series
    of spherical Bessel functions DLMF 30.11.3 up to the degree L, with z = i x and
    gamma = -i c for the oblate one, as in Flammer; and the relative size of its last term. The
    series of the second kind only converges for x > 1."""
    z = c * x
    f = _sw_sph(kf, L, z)
    s = ds = norm = t = mp.mpf(0)
    for j, l in enumerate(range(m + (n - m) % 2, L + 1, 2)):
        w = e[j] * math.perm(l + m, 2 * m)
        norm += w
        w *= (-1) ** ((l - n) // 2)
        t = w * f[l]
        s += t
        ds += w * (l / z * f[l] - f[l + 1])  # DLMF 10.51.2
    pre = (1 - kd / (x * x)) ** (mp.mpf(m) / 2) / norm
    return pre * s, pre * (c * ds + m * kd / (x * (x * x - kd)) * s), abs(t / s)


def _sw_taylor(m, c, kd, cv, x0, g0, dg0, x1):
    """Continue a solution g of p g'' + 2 (m + 1) x g' - (cv - m (m + 1) - c^2 x^2) g = 0 with
    p = x^2 - kd from x0 to x1 by its Taylor series at x0, which converges within the distance
    to the singularity (1, or +-i for kd = -1); returns g(x1) and g'(x1). The radial functions
    are p^(m/2) g."""
    h, a, c2 = x1 - x0, 2 * x0, c * c
    q, b, r = c2 * x0 * x0 + m * (m + 1) - cv, a * c2, -1 / (x0 * x0 - kd)
    u = [mp.mpf(0), mp.mpf(0), g0, dg0]  # u[k + 2] = the coefficient of (x - x0)^k
    g, dg, hk = g0 + dg0 * h, dg0, h  # hk = h^(k + 1)
    small, k = 0, 0
    while small < 5:
        t = (k + 1) * (k + m + 1) * a * u[k + 3] + (k * (k + 2 * m + 1) + q) * u[k + 2]
        t += b * u[k + 1] + c2 * u[k]
        u.append(r * t / ((k + 2) * (k + 1)))
        dt = (k + 2) * u[-1] * hk
        t = dt * h / (k + 2)
        g, dg, hk = g + t, dg + dt, hk * h
        # the last terms are negligible, by their binary exponents
        tiny = (
            mp.mag(dt) + mp.mag(h)
            < min(mp.mag(g), mp.mag(dg) + mp.mag(h)) - mp.mp.prec - 15
        )
        small = small + 1 if tiny else 0
        k += 1
    return g, dg


def _sw_rad(kf, m, n, c, x, kd):
    """The radial function of the first or second kind (kf), its derivative, and the
    characteristic value. The second kind for x < 2 is continued from its series at x = 2 with
    `_sw_taylor`, in steps of half the distance to the singularity. The extra digits are for the
    cancellation, which grows like e^c."""
    m, n = int(m), int(n)
    with wdps(10 + int(c) // 2):
        cv, e, L = _sw_cv(m, n, kd * c * c)
        if kf == 1:
            return *_sw_series(1, m, n, c, x, kd, e, L)[:2], cv
        x0 = max(x, mp.mpf(2))
        tol = mp.mpf(10) ** -(mp.mp.dps + 5)
        # the terms of the series decrease like x0^-l
        k = max(L - n, int((mp.mp.dps + 8) / mp.log10(x0))) // 2
        while True:
            L = n + 2 * k
            e, _ = _sw_coef(_sw_rows(m, n, kd * c * c, L), (n - m) // 2, cv)
            y, dy, last = _sw_series(2, m, n, c, x0, kd, e, L)
            if last <= tol:
                break
            k *= 2
        if x0 > x:
            pm = (x0 * x0 - kd) ** (mp.mpf(m) / 2)
            g = y / pm
            dg = dy / pm - m * x0 * g / (x0 * x0 - kd)
            while x0 > x:
                rho = x0 - 1 if kd == 1 else mp.sqrt(x0 * x0 + 1)
                x1 = max(x, x0 - rho / 2)
                g, dg = _sw_taylor(m, c, kd, cv, x0, g, dg, x1)
                x0 = x1
            pm = (x * x - kd) ** (mp.mpf(m) / 2)
            y, dy = pm * g, pm * (dg + m * x * g / (x * x - kd))
        return y, dy, cv


def _sw_radial(kf, m, n, c, x, kd):
    """The radial function of the first or second kind (kf) and its derivative; where they
    oscillate, with the modulus of the first and second kind as the scale."""
    y, dy, cv = _sw_rad(kf, m, n, c, x, kd)
    if (c * x) ** 2 <= cv + kd * m * m / (x * x - kd):
        return y, dy
    with mp.workdps(20):
        z, dz, _ = _sw_rad(3 - kf, m, n, c, x, kd)
    return (y, mp.hypot(y, z)), (dy, mp.hypot(dy, dz))


def _sw_rad1(m, n, c, x, kd):
    return _sw_radial(1, m, n, c, x, kd)


def _sw_rad2(m, n, c, x, kd):
    return _sw_radial(2, m, n, c, x, kd)


# ---------------------------------------------------------------- zeta functions


def r_gen_harmonic(n, a):
    n = int(n)
    with wdps(20):
        if n <= 10000:
            return mp.fsum(mp.mpf(k) ** -a for k in range(1, n + 1))
        return mp.zeta(a) - r_hurwitz_zeta(a, mp.mpf(n + 1))


def r_normalized_gen_harmonic(j, k, n, a):
    with wdps(20):
        p = [mp.mpf(i) ** -a for i in range(1, int(n) + 1)]
        return mp.fsum(p[int(j) - 1 : int(k)]) / mp.fsum(p)


# ---------------------------------------------------------------- exponential and trigonometric


def r_exp1(x):
    with wdps(xd(x)):
        return mp.e1(x)


def r_scaled_exp1(x):
    with wdps(xd(x)):
        return x * mp.exp(x) * mp.e1(x)


def r_si(x):
    with wdps(xd(x)):
        return mp.si(x)


def r_ci(x):
    # cephes returns Ci(|x|), the real part of Ci(x)
    with wdps(xd(x)):
        return mp.ci(abs(x))


def r_shi(x):
    with wdps(xd(x)):
        return mp.shi(x)


def r_chi(x):
    with wdps(xd(x)):
        return mp.chi(abs(x))


def r_fresnel_s(x):
    with wdps(2 * xd(x)):
        return mp.fresnels(x)


def r_fresnel_c(x):
    with wdps(2 * xd(x)):
        return mp.fresnelc(x)


# ---------------------------------------------------------------- inverse incomplete gamma, beta


def r_gammaincinv(a, p):
    return _invert_pos(lambda x: _gammainc_pq(a, x, False), p, a)


def r_gammainccinv(a, q):
    return _invert_pos(lambda x: _gammainc_pq(a, x, True), q, a, decreasing=True)


def r_betaincinv(a, b, y):
    return _invert_unit(lambda x, z: _betainc(a, b, x, z), y, a / (a + b))


def r_betainccinv(a, b, y):
    return _invert_unit(
        lambda x, z: _betainc(b, a, z, x), y, a / (a + b), decreasing=True
    )


# ---------------------------------------------------------------- elliptic


def r_ellipkm1(p):
    # K(1 - p) = pi / (2 agm(1, sqrt(p))), without the cancellation in 1 - p
    return mp.pi / (2 * mp.agm(1, mp.sqrt(p)))


def _converged(f, *a):
    """f(*a), with extra digits until two precision levels agree; mpmath's Carlson integrals
    can lose many digits (or overflow) for arguments of very different magnitudes."""
    prev = None
    for extra in (10, 40, 100, 300, 800):
        with wdps(extra):
            y = f(*a)
        if (
            prev is not None
            and mp.isfinite(y)
            and abs(y - prev) <= abs(y) * mp.mpf(10) ** (-mp.mp.dps - 5)
        ):
            return y
        prev = y
    raise ArithmeticError("no convergence")


def r_elliprc(x, y):
    return _converged(mp.elliprc, x, y)


def r_elliprd(x, y, z):
    return _converged(mp.elliprd, x, y, z)


def r_elliprf(x, y, z):
    return _converged(mp.elliprf, x, y, z)


def r_elliprg(x, y, z):
    return _converged(mp.elliprg, x, y, z)


def r_elliprj(x, y, z, p):
    return _converged(mp.elliprj, x, y, z, p)


def _ellipj(u, m):
    """sn, cn, dn, and the amplitude am(u | m) = n pi + atan2(sn(u'), cn(u')), where
    u' = u - 2 n K(m) in [-K, K), so that sn(u') = (-1)^n sn(u) and cn(u') = (-1)^n cn(u) >= 0."""
    with wdps(xd(u)):
        # complex with a zero imaginary part for m < 0
        sn, cn, dn = (mp.re(mp.ellipfun(k, u, m=m)) for k in ("sn", "cn", "dn"))
        if m == 1:
            return sn, cn, dn, mp.asin(mp.tanh(u))
        n = int(mp.floor((u + mp.ellipk(m)) / (2 * mp.ellipk(m))))
        s = (-1) ** n
        return sn, cn, dn, n * mp.pi + mp.atan2(s * sn, s * cn)


# ---------------------------------------------------------------- hypergeometric


def r_hyp0f0(z):
    return mp.exp(z)


def r_hyp1f0(a, z):
    with wdps(10):
        return (1 - z) ** -a


def r_hyp1f1(a, b, x):
    # For b < 0 the terms can dip far below the partial sum before growing again near n = -b,
    # where mpmath stops early, e.g. hyp1f1(2.5097844340670834, -424.1216104543897,
    # 127.87035219862966) is -1.635e24, but 0.516 at 80 digits. So add the depth of such a dip.
    extra = 20
    if b < 0 and b != int(b) and x != 0:
        logs, top, n = [0.0], 0.0, 0
        while n < 10**6:
            r = abs((a + n) * x / ((n + 1) * (b + n)))
            if r == 0 or (n > -b and r < 1 and b + n + 1 > abs(x)):
                break
            logs.append(logs[-1] + math.log10(r))
            n += 1
        peak = max(logs)
        last = max(i for i, v in enumerate(logs) if v > peak - mp.mp.dps - 20)
        for v in logs[: last + 1]:
            top = max(top, v)
            extra = max(extra, int(top - v) + 20)
    with wdps(extra):
        return mp.hyp1f1(a, b, x, maxterms=10**6)


def r_hypu(a, b, x):
    with wdps(20):
        return mp.hyperu(a, b, x)


def _hyp2f1_series(a, b, c, z):
    """The power series of hyp2f1 for |z| < 1 (in the working precision), and its largest term.

    For k > kmin, the term ratios z (a + k) / (k + 1) * (b + k) / (c + k) have factors of constant
    sign that tend monotonically to 1, so they are at most rho (below), which bounds the rest
    geometrically. Before that, e.g. for c < 0, the terms can be tiny before k = -c and huge after
    it."""
    s, t, top, k = mp.mpf(0), mp.mpf(1), mp.mpf(1), 0
    kmin = max(0, -a, -b, -c)
    while True:
        s += t
        top = max(top, abs(t))
        t *= (a + k) * (b + k) / ((c + k) * (k + 1)) * z
        k += 1
        if k > kmin:
            rho = abs(z) * max(1, (a + k) / (k + 1)) * max(1, (b + k) / (c + k))
            if rho < 1 and abs(t) * rho / (1 - rho) <= mp.eps * abs(s):
                return s + t, top
        if k > 10**6:
            raise NoConvergence


def r_hyp2f1(a, b, c, x):
    # mpmath's hyp2f1 can be far off for large parameters (mpmath/mpmath#1216), so for x <= 0.95
    # sum the series directly (for x < 0 the series of the Pfaff transformation in x / (x - 1),
    # DLMF 15.8.1), with 20 more digits than it loses to cancellation, up to 300 digits (else use
    # Arb).
    w = x / (x - 1) if x < 0 else x
    if w > 0.95 or a == int(a) <= 0 or b == int(b) <= 0 or c == int(c) <= 0:
        with wdps(20):
            return mp.hyp2f1(a, b, c, x, maxterms=10**6)
    extra = 20
    while True:
        with wdps(extra):
            if x < 0:
                # (c - b exactly, since e.g. c - b + 1 can be tiny)
                s, top = _hyp2f1_series(a, mp.fsub(c, b, exact=True), c, x / (x - 1))
            else:
                s, top = _hyp2f1_series(a, b, c, x)
            if s != 0 and top / abs(s) < mp.mpf(10) ** (extra - 15):
                return s * mp.exp(-a * mp.log1p(-x)) if x < 0 else s
        extra = int(mp.log10(top / abs(s))) + 20 if s != 0 else 2 * extra
        if mp.mp.dps + extra > 300:
            return _arb("hyp2f1", a, b, c, x)


# the most terms of `_pfq_series`; z = 1 - 2.2e-5 in the hyp3f2 domain needs about 5 million
PFQ_MAX_TERMS: Final = 10**7


def _pfq_series(a, b, z):
    """The power series of pFq for p = q + 1 as an Arb ball, for |z| < 1 or a polynomial; `b`
    contains the 1 of n!.

    For k > kmin, the factors (a_i + k) / (b_i + k) of the term ratios are positive and tend
    monotonically to 1, so the ratios are at most rho (below), which bounds the rest
    geometrically."""
    kmin = max(0, *(-v for v in (*a, *b)))
    a, b = ([flint.arb(v) for v in sorted(p)] for p in (a, b))
    z, tol = flint.arb(z), flint.arb(2) ** (-mp.mp.prec - 10)
    s, t, k = flint.arb(0), flint.arb(1), 0
    while not t.is_zero():
        s += t
        t *= z * math.prod(v + k for v in a) / math.prod(v + k for v in b)
        k += 1
        if k > kmin:
            rho = abs(z) * math.prod(((x + k) / (y + k)).max(1) for x, y in zip(a, b))
            rest = abs(t) * rho / (1 - rho)
            # (the upper bound of |s|, so that this also ends if s has no accuracy left)
            if rho < 1 and rest <= tol * abs(s).upper():
                return s + t + flint.arb(0, rest.upper())
        if k > PFQ_MAX_TERMS:
            raise NoConvergence
    return s


def r_hyppfq(*args):
    """pFq for p = q + 1, of a_1, ..., a_p, b_1, ..., b_q and z."""
    # mpmath's hyper can be far off for b < 0, like its hyp2f1 (mpmath/mpmath#1216), and Arb's
    # needs too many bits near |z| = 1, so sum the series directly, with more bits until it is
    # accurate. Only mpmath continues it to z <= -1, where it has to agree with itself at more
    # digits.
    p = len(args) // 2
    a, b, z = args[:p], args[p:-1], args[-1]
    if z <= -1 and not any(v == int(v) <= 0 for v in a):
        return _converged(mp.hyper, a, b, z)
    prec = 2 * mp.mp.prec
    while prec <= 2**14:
        with flint.ctx.workprec(prec):
            y = _pfq_series(a, (*b, 1), z)
        if y.rel_accuracy_bits() >= mp.mp.prec:
            return mp.mpf(y)
        prec *= 2
    raise NoConvergence


def r_hyp2f2(a1, a2, b1, b2, z):
    # Arb sums the series with a rigorous error bound, with more bits where its terms cancel
    return _arb("hyp2f2", a1, a2, b1, b2, z)


def r_hyp1f2(a, b1, b2, z):
    return _arb("hyp1f2", a, b1, b2, z)


def r_hyppfq_arb(name):
    """The reference of the spec `name` of pFq for p <= q or of 2F0, as for `r_hyp2f2`."""
    return lambda *args: _arb(name, *args)


# ---------------------------------------------------------------- orthogonal polynomials


def _poly(f, *a):
    with wdps(xd(*a[:-1]) + 20):
        return f(*a)


def r_jacobi(n, a, b, x):
    return _poly(mp.jacobi, n, a, b, x)


def r_sh_jacobi(n, a, b, x):
    return _poly(mp.jacobi, n, a, b, 2 * x - 1)


def r_gengegenbauer(n, p, q, x):
    with wdps(xd(n, p) + 20):
        return mp.jacobi(n, p - q, q - 1, 2 * x - 1) / mp.binomial(2 * n + p - 1, n)


def r_gegenbauer(n, a, x):
    return _poly(mp.gegenbauer, n, a, x)


def r_chebyt(n, x):
    return _poly(mp.chebyt, n, x)


def r_chebyu(n, x):
    return _poly(mp.chebyu, n, x)


def r_chebyc(n, x):
    return 2 * _poly(mp.chebyt, n, x / 2)


def r_chebys(n, x):
    return _poly(mp.chebyu, n, x / 2)


def r_sh_chebyt(n, x):
    return _poly(mp.chebyt, n, 2 * x - 1)


def r_sh_chebyu(n, x):
    return _poly(mp.chebyu, n, 2 * x - 1)


def r_legendre(n, x):
    return _poly(mp.legendre, n, x)


def r_sh_legendre(n, x):
    return _poly(mp.legendre, n, 2 * x - 1)


def r_genlaguerre(n, a, x):
    return _poly(mp.laguerre, n, a, x)


def r_laguerre(n, x):
    return _poly(mp.laguerre, n, 0, x)


def r_hermite_h(n, x):
    return _poly(mp.hermite, n, x)


def r_hermite_he(n, x):
    with wdps(xd(n) + 20):
        return 2 ** (-n / 2) * mp.hermite(n, x / mp.sqrt(2))


def r_hahn_q(n, big_n, a, b, x):
    """The terminating 3F2 sum (DLMF 18.20.5), in exact rational arithmetic."""
    a, b, x = (flint.fmpq(*float(v).as_integer_ratio()) for v in (a, b, x))
    n, big_n = int(n), int(big_n)
    s = flint.fmpq(1)
    for k in reversed(range(n)):
        s = 1 + s * (k - n) * (n + a + b + 1 + k) * (k - x) / (
            (a + 1 + k) * (k - big_n) * (k + 1)
        )
    return mp.mpf(int(s.p)) / int(s.q)


# ---------------------------------------------------------------- Legendre functions


def _legenp(n, m, x):
    """The Ferrers function P_n^m(x), with the Condon-Shortley phase, for integer m. For integer n
    and x < 0, use P_n^m(x) = (-1)^(n+m) P_n^m(-x), since mpmath may not converge near x = -1."""
    with wdps(xd(n) + 20):
        if x < 0 and n == int(n):
            return (-1) ** int(n + m) * mp.legenp(n, m, -x, type=2)
        return mp.legenp(n, m, x, type=2)


def r_assoc_legendre_p(n, m, x):
    return _legenp(n, m, x)


def r_assoc_legendre_p_norm(n, m, x):
    with wdps(xd(n) + 20):
        c = (2 * n + 1) / 2 * mp.factorial(n - m) / mp.factorial(n + m)
        return mp.sqrt(c) * _legenp(n, m, x)


def r_sph_legendre_p(n, m, theta):
    with wdps(xd(n) + 20):
        c = (2 * n + 1) / (4 * mp.pi) * mp.factorial(n - m) / mp.factorial(n + m)
        return mp.sqrt(c) * _legenp(n, m, mp.cos(theta))


def r_pmv(m, v, x):
    return _legenp(v, m, x)


def _legenq(n, m, x):
    """Q_n^m(x): the Ferrers function (mpmath type 2) for |x| < 1, and type 3 for x > 1. For
    x < 0, Q_n^m(x) = (-1)^(n+m+1) Q_n^m(-x) (DLMF 14.7.18), which for x < -1 is the real
    extension that `lqn` and `lqmn` compute."""
    with wdps(xd(n) + 20):
        sign = 1
        if x < 0:
            sign, x = (-1) ** int(n + m + 1), -x
        return sign * mp.legenq(n, m, x, type=2 if x < 1 else 3)


def r_legendre_q(n, x):
    return _legenq(n, 0, x)


def r_legendre_qp(n, x):
    return mp.diff(lambda t: _legenq(n, 0, t), x)


def r_assoc_legendre_q(n, m, x):
    return _legenq(n, m, x)


def r_assoc_legendre_qp(n, m, x):
    return mp.diff(lambda t: _legenq(n, m, t), x)


# ---------------------------------------------------------------- distributions


def r_bdtr(k, n, p):
    if k >= n:
        return mp.mpf(1)
    return _betainc(n - k, k + 1, 1 - p, p)


def r_bdtrc(k, n, p):
    if k >= n:
        return mp.mpf(0)
    return _betainc(k + 1, n - k, p, 1 - p)


def r_bdtri(k, n, y):
    # p such that bdtr(k, n, p) = y
    return _invert_unit(
        lambda p, q: _betainc(n - k, k + 1, q, p), y, (k + 1) / (n + 1), True
    )


def r_bdtrik(y, n, p):
    # real k in [0, n] such that bdtr(k, n, p) = y, or 0 if y lies below bdtr(0, n, p)
    if y <= (1 - p) ** n:
        return mp.mpf(0)
    return n * _invert_unit(lambda u, v: _betainc(n * v, n * u + 1, 1 - p, p), y, 0.5)


def r_bdtrin(k, y, p):
    # real n >= k such that bdtr(k, n, p) = y, found for n - k
    return k + _invert_pos(
        lambda nk: _betainc(nk, k + 1, 1 - p, p), y, 1, decreasing=True
    )


def r_btdtria(p, b, x):
    # a such that betainc(a, b, x) = p
    return _invert_pos(lambda a: _betainc(a, b, x, 1 - x), p, 1, decreasing=True)


def r_btdtrib(a, p, x):
    # b such that betainc(a, b, x) = p
    return _invert_pos(lambda b: _betainc(a, b, x, 1 - x), p, 1)


def r_nbdtr(k, n, p):
    return _betainc(n, k + 1, p, 1 - p)


def r_nbdtrc(k, n, p):
    return _betainc(k + 1, n, 1 - p, p)


def r_nbdtri(k, n, y):
    # p such that nbdtr(k, n, p) = y
    return _invert_unit(lambda p, q: _betainc(n, k + 1, p, q), y, n / (n + k + 1))


def r_nbdtrik(y, n, p):
    # real k >= 0 such that nbdtr(k, n, p) = y, or 0 if y lies below nbdtr(0, n, p)
    if y <= p**n:
        return mp.mpf(0)
    return _invert_pos(lambda k: _betainc(n, k + 1, p, 1 - p), y, 1)


def r_nbdtrin(k, y, p):
    # real n such that nbdtr(k, n, p) = y
    return _invert_pos(lambda n: _betainc(n, k + 1, p, 1 - p), y, 1, decreasing=True)


def r_pdtr(k, m):
    return _gammainc_pq(k + 1, m, True)


def r_pdtrc(k, m):
    return _gammainc_pq(k + 1, m, False)


def r_pdtri(k, y):
    # m such that pdtr(k, m) = y
    return _invert_pos(
        lambda m: _gammainc_pq(k + 1, m, True), y, k + 1, decreasing=True
    )


def r_pdtrik(p, m):
    # real k >= 0 such that pdtr(k, m) = p, or 0 if p lies below pdtr(0, m)
    if p <= mp.exp(-m):
        return mp.mpf(0)
    return _invert_pos(lambda k: _gammainc_pq(k + 1, m, True), p, m + 1)


def r_gdtr(a, b, x):
    return _gammainc_pq(b, a * x, False)


def r_gdtrc(a, b, x):
    return _gammainc_pq(b, a * x, True)


def r_gdtria(p, b, x):
    # a such that gdtr(a, b, x) = p
    return _invert_pos(lambda a: _gammainc_pq(b, a * x, False), p, b / x)


def r_gdtrib(a, p, x):
    # b such that gdtr(a, b, x) = p
    return _invert_pos(
        lambda b: _gammainc_pq(b, a * x, False), p, a * x, decreasing=True
    )


def r_gdtrix(a, b, p):
    # x such that gdtr(a, b, x) = p
    return _invert_pos(lambda x: _gammainc_pq(b, a * x, False), p, b / a)


def r_chdtr(v, x):
    return _gammainc_pq(v / 2, x / 2, False)


def r_chdtrc(v, x):
    return _gammainc_pq(v / 2, x / 2, True)


def r_chdtri(v, y):
    # x such that chdtrc(v, x) = y
    return _invert_pos(
        lambda x: _gammainc_pq(v / 2, x / 2, True), y, v, decreasing=True
    )


def r_chdtriv(p, x):
    # v such that chdtr(v, x) = p
    return _invert_pos(
        lambda v: _gammainc_pq(v / 2, x / 2, False), p, x, decreasing=True
    )


def _chndtr_series(x, v, nc):
    """sum_j Pois(j; nc/2) P(v/2 + j, x/2) (A&S 26.4.25), from the largest Poisson weight both
    ways, with P(a, y) - P(a + 1, y) = y^a e^-y / Gamma(a + 1) (A&S 6.5.21)."""
    if x == 0:
        return mp.mpf(0)
    lam, y = mp.mpf(nc) / 2, mp.mpf(x) / 2
    j0 = int(lam)
    a0 = mp.mpf(v) / 2 + j0
    w0 = mp.exp(j0 * mp.log(lam) - lam - mp.loggamma(j0 + 1))
    p0 = _gammainc_pq(a0, y, False)
    t0 = mp.exp(a0 * mp.log(y) - y - mp.loggamma(a0 + 1))
    s, tol = w0 * p0, mp.eps / 2**10
    w, p, a, t, j = w0, p0, a0, t0, j0
    while j > 0:
        w, t, j, a = w * j / lam, t * a / y, j - 1, a - 1
        p += t
        s += w * p
        if w * p < tol * s:
            break
    w, p, a, t, j = w0, p0, a0, t0, j0
    while True:
        p -= t
        w, a, j = w * lam / (j + 1), a + 1, j + 1
        t *= y / a
        s += w * p
        if j > lam and w * p < tol * s:
            return s


def r_chndtr(x, v, nc):
    # the forward recurrence cancels, and the other way loses digits to the terms' exponents
    return _converged(_chndtr_series, x, v, nc)


def r_chndtri(p, v, nc):
    # x such that chndtr(x, v, nc) = p
    return _invert_pos(lambda x: _chndtr_series(x, v, nc), p, v + nc)


def r_chndtridf(x, p, nc):
    # v such that chndtr(x, v, nc) = p
    return _invert_pos(
        lambda v: _chndtr_series(x, v, nc), p, max(x - nc, 1), decreasing=True
    )


def r_chndtrinc(x, v, p):
    # nc such that chndtr(x, v, nc) = p
    return _invert_pos(
        lambda nc: _chndtr_series(x, v, nc), p, max(x - v, 1), decreasing=True
    )


def _ncfdtr_series(a, b, nc, x, k=0):
    """sum_j Pois(j + k; nc/2) I_y(a/2 + j, b/2) (A&S 26.6.20 if k = 0) with y = a x / (b + a x)
    and Pois(k; l) = e^-l l^k / Gamma(k + 1), from the largest Poisson weight both ways, with
    I_y(c, d) - I_y(c + 1, d) = y^c (1 - y)^d / (c B(c, d)) (A&S 26.5.16)."""
    a, b, x = mp.mpf(a), mp.mpf(b), mp.mpf(x)
    lam, y, z, d = mp.mpf(nc) / 2, a * x / (b + a * x), b / (b + a * x), b / 2
    j0 = int(lam)
    c0 = a / 2 + j0
    w0 = (
        mp.exp((j0 + k) * mp.log(lam) - lam - mp.loggamma(j0 + k + 1))
        if lam
        else mp.mpf(1)
    )
    i0 = _betainc(c0, d, y, z)
    t0 = mp.exp(c0 * mp.log(y) + d * mp.log(z) - mp.log(c0) - mp.log(mp.beta(c0, d)))
    s, tol = w0 * i0, mp.eps / 2**10
    w, i, c, t, j = w0, i0, c0, t0, j0
    while j > 0:
        t *= c / (y * (c + d - 1))
        w, j, c = w * (j + k) / lam, j - 1, c - 1
        i += t
        s += w * i
        if w * i < tol * s:
            break
    w, i, c, t, j = w0, i0, c0, t0, j0
    while True:
        i -= t
        t *= y * (c + d) / (c + 1)
        w, c, j = w * lam / (j + k + 1), c + 1, j + 1
        s += w * i
        if j > lam and w * i < tol * s:
            return s


def r_ncfdtr(a, b, nc, x):
    if x == 0:
        return mp.mpf(0)
    return _converged(_ncfdtr_series, a, b, nc, x)


def r_ncfdtri(a, b, nc, p):
    # x such that ncfdtr(a, b, nc, x) = p
    return _invert_pos(lambda x: _ncfdtr_series(a, b, nc, x), p, 1)


def r_fdtri(a, b, y):
    # x such that fdtr(a, b, x) = y
    return _invert_pos(lambda x: r_fdtr(a, b, x), y, 1)


def r_fdtridfd(a, p, x):
    # b such that fdtr(a, b, x) = p
    return _only_root_pos(lambda b: r_fdtr(a, b, x), p)


def r_ncfdtridfd(a, p, nc, x):
    # b such that ncfdtr(a, b, nc, x) = p
    return _only_root_pos(lambda b: _ncfdtr_series(a, b, nc, x), p)


def r_ncfdtridfn(p, b, nc, x):
    # a such that ncfdtr(a, b, nc, x) = p
    return _only_root_pos(lambda a: _ncfdtr_series(a, b, nc, x), p)


def r_ncfdtrinc(a, b, p, x):
    # nc such that ncfdtr(a, b, nc, x) = p, which decreases in nc
    return _invert_pos(lambda nc: _ncfdtr_series(a, b, nc, x), p, 1, decreasing=True)


def _stdtr_tail(v, t):
    """P(T <= -|t|) = I_x(v/2, 1/2) / 2, with x = v / (v + t^2)."""
    return _betainc(v / 2, mp.mpf(1) / 2, v / (v + t * t), t * t / (v + t * t)) / 2


def r_stdtr(v, t):
    with wdps(xd(v, t)):
        tail = _stdtr_tail(v, t)
        return tail if t < 0 else 1 - tail


def _nctdtr_sum(v, nc, t):
    """P(T <= t) = Phi(-nc) + (sgn(t) S_0 + sgn(nc) S_1/2) / 2, with
    S_k = sum_j Pois(j + k; nc^2/2) I_x(j + k + 1/2, v/2) and x = t^2 / (v + t^2), where S_0 is
    the noncentral F CDF of T^2."""
    if t == 0:
        return mp.ncdf(-nc)
    v, nc, t = mp.mpf(v), mp.mpf(nc), mp.mpf(t)
    s = mp.sign(t) * _ncfdtr_series(1, v, nc * nc, t * t)
    if nc != 0:
        s += mp.sign(nc) * _ncfdtr_series(2, v, nc * nc, t * t / 2, mp.mpf(1) / 2)
    return mp.ncdf(-nc) + s / 2


def r_nctdtr(v, nc, t):
    # cancels if t < 0 < nc
    return _converged(_nctdtr_sum, v, nc, t)


def r_nctdtrit(v, nc, p):
    # in `asinh(t)`, since the tails can be as heavy as `|t|^-v`
    return mp.sinh(_invert(lambda s: r_nctdtr(v, nc, mp.sinh(s)) - p, mp.asinh(nc)))


def r_stdtridf(p, t):
    # v such that stdtr(v, t) = p, via the tail P(T <= -|t|), which is p or 1 - p
    return _only_root_pos(lambda v: _stdtr_tail(v, t), p if t < 0 else 1 - p)


def r_nctdtridf(p, nc, t):
    # v such that nctdtr(v, nc, t) = p
    return _only_root_pos(lambda v: r_nctdtr(v, nc, t), p)


def r_nctdtrinc(v, p, t):
    # nc such that nctdtr(v, nc, t) = p, which decreases from 1 to 0 in nc
    return _invert(lambda nc: p - r_nctdtr(v, nc, t), t)


def r_stdtri(v, p):
    if p == 0.5:
        return mp.mpf(0)
    q = min(p, 1 - p)
    t = _invert_pos(lambda t: _stdtr_tail(v, t), q, 1, decreasing=True)
    return -t if p < 0.5 else t


def r_tukeylambdacdf(x, lmbda):
    # p such that the quantile function (p^lambda - (1 - p)^lambda) / lambda is x
    if lmbda == 0:
        return 1 / (1 + mp.exp(-x))
    if lmbda > 0 and abs(x) >= 1 / lmbda:
        return mp.mpf(1 if x > 0 else 0)
    return _invert_unit(lambda p, q: (p**lmbda - q**lmbda) / lmbda, x, 0.5)


def r_von_mises_cdf(k, x):
    # int_-pi^x exp(k (cos t - c)) dt exp(k (c - 1)) / (2 pi I_0(k) e^-k), with the maximum
    # c = cos(min(x, 0)) of the integrand factored out, since the tail can underflow the
    # absolute tolerance of the quadrature
    with wdps(10):
        c = mp.cos(min(x, 0))
        f = lambda t: mp.exp(k * (mp.cos(t) - c))
        w = 1 / mp.sqrt(k)
        pts = (
            [-mp.pi] + [t for t in (-10 * w, -w, 0, w, 10 * w) if -mp.pi < t < x] + [x]
        )
        return (
            mp.quad(f, pts)
            * mp.exp(k * (c - 1))
            / (2 * mp.pi * mp.besseli(0, k) * mp.exp(-k))
        )


def _geninvgauss_logpdf(x, p, b):
    return (p - 1) * mp.log(x) - mp.log(2 * mp.besselk(p, b)) - b * (x + 1 / x) / 2


def r_geninvgauss_pdf(x, p, b):
    with wdps(xd(x, b)):
        return mp.exp(_geninvgauss_logpdf(x, p, b))


def r_geninvgauss_logpdf(x, p, b):
    with wdps(xd(x, b)):
        return _geninvgauss_logpdf(x, p, b)


def _genhyperbolic_logpdf(x, p, a, b):
    g, s = mp.sqrt(a * a - b * b), mp.sqrt(1 + x * x)
    c = p / 2 * mp.log(g * g) - mp.log(mp.sqrt(2 * mp.pi)) - (p - 0.5) * mp.log(a)
    c -= mp.log(mp.besselk(p, g))
    return c + b * x + mp.log(mp.besselk(p - 0.5, a * s)) + (p - 0.5) * mp.log(s)


def r_genhyperbolic_pdf(x, p, a, b):
    with wdps(xd(x, a)):
        return mp.exp(_genhyperbolic_logpdf(x, p, a, b))


def r_genhyperbolic_logpdf(x, p, a, b):
    with wdps(xd(x, a)):
        return _genhyperbolic_logpdf(x, p, a, b)


def r_cdf_cvm_inf(x):
    # Csorgo and Faraway (1996), eq. 1.2
    with wdps(10):
        s, k = mp.mpf(0), 0
        while True:
            y = (4 * k + 1) ** 2 / (16 * x)
            t = (
                mp.gamma(k + 0.5)
                / (mp.gamma(0.5) * mp.factorial(k))
                * mp.sqrt(4 * k + 1)
            )
            t *= mp.exp(-y) * mp.besselk(0.25, y)
            s += t
            if t < mp.eps * s:
                return s / (mp.pi * mp.sqrt(x))
            k += 1


# ---------------------------------------------------------------- Kolmogorov-Smirnov


def _smirnov_terms(n, x, js):
    """The terms binom(n, j) (x + j/n)^(j - 1) (1 - x - j/n)^(n - j) of the Birnbaum-Tingey sum,
    for j in js."""
    return [
        mp.binomial(n, j)
        * (x + mp.mpf(j) / n) ** (j - 1)
        * (mp.mpf(n - j) / n - x) ** (n - j)
        for j in js
    ]


def _smirnov_jmax(n, x):
    """floor(n (1 - x)), exactly, also for an x with more digits than an f64."""
    return math.floor(n * (1 - Fraction(*x.as_integer_ratio())))


def _smirnov(n, x):
    """P(D_n^+ >= x), with the exact Birnbaum-Tingey sum (all its terms are positive)."""
    n = int(n)
    return x * mp.fsum(_smirnov_terms(n, x, range(_smirnov_jmax(n, x) + 1)))


def _smirnovc(n, x):
    """1 - smirnov(n, x). Where that cancels, use the rest of the sum, since the full sum over
    j = 0, ..., n is 1 (Abel's identity); its terms alternate, but there are only about nx."""
    c = 1 - _smirnov(n, x)
    if c > 1e-3:
        return c
    n = int(n)
    with wdps(40):
        return x * mp.fsum(_smirnov_terms(n, x, range(_smirnov_jmax(n, x) + 1, n + 1)))


def r_smirnov(n, x):
    with wdps(20):
        return _smirnov(n, x)


def r_smirnovc(n, x):
    with wdps(20):
        return _smirnovc(n, x)


def r_smirnovp(n, x):
    # the derivative of the sum of smirnov, term by term, with a = x + j/n and b = 1 - a; the
    # term j = 0 is (1 - x)^n
    n = int(n)
    with wdps(20):
        d = -n * (1 - x) ** (n - 1)
        for j in range(1, _smirnov_jmax(n, x) + 1):
            a, b = x + mp.mpf(j) / n, mp.mpf(n - j) / n - x
            t = a ** (j - 1) * b ** (n - j) + x * (j - 1) * a ** (j - 2) * b ** (n - j)
            d += mp.binomial(n, j) * (t - x * (n - j) * a ** (j - 1) * b ** (n - j - 1))
        return d


def r_smirnovi(n, q):
    return _invert_unit(lambda x, _: _smirnov(n, x), q, 0.5 / mp.sqrt(n), True)


def r_smirnovci(n, p):
    return _invert_unit(lambda x, _: _smirnovc(n, x), p, 0.5 / mp.sqrt(n))


def _kolmogorov(x, deriv=False):
    """P(K >= x) = 2 sum_k>=1 (-1)^(k-1) exp(-2 k^2 x^2), or, for x < 1, its complement
    sqrt(2 pi) / x sum_k>=1 exp(-(2k - 1)^2 pi^2 / (8 x^2)); returns (sf, cdf), or their
    derivatives."""
    if x >= 1:
        s = mp.nsum(lambda k: (-1) ** (k - 1) * mp.exp(-2 * k * k * x * x), [1, mp.inf])
        if deriv:
            s = mp.nsum(
                lambda k: (-1) ** k * 4 * k * k * x * mp.exp(-2 * k * k * x * x),
                [1, mp.inf],
            )
        return 2 * s, 1 - 2 * s if not deriv else -2 * s
    c = mp.sqrt(2 * mp.pi)

    def a(k):
        return (2 * k - 1) ** 2 * mp.pi**2 / 8

    s = c / x * mp.nsum(lambda k: mp.exp(-a(k) / x**2), [1, mp.inf])
    if deriv:
        s = (
            mp.nsum(
                lambda k: mp.exp(-a(k) / x**2) * (2 * a(k) / x**4 - 1 / x**2),
                [1, mp.inf],
            )
            * c
        )
    return (1 - s, s) if not deriv else (-s, s)


def r_kolmogorov(x):
    with wdps(20):
        return _kolmogorov(x)[0]


def r_kolmogc(x):
    with wdps(20):
        return _kolmogorov(x)[1]


def r_kolmogp(x):
    with wdps(20):
        return _kolmogorov(x, True)[0]


def r_kolmogi(q):
    return _invert_pos(lambda x: _kolmogorov(x)[0], q, 1, decreasing=True)


def r_kolmogci(p):
    return _invert_pos(lambda x: _kolmogorov(x)[1], p, 1)


# ---------------------------------------------------------------- other


def _wright_bessel(a, b, x):
    """sum_k x^k / (k! Gamma(a k + b)); all terms are nonnegative for a, b, x >= 0."""
    s, t, k, peaked, last = mp.mpf(0), mp.mpf(1), 0, False, mp.mpf(0)
    while True:
        term = t * mp.rgamma(a * k + b)
        s += term
        if k > 0 and term <= last:
            peaked = True
        if peaked and term <= mp.eps * s:
            return s
        last = term
        k += 1
        t *= x / k


def r_wright_bessel(a, b, x):
    with wdps(10):
        return _wright_bessel(a, b, x)


def r_log_wright_bessel(a, b, x):
    with wdps(10):
        return mp.log(_wright_bessel(a, b, x))


def r_wrightomega(x):
    # the solution w of w + log(w) = x, by Newton's method
    with wdps(10):
        w = x - mp.log(x) if x > 1 else mp.exp(x)
        for _ in range(1000):
            dw = (w + mp.log(w) - x) / (1 + 1 / w)
            w -= dw
            if abs(dw) <= mp.eps * abs(w):
                break
        return w


# ---------------------------------------------------------------- complex


def _c(f):
    """A reference for f(z), with z = re + i im."""

    def ref(re, im):
        with wdps(xd(re, im)):
            return f(mp.mpc(re, im))

    ref.__name__ = f"c_{f.__name__}"
    return ref


def c_wofz(z):
    return mp.exp(-z * z) * mp.erfc(-1j * z)


def c_lambertw(z):
    return mp.lambertw(z)


def c_hyp1f1(a, b, re, im):
    with wdps(20):
        return mp.hyp1f1(a, b, mp.mpc(re, im), maxterms=10**6)


def c_hyp2f1(a, b, c, re, im):
    with wdps(20):
        return mp.hyp2f1(a, b, c, mp.mpc(re, im), maxterms=10**6)


def _cbessel(f, v, re, im, osc=False):
    with wdps(xd(v, re, im)):
        z = mp.mpc(re, im)
        y = f(v, z)
        if osc:
            return y, mp.sqrt(abs(mp.besselj(v, z)) ** 2 + abs(mp.bessely(v, z)) ** 2)
        return y


def c_besselj(v, re, im):
    return _cbessel(mp.besselj, v, re, im, True)


def c_bessely(v, re, im):
    return _cbessel(mp.bessely, v, re, im, True)


def c_besseli(v, re, im):
    return _cbessel(mp.besseli, v, re, im)


def c_besselk(v, re, im):
    return _cbessel(mp.besselk, v, re, im)


def c_hankel1(v, re, im):
    return _cbessel(mp.hankel1, v, re, im)


def c_hankel2(v, re, im):
    return _cbessel(mp.hankel2, v, re, im)


def c_hankel1e(v, re, im):
    # H1_v(z) e^(-iz)
    return _cbessel(lambda v, z: mp.hankel1(v, z) * mp.exp(-1j * z), v, re, im)


def c_hankel2e(v, re, im):
    # H2_v(z) e^(iz)
    return _cbessel(lambda v, z: mp.hankel2(v, z) * mp.exp(1j * z), v, re, im)


def _hankel_prime(f, v, re, im, n):
    """The n-th derivative 2^-n sum_k (-1)^k binom(n, k) H_(v - n + 2k)(z), DLMF 10.6.7."""
    n = int(n)
    terms = [
        (-1) ** k * mp.binomial(n, k) * f(v - n + 2 * k, mp.mpc(re, im))
        for k in range(n + 1)
    ]
    return mp.fsum(terms) / 2**n


def c_hankel1p(v, re, im, n):
    with wdps(xd(v, re, im) + 10):
        return _hankel_prime(mp.hankel1, v, re, im, n)


def c_hankel2p(v, re, im, n):
    with wdps(xd(v, re, im) + 10):
        return _hankel_prime(mp.hankel2, v, re, im, n)


def c_modified_fresnel(x):
    """F+, K+, F-, K-, with F+(x) = int_x^inf exp(i t^2) dt
    = sqrt(pi/2) (1/2 - C(u) + i (1/2 - S(u))), u = x sqrt(2/pi), and F- its conjugate;
    K+(x) = exp(-i (x^2 + pi/4)) F+(x) / sqrt(pi), and K-(x) its conjugate."""
    with wdps(2 * xd(x)):
        u = x * mp.sqrt(2 / mp.pi)
        fp = mp.sqrt(mp.pi / 2) * mp.mpc(0.5 - mp.fresnelc(u), 0.5 - mp.fresnels(u))
        kp = mp.exp(-1j * (x * x + mp.pi / 4)) * fp / mp.sqrt(mp.pi)
        return fp, kp, mp.conj(fp), mp.conj(kp)


def c_sph_harm_y(n, m, theta, phi):
    with wdps(xd(n) + 20):
        return mp.spherharm(int(n), int(m), theta, phi)


# ---------------------------------------------------------------- functions of arrays
# The scalar arguments come first, followed by the elements of the array.


def _exact(x):
    return Fraction(*x.as_integer_ratio())


def r_logsumexp(*xs):
    """ln(sum e^x); near 0, ln(1 + (e^m - 1) + the sum of the other e^x) with m = max x, whose
    two terms are accurate, so that they cancel only once (unlike the e^x - 1 of all x)."""
    with wdps(20):
        s = mp.fsum(map(mp.exp, xs))
        if not 0.5 < s < 2:
            return mp.log(s)
        i = xs.index(max(xs))
        rest = mp.fsum(mp.exp(x) for k, x in enumerate(xs) if k != i)
        return mp.log1p(mp.expm1(xs[i]) + rest)


def _softmax_sum(j, xs):
    """sum e^(x_k - x_j) over k != j, with the exact differences, and without the x_k = -inf."""
    j = int(j)
    return mp.fsum(
        mp.exp(mp.fsub(x, xs[j], exact=True))
        for k, x in enumerate(xs)
        if k != j and x != -mp.inf
    )


def r_softmax(j, *xs):
    with wdps(20):
        return 1 / (1 + _softmax_sum(j, xs))


def r_log_softmax(j, *xs):
    with wdps(20):
        return -mp.log1p(_softmax_sum(j, xs))


def c_cevalpoly(re, im, *c):
    """The exact value, with Horner's rule in Gaussian rationals; `c` starts at the highest
    degree."""
    zr, zi = _exact(re), _exact(im)
    yr = yi = Fraction(0)
    for ck in map(_exact, c):
        yr, yi = yr * zr - yi * zi + ck, yr * zi + yi * zr
    return mp.mpc(
        mp.fdiv(yr.numerator, yr.denominator), mp.fdiv(yi.numerator, yi.denominator)
    )


def _poisson_binom(p):
    """The exact pmf of the number of successes as integers c_k / d: the coefficients of
    prod (1 - p_i + p_i x), with the p_i = a_i / 2^e on a common denominator 2^e."""
    ps = list(map(_exact, p))
    e = max(f.denominator for f in ps)
    c = [1]
    for f in ps:
        a = f.numerator * (e // f.denominator)
        c = [x * (e - a) + y * a for x, y in zip([*c, 0], [0, *c], strict=True)]
    return c, e ** len(ps)


def r_poisson_binom_pmf(k, *p):
    c, d = _poisson_binom(p)
    return mp.fdiv(c[int(k)], d)


def r_poisson_binom_cdf(k, *p):
    c, d = _poisson_binom(p)
    return mp.fdiv(sum(c[: int(k) + 1]), d)


# the numbers of subsets of {1, ..., n} by their sum, for n = 0, 1, ...
_WILCOXON_COUNTS = [[1]]


def _wilcoxon_counts(n):
    """The number of subsets of {1, ..., n} with the sum k, for k = 0, ..., n (n + 1) / 2: the
    coefficients of prod (1 + x^j), so 2^n times the pmf of the signed-rank statistic."""
    while len(_WILCOXON_COUNTS) <= n:
        m = len(_WILCOXON_COUNTS)
        c = [*_WILCOXON_COUNTS[-1], *[0] * m]
        _WILCOXON_COUNTS.append(
            [x + (c[k - m] if k >= m else 0) for k, x in enumerate(c)]
        )
    return _WILCOXON_COUNTS[n]


def r_wilcoxon_pmf(n, k):
    return mp.ldexp(_wilcoxon_counts(int(n))[int(k)], -int(n))


def r_wilcoxon_cdf(n, k):
    return mp.ldexp(sum(_wilcoxon_counts(int(n))[: int(k) + 1]), -int(n))


def r_wilcoxon_sf(n, k):
    """P(T >= k)."""
    return mp.ldexp(sum(_wilcoxon_counts(int(n))[int(k) :]), -int(n))


def _osc_cz(y, v, z, s=1):
    """`osc_cz`: y with the scale s sqrt(|J_v(z)|^2 + |Y_v(z)|^2) for |z| > |v|."""
    if abs(z) > abs(v):
        return y, s * mp.sqrt(abs(mp.besselj(v, z)) ** 2 + abs(mp.bessely(v, z)) ** 2)
    return y


def c_besselje(v, re, im):
    # J_v(z) e^-|Im z|
    with wdps(xd(v, re, im)):
        z, s = mp.mpc(re, im), mp.exp(-abs(im))
        return _osc_cz(mp.besselj(v, z) * s, v, z, s)


def c_besselye(v, re, im):
    # Y_v(z) e^-|Im z|
    with wdps(xd(v, re, im)):
        z, s = mp.mpc(re, im), mp.exp(-abs(im))
        return _osc_cz(mp.bessely(v, z) * s, v, z, s)


def c_besselie(v, re, im):
    # I_v(z) e^-|Re z|
    return _cbessel(lambda v, z: mp.besseli(v, z) * mp.exp(-abs(z.real)), v, re, im)


def c_besselke(v, re, im):
    # K_v(z) e^z
    return _cbessel(lambda v, z: mp.besselk(v, z) * mp.exp(z), v, re, im)


def c_jvp(v, re, im, n):
    with wdps(xd(v, re, im)):
        z = mp.mpc(re, im)
        return _osc_cz(_bessel_diff(mp.besselj, v, z, n, -1), v, z)


def c_yvp(v, re, im, n):
    with wdps(xd(v, re, im)):
        z = mp.mpc(re, im)
        return _osc_cz(_bessel_diff(mp.bessely, v, z, n, -1), v, z)


def c_ivp(v, re, im, n):
    with wdps(xd(v, re, im)):
        return _bessel_diff(mp.besseli, v, mp.mpc(re, im), n, 1)


def c_kvp(v, re, im, n):
    with wdps(xd(v, re, im)):
        return _bessel_diff(mp.besselk, v, mp.mpc(re, im), n, 1)


def c_airy(re, im, scaled=False):
    """Ai, Ai', Bi, Bi'; scaled like `airy_scaled`, i.e. Ai and Ai' by e^zeta, and Bi and Bi' by
    e^-|Re zeta|, with zeta = 2/3 z^(3/2)."""
    with wdps(xd(re, im) + 10):
        z = mp.mpc(re, im)
        ai, aip = mp.airyai(z), mp.airyai(z, derivative=1)
        bi, bip = mp.airybi(z), mp.airybi(z, derivative=1)
        if not scaled:
            return [ai, aip, bi, bip]
        zeta = 2 * z * mp.sqrt(z) / 3
        ea, eb = mp.exp(zeta), mp.exp(-abs(zeta.real))
        return [ai * ea, aip * ea, bi * eb, bip * eb]


def c_erfcx(z):
    return mp.exp(z * z) * mp.erfc(z)


def c_dawsn(z):
    return mp.sqrt(mp.pi) / 2 * mp.exp(-z * z) * mp.erfi(z)


def c_ndtr(z):
    return mp.erfc(-z / mp.sqrt(2)) / 2


def c_log_ndtr(z):
    # the principal branch; log1p for Re z > 0, where ndtr(z) can be close to 1
    if z.real > 0:
        return mp.log1p(-mp.erfc(z / mp.sqrt(2)) / 2)
    return mp.log(c_ndtr(z))


def c_fresnel(re, im):
    with wdps(2 * xd(re, im)):
        z = mp.mpc(re, im)
        return mp.fresnels(z), mp.fresnelc(z)


def c_hyp1f0(a, re, im):
    with wdps(10):
        return (1 - mp.mpc(re, im)) ** -a


def c_hyp0f1(b, re, im):
    # Arb, since mpmath's hyp0f1 can be far off for b < 0 (mpmath/mpmath#1216), with increasing
    # precision until it is accurate to the working precision. For 4|z| > (b - 1)^2, the scale of
    # 0F1(; b; z) = Gamma(b) (-z)^(-v/2) J_v(x), with v = b - 1 and x = 2 sqrt(-z).
    prec = mp.mp.prec
    while True:
        with flint.ctx.workprec(prec):
            y = flint.acb(re, im).hypgeom_0f1(b)
        if y.rel_accuracy_bits() >= mp.mp.prec:
            break
        if prec > 2**14:
            raise NoConvergence
        prec *= 2
    z = mp.mpc(re, im)
    if 4 * abs(z) <= (b - 1) ** 2:
        return mp.mpc(y)
    with wdps(xd(b, re, im)):
        v, x = b - 1, 2 * mp.sqrt(-z)
        env = mp.sqrt(abs(mp.besselj(v, x)) ** 2 + abs(mp.bessely(v, x)) ** 2)
        return mp.mpc(y), abs(mp.gamma(b)) * abs(z) ** (-v / 2) * env


# ---------------------------------------------------------------- zeros and Gauss quadrature


def _airy_zeros(i, k, bi=False):
    """Output i of `airy_ai_zeros` (or `airy_bi_zeros`) for zero k: the zeros of Ai and Ai', the
    value of Ai at the zero of Ai', and the value of Ai' at the zero of Ai."""
    zero, f = (mp.airybizero, mp.airybi) if bi else (mp.airyaizero, mp.airyai)
    k = int(k)
    if i < 2:
        return zero(k, i)
    return f(zero(k, 1)) if i == 2 else f(zero(k), derivative=1)


def _bessel_zero(i, v, k):
    """Output i of `bessel_zeros`: zero k of J_v, J_v', Y_v, or Y_v'. Unlike `bessel_zeros`,
    mpmath counts x = 0 as the first zero of J_0'."""
    v, k = int(v), int(k)
    if i < 2:
        return mp.besseljzero(v, k + (i == 1 and v == 0), i)
    return mp.besselyzero(v, k, i - 2)


def _y_zero(i, kind, c, k):
    """Output i of `bessel_y0_zeros`, `bessel_y1_zeros`, or `bessel_y1_prime_zeros` (kind 0 to
    2): zero k of Y_0, Y_1, or Y_1', or the value there of Y_1, Y_0, or Y_1.

    The complex zeros (c = 1) are zeros of Y_n(-z) = (-1)^n (Y_n(z) + 2i J_n(z)) (DLMF 10.11.6),
    so by the leading terms of DLMF 10.17.3 and 10.17.4, they are near
    -(k - (2 kind + 1) / 4) pi + i ln(3) / 2, except for the first ones of Y_1 and Y_1'. They are
    found by Newton's method from there, and must not move by more than 1 (the spacing is pi)."""
    n, d, k = int(kind == 1), int(kind == 2), int(k)
    if not c:
        z = mp.besselyzero(n + d, k, d)
    else:
        if kind and k == 1:
            z0 = mp.mpc(-0.5, 0.8) if kind == 1 else mp.mpc(0.6, 0.9)
        else:
            z0 = mp.mpc(-(k - mp.mpf(2 * kind + 1) / 4) * mp.pi, mp.log(3) / 2)
        z = mp.findroot(lambda z: mp.bessely(n + d, z, d), z0)
        if abs(z - z0) > 1:
            raise ArithmeticError(f"zero {z} is not near {z0}")
    return z if i == 0 else mp.bessely(1 - n, z)


# the zeros of `_jnjnp_zeros`, computed once per process
_JNJNP: list[tuple[Any, int, int, int]] = []


def _jnjnp_zeros():
    """(x, n, m, p) of the zeros x of all J_n (p = 0) and J_n' (p = 1) below 72, in ascending
    order, where x is zero m of J_n or J_n', counting x = 0 as zero 0 of J_0'.

    The zeros are the sign changes on a grid with step 1/4, which is less than their spacing
    (> 2), polished with `_solve` at 30 digits. The orders n < 72 have all of them, since the
    zeros of J_n and J_n' are greater than n (DLMF 10.21.3). J_0' = -J_1 (DLMF 10.6.3), so they
    share their zeros, which are listed as those of J_1 first."""
    if not _JNJNP:
        big_x = 72
        with mp.workdps(30):
            zeros = [(mp.mpf(0), 0, 0, 1)]
            for n in range(big_x):
                for p in range(int(n > 0) + 1):
                    m = 0
                    f = lambda x, n=n, p=p: mp.besselj(n, x, p)
                    a, fa = mp.mpf(1) / 4, f(mp.mpf(1) / 4)
                    while a < big_x:
                        b = a + mp.mpf(1) / 4
                        fb = f(b)
                        if (fa < 0) != (fb < 0):
                            m += 1
                            zeros.append((_solve(f, a, b, fa, fb), n, m, p))
                        a, fa = b, fb
            zeros += [(x, 0, m, 1) for x, n, m, p in zeros if (n, p) == (1, 0)]
            _JNJNP.extend(sorted(zeros, key=lambda z: (z[0], z[3])))
    return _JNJNP


def r_jnjnp(i, nt, k):
    """Output i of `bessel_j_jprime_zeros`, for zero k; the count `nt` doesn't change them."""
    zeros = _jnjnp_zeros()
    if len(zeros) < 1200:
        raise ArithmeticError("too few zeros")
    return mp.mpf(zeros[int(k) - 1][i])


def _kelvin_zero(kind, k):
    """Zero k of ber, bei, ker, kei, ber', bei', ker', or kei' (`_kelvin` kind 0 to 7), from the
    sign change within 1 of sqrt(2) pi (k + c / 8), by the leading terms of the expansions of
    DLMF 10.67 (the spacing is sqrt(2) pi)."""
    c = [-3, 1, -5, -1, 3, -1, -3, 1][kind]
    x0 = mp.sqrt(2) * mp.pi * (k + mp.mpf(c) / 8)
    f = lambda x: _kelvin(kind, x)[0]
    a, b = x0 - 1, x0 + 1
    fa, fb = f(a), f(b)
    if (fa < 0) == (fb < 0):
        raise ArithmeticError(f"no zero near {x0}")
    return _solve(f, a, b, fa, fb)


def _erf_zero(k):
    """Zero k of erf in the first quadrant, by Newton's method from DLMF 7.13.1 and 7.13.2; it
    must not move by more than a quarter of the spacing, about pi / (sqrt(2) lambda)."""
    lam = mp.sqrt((k - mp.mpf(1) / 8) * mp.pi)
    mu = mp.log(lam * mp.sqrt(2 * mp.pi))
    t = (1 - mu + mu**2 / 2) / (16 * lam**3)
    z0 = mp.mpc(lam - mu / (4 * lam) + t, lam + mu / (4 * lam) + t)
    z = mp.findroot(
        mp.erf, z0, solver="newton", df=lambda z: 2 / mp.sqrt(mp.pi) * mp.exp(-(z**2))
    )
    if abs(z - z0) > mp.pi / (4 * mp.sqrt(2) * lam):
        raise ArithmeticError(f"zero {z} is not near {z0}")
    return z


def _fresnel_zero(i, k):
    """Zero k of S (i = 0) or C (i = 1) in the first quadrant, by Newton's method from
    DLMF 7.13.5 to 7.13.7; it must not move by more than a quarter of the spacing, about
    2 / lambda."""
    lam = 2 * mp.sqrt(k) if i == 0 else mp.sqrt(4 * k - 1)
    alpha = 2 / mp.pi * mp.log(mp.pi * lam)
    z0 = mp.mpc(
        lam + alpha * (alpha * mp.pi - 4) / (8 * mp.pi * lam**3), alpha / (2 * lam)
    )
    f, df = (mp.fresnels, mp.sin) if i == 0 else (mp.fresnelc, mp.cos)
    z = mp.findroot(f, z0, solver="newton", df=lambda z: df(mp.pi * z**2 / 2))
    if abs(z - z0) > 1 / (2 * lam):
        raise ArithmeticError(f"zero {z} is not near {z0}")
    return z


def _gauss_rule(base, n, k, a, b):
    """Node and weight k (from 1, ascending) of the n-point Gauss rule of the Jacobi (base "j",
    with alpha = a and beta = b), generalized Laguerre ("l", with alpha = a), or Hermite ("h")
    polynomials.

    The node is found by bisection to 30 digits with Sturm sequences: the number of zeros of p_n
    above x is the number of sign changes of p_0(x), ..., p_n(x), with the recurrences of
    DLMF 18.9.2 and Table 18.9.1 (normalized to positive leading coefficients). The weights are
    c / ((1 - x^2) p_n'(x)^2) (Jacobi, with p_n' from DLMF 18.9.15),
    Gamma(n + a + 1) / (n! x p_n'(x)^2) (Laguerre, with 18.9.23), and
    2^(n-1) n! sqrt(pi) / (n p_(n-1)(x))^2 (Hermite)."""
    one, a, b = mp.mpf(1), mp.mpf(a), mp.mpf(b)
    if base == "j":
        lo, hi = -one, one
        p1 = lambda x: a + 1 + (a + b + 2) * (x - 1) / 2
        s, t = a + b, a * a - b * b

        def rec(j):
            u = 2 * (j + 1) * (j + s + 1) * (2 * j + s)
            return (
                (2 * j + s) * (2 * j + s + 1) * (2 * j + s + 2) / u,
                t * (2 * j + s + 1) / u,
                2 * (j + a) * (j + b) * (2 * j + s + 2) / u,
            )

    elif base == "l":
        lo, hi = mp.mpf(0), 4 * n + 2 * a + 10
        p1 = lambda x: x - a - 1
        rec = lambda j: (one, -(2 * j + a + 1), j * (j + a))
    else:
        lo, hi = -mp.sqrt(2 * n + 2), mp.sqrt(2 * n + 2)
        p1 = lambda x: 2 * x
        rec = lambda j: (2 * one, mp.mpf(0), 2 * one * j)
    coef = [rec(j) for j in range(1, n)]

    def above(x):
        """The number of zeros of p_n above x."""
        p, q, count = p1(x), one, 0
        if p < 0:
            count += 1
        for aj, bj, cj in coef:
            p, q = (aj * x + bj) * p - cj * q, p
            if (p < 0) != (q < 0):
                count += 1
        return count

    if (base == "h" or (base == "j" and a == b)) and 2 * k - 1 == n:
        x = mp.mpf(0)  # the middle node of a symmetric rule
    else:
        while hi - lo > mp.mpf(10) ** -30 * max(abs(lo), abs(hi)):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if above(mid) > n - k else (lo, mid)
        x = (lo + hi) / 2
    if base == "j":
        c = 2 ** (a + b + 1) * mp.gamma(n + a + 1) * mp.gamma(n + b + 1)
        c /= mp.gamma(n + a + b + 1) * mp.factorial(n)
        dp = (n + a + b + 1) / 2 * mp.jacobi(n - 1, a + 1, b + 1, x)
        w = c / ((1 - x) * (1 + x) * dp**2)
    elif base == "l":
        w = mp.gamma(n + a + 1) / (
            mp.factorial(n) * x * mp.laguerre(n - 1, a + 1, x) ** 2
        )
    else:
        w = (
            2 ** (n - 1)
            * mp.factorial(n)
            * mp.sqrt(mp.pi)
            / (n * mp.hermite(n - 1, x)) ** 2
        )
    return x, w


def _chebyshev_rule(base, n, k):
    """Node and weight k (from 1, ascending) of the n-point Gauss rule of T_n (base "t") or U_n
    ("u"), in closed form: the zeros of DLMF 18.5.1 and 18.5.2, and the weights of DLMF 3.5.23
    for T_n, and pi / (n + 1) sin(k pi / (n + 1))^2 for U_n."""
    if base == "t":
        th = (2 * k - 1) * mp.pi / (2 * n)
        return -mp.cos(th) if 2 * k - 1 != n else mp.mpf(0), mp.pi / n
    th = k * mp.pi / (n + 1)
    return -mp.cos(th) if 2 * k != n + 1 else mp.mpf(0), mp.pi / (n + 1) * mp.sin(
        th
    ) ** 2


def _gauss(i, n, k, base, a=0, b=0, scale=1, shift=False):
    """Node (i = 0) or weight (i = 1) k of the n-point rule of `_gauss_rule` or `_chebyshev_rule`,
    for x = scale t, or for x = (1 + t) / 2 with shift, which multiply the weights by scale or by
    2^(-a - b - 1) (for the Jacobi weight function)."""
    n, k = int(n), int(k)
    with wdps(xd(n) + 10):
        if base in "tu":
            a = b = mp.mpf(-1 if base == "t" else 1) / 2
            t, w = _chebyshev_rule(base, n, k)
        else:
            t, w = _gauss_rule(base, n, k, a, b)
        if shift:
            return (1 + t) / 2 if i == 0 else w * 2 ** (-a - b - 1)
        return scale * (t if i == 0 else w)


def r_bernoulli(n):
    return mp.bernoulli(int(n))


def r_euler(n):
    return mp.eulernum(int(n))


def r_gauss(f, i, n, k, *p):
    """Node (i = 0) or weight (i = 1) k of the n-point rule of `roots_{f}`, with its parameters
    p, as a transformed Jacobi, Laguerre, or Hermite rule."""
    p = [mp.mpf(x) for x in p]
    match f:
        case "jacobi" | "sh_jacobi":
            return _gauss(i, n, k, "j", *p, shift=f.startswith("sh_"))
        case "gengegenbauer":
            return _gauss(i, n, k, "j", p[0] - p[1], p[1] - 1, shift=True)
        case "legendre" | "sh_legendre":
            return _gauss(i, n, k, "j", shift=f.startswith("sh_"))
        case "gegenbauer":
            return _gauss(i, n, k, "j", p[0] - 0.5, p[0] - 0.5)
        case "genlaguerre" | "laguerre":
            return _gauss(i, n, k, "l", *p)
        case "hermite_h" | "hermite_he":
            return _gauss(i, n, k, "h", scale=1 if f == "hermite_h" else mp.sqrt(2))
    # the Chebyshev rules: C and S are 2 T and U of x / 2
    base, scale = {"t": ("t", 1), "u": ("u", 1), "c": ("t", 2), "s": ("u", 2)}[f[-1]]
    return _gauss(i, n, k, base, scale=scale, shift=f.startswith("sh_"))


# ---------------------------------------------------------------- specs

SPECS: list[Spec] = []
# section titles, in order, with their heading level: 2 for the groups of the crate docs
# (src/lib.rs), and 3 for their subgroups
SECTIONS: dict[str, int] = {}


def section(title: str, level: int = 2) -> None:
    SECTIONS[title] = level


def add(name: str, args: str, ref: FunctionType, *domains: Domain, **kw: Any) -> None:
    """Add a spec for the function `name` with the space-separated argument names `args` to the
    last section."""
    SPECS.append(
        Spec(name, tuple(args.split()), ref, list(domains), list(SECTIONS)[-1], **kw)
    )


def root_domain(
    fwd: str,
    F: Callable[..., Any],
    args: str,
    dists: Mapping[str, Dist],
    where: tuple[str, Callable[..., bool]] | None = None,
    n: int = 500,
) -> Domain:
    """Domain of the inverse of `fwd` with respect to one of its arguments: the arguments of `fwd`
    are sampled from `dists` (in order, and if given, until `where` holds for them), and in
    `args`, p = F(...) at 20 digits, rounded to f64 (and resampled if 0 or 1), replaces the
    argument, so that the root is near its sample."""
    names = list(dists)
    call = f"p = {fwd}({', '.join(names)})"
    cond = f"{where[0]}, " if where else ""
    pretty_cond = cond.replace(">=", "\u2265")
    pretty_call = f"p = {fwd}({', '.join(GREEK.get(a, a) for a in names)})"

    def sampler(r: Random) -> list[float]:
        while True:
            v = {k: d.sample(r) for k, d in dists.items()}
            if where and not where[1](*v.values()):
                continue
            with mp.workdps(20):
                v["p"] = float(F(*map(mp.mpf, v.values())))
            if 0 < v["p"] < 1:
                return [v[a] for a in args.split()]

    dom = D(*dists.values())
    return D(
        sampler=sampler,
        label=f"{dom.desc(names)}, {cond}{call}",
        pretty=f"{dom.pretty(names)}<br>{pretty_cond}{pretty_call}",
        n=n,
    )


INV_E = float(-1 / mp.e)  # the smallest f64 > -1/e
if INV_E < -1 / mp.e:
    INV_E = math.nextafter(INV_E, 0)


def a_pm_sqrt_a(r: Random) -> tuple[float, float]:
    """(a, x) with a in [10, 1e8] (log), and x in a +- 6 sqrt(a): the transition region."""
    a = 10 ** r.uniform(1, 8)
    return a, max(0.0, a + r.uniform(-6, 6) * a**0.5)


A_PM_SQRT_A = D(
    sampler=a_pm_sqrt_a,
    label="a in [10, 1e8] (log), x in a +- 6 sqrt(a)",
    pretty=f"{nobreak('a', IN, LU(10, 1e8).pretty())} (log), "
    + nobreak("x", IN, f"a {PM} 6{SQRT}a"),
)


def ab_mu_pm_sigma(r: Random) -> tuple[float, float, float]:
    """(a, b, x) with a, b in [100, 1e8] (log), and x in mu +- 6 sigma: the transition region."""
    a, b = 10 ** r.uniform(2, 8), 10 ** r.uniform(2, 8)
    mu = a / (a + b)
    return a, b, mu + r.uniform(-6, 6) * (mu * (1 - mu) / (a + b + 1)) ** 0.5


AB_MU_PM_SIGMA = D(
    sampler=ab_mu_pm_sigma,
    label="a, b in [100, 1e8] (log), x in mu +- 6 sigma",
    pretty=f"{nobreak('a, b', IN, LU(100, 1e8).pretty())} (log), "
    + nobreak("x", IN, f"{MU} {PM} 6{SIGMA}"),
)


def v_near_half_int(r: Random) -> tuple[float, float]:
    """(v, x) with v in -k - 1/2 +- [1e-12, 1e-3] (log) for k in {0, ..., 20}, and x in [0, 100]."""
    t = 10 ** r.uniform(-12, -3)
    return -r.randint(0, 20) - 0.5 + r.choice((-1, 1)) * t, r.uniform(0, 100)


V_NEAR_HALF_INT = D(
    sampler=v_near_half_int,
    label="v in -k - 1/2 +- [1e-12, 1e-3] (log), k in {0, ..., 20}, x in [0, 100]",
    pretty=nobreak("v", IN, f"{MINUS}k {MINUS} {HALF} {PM} {LU(1e-12, 1e-3).pretty()}")
    + f" (log), {nobreak('k', IN, Int(0, 20).pretty())}, {nobreak('x', IN, U(0, 100).pretty())}",
)

CZ = D(U(-20, 20), U(-20, 20), n=5000)

# the number of zeros or numbers of the functions with a const `N`; must match src/main.rs and
# harness.cpp
ZEROS_N = 500


def each_k(hi: int = ZEROS_N, lo: int = 1) -> Domain:
    """Each k in {lo, ..., hi}, once."""
    return D(Each(lo, hi), n=hi - lo + 1)


def nk_domain(nd: Dist, params: Mapping[str, Dist], n: int = 1000) -> Domain:
    """n from `nd`, k in {1, ..., n}, and the `params`."""

    def sample(r: Random) -> tuple[float, ...]:
        m = int(nd.sample(r))
        return (m, r.randint(1, m), *(d.sample(r) for d in params.values()))

    log = " (log)" if _is_log(nd) else ""
    return D(
        sampler=sample,
        label=", ".join(
            [f"n in {nd.desc()}", "k in {1, ..., n}"]
            + [f"{a} in {d.desc()}" for a, d in params.items()]
        ),
        pretty=", ".join(
            [
                nobreak("n", IN, nd.pretty()) + log,
                nobreak("k", IN, f"{{1, {ELLIPSIS}, n}}"),
            ]
            + [nobreak(GREEK.get(a, a), IN, d.pretty()) for a, d in params.items()]
        ),
        n=n,
    )


section("Airy functions")
for i, f in enumerate(["Ai", "Ai'", "Bi", "Bi'"]):
    add(
        f"airy[{i}]",
        "x",
        lambda x, i=i: _airy(x)[i],
        D(U(-20, 20)),
        D(LU(1e-300, 1e6, 0)),
        note="airy",
    )
for i in range(4):
    dom = D(U(0, 20)) if i < 2 else D(U(-20, 20))
    dom2 = D(LU(1e-300, 1e6)) if i < 2 else D(LU(1e-300, 1e6, 0))
    add(
        f"airy_scaled[{i}]",
        "x",
        lambda x, i=i: _airy(x, True)[i],
        dom,
        dom2,
        note="airy",
    )
for f, scaled in [("airy", False), ("airy_scaled", True)]:
    for i in range(4):
        add(
            f"{f}[{i}]#c",
            "re im",
            lambda re, im, i=i, scaled=scaled: c_airy(re, im, scaled)[i],
            CZ,
            complex=True,
        )
for i in range(4):
    add(
        f"airy_integrals[{i}]",
        "x",
        lambda x, i=i: _airy_int(x)[i],
        D(U(0, 30)),
        D(LU(1e-300, 1e3)),
        note="zeros_c" if i == 3 else "",
    )

for f in ["ai", "bi"]:
    for i in range(4):
        add(
            f"airy_{f}_zeros.{i}",
            "k",
            lambda k, i=i, bi=f == "bi": _airy_zeros(i, k, bi),
            each_k(),
            note="zeros_n",
        )

section("Elliptic functions and integrals")
for f, r in [("ellipk", r_ellipk), ("ellipe", r_ellipe)]:
    add(f, "m", r, D(U(0, 1)), D(Offset(1, LU(1e-16, 1), -1)), D(LU(1e-300, 1e300, -1)))
for f, r in [("ellipkinc", r_ellipkinc), ("ellipeinc", r_ellipeinc)]:
    add(
        f,
        "phi m",
        r,
        D(U(0, math.pi / 2), U(0, 1)),
        D(LU(1e-3, 1e3, 0), Offset(1, LU(1e-3, 1e5), -1)),
        D(LU(1e-3, 1e300, 0), Offset(1, LU(1e-3, 1e300), -1)),
    )
add("ellipkm1", "p", r_ellipkm1, D(U(0, 1)), D(LU(1e-300, 1)))
add(
    "elliprc",
    "x y",
    r_elliprc,
    D(LU(1e-3, 1e3), LU(1e-3, 1e3)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e300)),
)
for f, r in [("elliprd", r_elliprd), ("elliprf", r_elliprf), ("elliprg", r_elliprg)]:
    add(
        f,
        "x y z",
        r,
        D(U(0, 10), U(0, 10), U(0, 10)),
        D(LU(1e-300, 1e300), LU(1e-300, 1e300), LU(1e-300, 1e300)),
    )
add(
    "elliprj",
    "x y z p",
    r_elliprj,
    D(U(0, 10), U(0, 10), U(0, 10), U(0, 10)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e300), LU(1e-300, 1e300), LU(1e-300, 1e300)),
)
for i in range(4):
    add(
        f"ellipj.{i}",
        "u m",
        lambda u, m, i=i: _ellipj(u, m)[i],
        D(U(-10, 10), U(0, 1)),
        D(LU(1e-3, 1e3, 0), Offset(1, LU(1e-16, 1), -1)),
        note="zeros_c" if i != 2 else "",
    )

section("Ellipsoidal harmonics")
H2, K2 = "h\u00b2", "k\u00b2"
LAME_HK = LU(1e-2, 1e2)
LAME_S = LU(1e-12, 1e3)


def lame_args(r: Random) -> tuple[float, float, int, int]:
    """(h2, k2, n, p) with h2 and k2 - h2 sampled from `LAME_HK`, n in {0, ..., 10}, and p in
    {1, ..., 2n + 1}."""
    h2 = LAME_HK.sample(r)
    n = r.randint(0, 10)
    return h2, h2 + LAME_HK.sample(r), n, r.randint(1, 2 * n + 1)


def lame_s(r: Random, k2: float) -> float:
    """s = k (1 + t), with t sampled from `LAME_S`."""
    return math.sqrt(k2) * (1 + LAME_S.sample(r))


def lame_signs(r: Random) -> tuple[int, int]:
    return r.choice((-1, 1)), r.choice((-1, 1))


LAME = (
    f"{nobreak(H2, IN, LAME_HK.pretty())} (log), "
    + f"{nobreak(f'{K2} {MINUS} {H2}', IN, LAME_HK.pretty())} (log), "
    + f"{nobreak('n', IN, Int(0, 10).pretty())}, "
    + nobreak("p", IN, f"{{1, {ELLIPSIS}, 2n + 1}}")
)
LAME_LABEL = (
    "h2, k2 - h2 in [0.01, 100] (log), n in {0, ..., 10}, p in {1, ..., 2n + 1}"
)
LAME_S_PRETTY = f"{nobreak('s', IN, f'k (1 + {LAME_S.pretty()})')} (log)"
LAME_SIGNS = nobreak("signm, signn", IN, f"{{{PM}1}}")
add(
    "ellip_harm",
    "h2 k2 n p s signm signn",
    r_ellip_harm,
    D(
        sampler=lambda r: (
            *(a := lame_args(r)),
            r.uniform(0, 2 * math.sqrt(a[1])),
            *lame_signs(r),
        ),
        label=f"{LAME_LABEL}, s in [0, 2k], signm, signn in +-1",
        pretty=f"{LAME}, {nobreak('s', IN, '[0, 2k]')}, {LAME_SIGNS}",
    ),
    D(
        sampler=lambda r: (
            *(a := lame_args(r)),
            lame_s(r, a[1]),
            *lame_signs(r),
        ),
        label=f"{LAME_LABEL}, s in k (1 + [1e-12, 1e3]) (log), signm, signn in +-1",
        pretty=f"{LAME}, {LAME_S_PRETTY}, {LAME_SIGNS}",
    ),
    cpp=False,
    note="zeros_c lame",
)
add(
    "ellip_harm_2",
    "h2 k2 n p s",
    r_ellip_harm_2,
    D(
        sampler=lambda r: (*(a := lame_args(r)), lame_s(r, a[1])),
        label=f"{LAME_LABEL}, s in k (1 + [1e-12, 1e3]) (log)",
        pretty=f"{LAME}, {LAME_S_PRETTY}",
    ),
    cpp=False,
)
add(
    "ellip_normal",
    "h2 k2 n p",
    r_ellip_normal,
    D(sampler=lame_args, label=LAME_LABEL, pretty=LAME, n=2000),
    cpp=False,
)

section("Bessel functions")
add(
    "bessel_j",
    "v x",
    r_jv,
    D(U(0, 20), U(0, 50)),
    D(U(-100, 100), LU(1e-300, 1e4)),
    note="osc",
)
add(
    "bessel_y",
    "v x",
    r_yv,
    D(U(0, 20), U(0, 50)),
    D(U(-100, 100), LU(1e-300, 1e4)),
    note="osc",
)
add("bessel_i", "v x", r_iv, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 700)))
add("bessel_k", "v x", r_kv, D(U(0, 20), U(0, 50)), D(U(0, 100), LU(1e-300, 700)))
add(
    "bessel_je",
    "v x",
    r_jv,
    D(U(0, 20), U(0, 50)),
    D(U(-100, 100), LU(1e-300, 1e8)),
    note="osc",
)
add(
    "bessel_ye",
    "v x",
    r_yv,
    D(U(0, 20), U(0, 50)),
    D(U(-100, 100), LU(1e-300, 1e8)),
    note="osc",
)
add("bessel_ie", "v x", r_ie, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 1e8)))
add("bessel_ke", "v x", r_ke, D(U(0, 20), U(0, 50)), D(U(0, 100), LU(1e-300, 1e8)))
for f, r in [
    ("bessel_j", c_besselj),
    ("bessel_y", c_bessely),
    ("bessel_i", c_besseli),
    ("bessel_k", c_besselk),
    ("hankel_1", c_hankel1),
    ("hankel_2", c_hankel2),
]:
    add(
        f"{f}#c",
        "v re im",
        r,
        D(U(-20, 20), U(-50, 50), U(-50, 50), n=5000),
        complex=True,
        note="osc_c" if f in ("bessel_j", "bessel_y") else "",
    )
for f, r in [("hankel_1e", c_hankel1e), ("hankel_2e", c_hankel2e)]:
    add(f, "v re im", r, D(U(-20, 20), U(-50, 50), U(-50, 50), n=5000), complex=True)
for f, r in [
    ("bessel_je", c_besselje),
    ("bessel_ye", c_besselye),
    ("bessel_ie", c_besselie),
    ("bessel_ke", c_besselke),
]:
    add(
        f"{f}#c",
        "v re im",
        r,
        D(U(-20, 20), U(-50, 50), U(-50, 50), n=5000),
        complex=True,
        note="osc_cz" if f in ("bessel_je", "bessel_ye") else "",
    )
add(
    "wright_bessel",
    "a b x",
    r_wright_bessel,
    D(U(0, 5), U(0, 5), U(0, 100), n=2000),
    timeout=60,
)
add(
    "log_wright_bessel",
    "a b x",
    r_log_wright_bessel,
    D(U(0, 5), U(0, 5), LU(1e-3, 1e3), n=2000),
    timeout=60,
    note="zeros_c",
)
add(
    "iv_ratio",
    "v x",
    r_iv_ratio,
    D(U(0.5, 20), U(0, 50)),
    D(LU(0.5, 1e4), LU(1e-300, 1e4)),
)
add(
    "iv_ratio_c",
    "v x",
    r_iv_ratio_c,
    D(U(0.5, 20), U(0, 50)),
    D(LU(0.5, 1e4), LU(1e-300, 1e4)),
)
add("iv_ratioinv", "v r", r_iv_ratioinv, D(U(0.5, 20), U(0, 1), n=2000), note="inv")
for i, r in enumerate([r_lmbda, r_lmbdap]):
    add(
        f"jahnke_emden_lambda.{i}",
        "v x",
        r,
        D(U(0, 20), U(0, 50)),
        D(U(0, 100), LU(1e-300, 1e4)),
        D(Int(0, 100), LU(1e-300, 1e4)),
        cpp=False,
        note="lmbda",
    )

section("Zeros of Bessel functions", 3)
for i in range(4):
    add(
        f"bessel_zeros.{i}",
        "v k",
        lambda v, k, i=i: _bessel_zero(i, v, k),
        D(Int(0, 100), Int(1, ZEROS_N), n=1000),
        note="zeros_n",
    )
for i in range(4):
    add(
        f"bessel_j_jprime_zeros.{i}",
        "N k",
        lambda nt, k, i=i: r_jnjnp(i, nt, k),
        D(
            sampler=lambda r: (nt := 100 * r.randint(1, 12), r.randint(1, nt)),
            label="N in {100, 200, ..., 1200}, k in {1, ..., N}",
            pretty=f"{nobreak('N', IN, f'{{100, 200, {ELLIPSIS}, 1200}}')}, "
            + nobreak("k", IN, f"{{1, {ELLIPSIS}, N}}"),
            n=1000,
        ),
        timeout=120,  # the first call per process finds all zeros
        note="jnjnp",
    )
for kind, f in enumerate(["y0", "y1", "y1_prime"]):
    for i in range(2):
        add(
            f"bessel_{f}_zeros.{i}",
            "c k",
            lambda c, k, i=i, kind=kind: _y_zero(i, kind, c, k),
            *(
                D(
                    Int(c, c),
                    Each(1, ZEROS_N),
                    n=ZEROS_N,
                    label=f"{label}, k in each of {{1, ..., {ZEROS_N}}}",
                    pretty=f"{label}, {nobreak('k', IN, Int(1, ZEROS_N).pretty())}",
                )
                for c, label in enumerate(["real", "complex"])
            ),
            complex=True,
            note="zeros_n y_zeros",
        )

section("Faster versions of common Bessel functions", 3)
for f, r in [("j0", r_j0), ("j1", r_j1), ("y0", r_y0), ("y1", r_y1)]:
    add(f"bessel_{f}", "x", r, D(U(0, 30)), D(LU(1e-300, 1e8)), note="osc")
for f, r in [("i0", r_i0), ("i1", r_i1)]:
    add(f"bessel_{f}", "x", r, D(U(-30, 30)), D(LU(1e-300, 700)))
for f, r in [("k0", r_k0), ("k1", r_k1)]:
    add(f"bessel_{f}", "x", r, D(U(0, 30)), D(LU(1e-300, 700)))
for f, r in [("i0e", r_i0e), ("i1e", r_i1e)]:
    add(f"bessel_{f}", "x", r, D(U(-30, 30)), D(LU(1e-300, 1e300, 0)))
for f, r in [("k0e", r_k0e), ("k1e", r_k1e)]:
    add(f"bessel_{f}", "x", r, D(U(0, 30)), D(LU(1e-300, 1e300)))

section("Integrals of Bessel functions", 3)
for f, r, hi in [
    ("it1j0y0.0", r_itj0, 1e4),
    ("it1j0y0.1", r_ity0, 1e4),
    ("it2j0y0.0", r_it2j0, 1e4),
    ("it2j0y0.1", r_it2y0, 1e3),
    ("it1i0k0.0", r_iti0, 700),
    ("it1i0k0.1", r_itk0, 700),
    ("it2i0k0.0", r_it2i0, 700),
    ("it2i0k0.1", r_it2k0, 700),
]:
    add(
        f,
        "x",
        r,
        D(U(0, 30)),
        D(LU(1e-300, hi)),
        note="zeros_c" if f[-1] == "1" else "",
    )
add(
    "besselpoly",
    "a lambda v",
    r_besselpoly,
    D(U(0, 10), U(0, 10), U(0, 10)),
    note="zeros_c",
)

section("Derivatives of Bessel functions", 3)
for f, r in [("j", r_jvp), ("y", r_yvp), ("i", r_ivp), ("k", r_kvp)]:
    add(
        f"bessel_{f}_prime",
        "v x n",
        r,
        D(U(0, 20), U(0, 50), Int(1, 3)),
        cpp=False,
        note="osc" if f in "jy" else "",
    )
for f, r in [("j", c_jvp), ("y", c_yvp), ("i", c_ivp), ("k", c_kvp)]:
    add(
        f"bessel_{f}_prime#c",
        "v re im n",
        r,
        D(U(-20, 20), U(-50, 50), U(-50, 50), Int(1, 3), n=5000),
        complex=True,
        cpp=False,
        note="osc_cz" if f in "jy" else "",
    )
for f, r in [("hankel_1_prime", c_hankel1p), ("hankel_2_prime", c_hankel2p)]:
    add(
        f,
        "v re im n",
        r,
        D(U(-20, 20), U(-50, 50), U(-50, 50), Int(1, 3), n=5000),
        complex=True,
        cpp=False,
    )

section("Spherical Bessel functions", 3)
for f, r, rp in [("j", r_sph_j, r_sph_jp), ("y", r_sph_y, r_sph_yp)]:
    for g, ref in [(f"sph_bessel_{f}", r), (f"sph_bessel_{f}_prime", rp)]:
        add(
            g,
            "n x",
            ref,
            D(Int(0, 50), U(0, 50)),
            D(Int(0, 100), LU(1e-300, 1e4)),
            note="osc_sph",
        )
for f, r, rp in [("i", r_sph_i, r_sph_ip), ("k", r_sph_k, r_sph_kp)]:
    for g, ref in [(f"sph_bessel_{f}", r), (f"sph_bessel_{f}_prime", rp)]:
        add(g, "n x", ref, D(Int(0, 50), U(0, 50)), D(Int(0, 100), LU(1e-300, 700)))

section("Riccati-Bessel functions", 3)
for f, r, rp in [("j", r_riccati_j, r_riccati_jp), ("y", r_riccati_y, r_riccati_yp)]:
    for i, ref in enumerate([r, rp]):
        add(
            f"riccati_{f}.{i}",
            "n x",
            ref,
            D(Int(0, 100), U(0, 50)),
            D(Int(0, 100), LU(1e-300, 1e4)),
            note="riccati",
        )

section("Struve functions")
add(
    "struve_h",
    "v x",
    r_struve_h,
    D(U(0, 10), U(0, 50)),
    D(U(-50, 50), LU(1e-300, 1e3)),
    note="osc",
    timeout=60,
)
add("struve_l", "v x", r_struve_l, V_NEAR_HALF_INT, timeout=60)
add("itstruve0", "x", r_itstruve0, D(U(0, 30)), D(LU(1e-300, 1e4)))
add("it2struve0", "x", r_it2struve0, D(U(0, 30)), D(LU(1e-300, 1e4)))
add("itmodstruve0", "x", r_itmodstruve0, D(U(0, 30)), D(LU(1e-300, 700)))

section("Raw statistical functions")

section("Binomial distribution", 3)


def kn_sampler(r: Random, n_hi: int = 1000) -> tuple[int, int]:
    n = int(10 ** r.uniform(0, math.log10(n_hi)))
    return r.randint(0, n), n


def kn_below_n(r: Random) -> tuple[int, int]:
    """`kn_sampler`, with k < n."""
    k, n = kn_sampler(r)
    return min(k, n - 1), n


def k_near_np(r: Random) -> tuple[int, int, float]:
    """(k, n, p) with n in [1e3, 2e9] (log), p in [1e-6, 1] (log), and k in np +- 6 sigma."""
    n, p = int(10 ** r.uniform(3, math.log10(2e9))), 10 ** r.uniform(-6, 0)
    k = round(n * p + r.uniform(-6, 6) * (n * p * (1 - p)) ** 0.5)
    return min(max(k, 0), n), n, p


KN = f"{nobreak('n', IN, '{1, ' + ELLIPSIS + ', 1000}')} (log), {nobreak('k', IN, '{0, ' + ELLIPSIS + ', n}')}"
for f, r in [("bdtr", r_bdtr), ("bdtrc", r_bdtrc)]:
    add(
        f,
        "k n p",
        r,
        D(
            sampler=lambda r: (*kn_sampler(r), r.uniform(0, 1)),
            label="n in {1, ..., 1000} (log), k in {0, ..., n}, p in [0, 1]",
            pretty=f"{KN}, {nobreak('p', IN, U(0, 1).pretty())}",
            n=5000,
        ),
        D(
            sampler=k_near_np,
            label="n in [1e3, 2e9] (log), p in [1e-6, 1] (log), k in np +- 6 sigma",
            pretty=f"{nobreak('n', IN, LU(1e3, 2e9).pretty())}, "
            + f"{nobreak('p', IN, LU(1e-6, 1).pretty())} (log), "
            + nobreak("k", IN, f"np {PM} 6{SIGMA}"),
            n=2000,
        ),
        timeout=60,
        note="binom_dist",
    )
add(
    "bdtri",
    "k n y",
    r_bdtri,
    D(
        sampler=lambda r: (*kn_below_n(r), r.uniform(0, 1)),
        label="n in {1, ..., 1000} (log), k in {0, ..., n - 1}, y in [0, 1]",
        pretty=f"{nobreak('n', IN, '{1, ' + ELLIPSIS + ', 1000}')} (log), "
        + f"{nobreak('k', IN, '{0, ' + ELLIPSIS + f', n {MINUS} 1}}')}, "
        + nobreak("y", IN, U(0, 1).pretty()),
        n=2000,
    ),
    timeout=120,
    note="inv",
)

add(
    "bdtrik",
    "y n p",
    r_bdtrik,
    D(U(0, 1), LU(1, 1000), U(0, 1), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)
add(
    "bdtrin",
    "k y p",
    r_bdtrin,
    D(LU(1e-2, 1000), U(0, 1), U(0, 1), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)

section("Beta distribution", 3)
add(
    "btdtria",
    "p b x",
    r_btdtria,
    D(U(0, 1), LU(1e-2, 100), U(0, 1), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)
add(
    "btdtrib",
    "a p x",
    r_btdtrib,
    D(LU(1e-2, 100), U(0, 1), U(0, 1), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)

section("F distribution", 3)
for f, r in [("fdtr", r_fdtr), ("fdtrc", r_fdtrc)]:
    add(f, "a b x", r, D(LU(1e-3, 1e6), LU(1e-3, 1e6), LU(1e-30, 1e30)), timeout=60)
add(
    "fdtri",
    "a b y",
    r_fdtri,
    D(LU(1e-3, 1e3), LU(1e-3, 1e3), U(0, 1), n=2000),
    timeout=120,
    note="inv",
)
add(
    "fdtridfd",
    "a p x",
    r_fdtridfd,
    # for x < 1, fdtr can have a maximum in b, and two roots
    root_domain(
        "fdtr",
        r_fdtr,
        "a p x",
        {"a": LU(1e-2, 1e3), "b": LU(1e-2, 1e3), "x": LU(1, 1e3)},
    ),
    cpp=False,
    timeout=120,
    note="inv root",
)

section("Noncentral F distribution", 3)
add(
    "ncfdtr",
    "a b nc x",
    r_ncfdtr,
    D(LU(1e-2, 1e3), LU(1e-2, 1e3), LU(1e-3, 1e3), LU(1e-3, 1e3), n=2000),
    cpp=False,
    timeout=120,
)
add(
    "ncfdtri",
    "a b nc p",
    r_ncfdtri,
    D(LU(1e-2, 1e3), LU(1e-2, 1e3), LU(1e-3, 1e3), U(0, 1), n=500),
    cpp=False,
    timeout=120,
    note="inv",
)
NCF = {"a": LU(1e-2, 1e3), "b": LU(1e-2, 1e3), "nc": LU(1e-3, 100), "x": LU(1e-2, 100)}
add(
    "ncfdtridfd",
    "a p nc x",
    r_ncfdtridfd,
    # for x < 1 + nc/a, ncfdtr can have a maximum in b, and two roots
    root_domain(
        "ncfdtr",
        _ncfdtr_series,
        "a p nc x",
        NCF | {"x": LU(1, 1e3)},
        where=("x >= 1 + nc/a", lambda a, b, nc, x: x >= 1 + nc / a),
    ),
    cpp=False,
    timeout=300,
    note="inv root",
)
add(
    "ncfdtridfn",
    "p b nc x",
    r_ncfdtridfn,
    root_domain("ncfdtr", _ncfdtr_series, "p b nc x", NCF),
    cpp=False,
    timeout=300,
    note="inv root",
)
add(
    "ncfdtrinc",
    "a b p x",
    r_ncfdtrinc,
    root_domain("ncfdtr", _ncfdtr_series, "a b p x", NCF),
    cpp=False,
    timeout=120,
    note="inv",
)

section("Gamma distribution", 3)


def gamma_trans(r: Random) -> tuple[float, float, float]:
    """(a, b, x) with a in [1e-5, 1e5], b in [10, 1e8] (log), and ax in b +- 6 sqrt(b)."""
    a, b = 10 ** r.uniform(-5, 5), 10 ** r.uniform(1, 8)
    return a, b, max(0.0, b + r.uniform(-6, 6) * b**0.5) / a


for f, r in [("gdtr", r_gdtr), ("gdtrc", r_gdtrc)]:
    add(
        f,
        "a b x",
        r,
        D(LU(1e-3, 1e3), LU(1e-3, 1e3), LU(1e-3, 1e3), n=5000),
        D(
            sampler=gamma_trans,
            label="a in [1e-5, 1e5] (log), b in [10, 1e8] (log), ax in b +- 6 sqrt(b)",
            pretty=f"{nobreak('a', IN, LU(1e-5, 1e5).pretty())} (log), "
            + f"{nobreak('b', IN, LU(10, 1e8).pretty())} (log), {nobreak('ax', IN, f'b {PM} 6{SQRT}b')}",
            n=2000,
        ),
        timeout=60,
    )
add(
    "gdtria",
    "p b x",
    r_gdtria,
    D(U(0, 1), LU(1e-2, 1e2), LU(1e-2, 1e2), n=2000),
    timeout=120,
    note="inv",
)
add(
    "gdtrib",
    "a p x",
    r_gdtrib,
    D(LU(1e-2, 1e2), U(0, 1), LU(1e-2, 1e2), n=2000),
    timeout=120,
    note="inv",
)
add(
    "gdtrix",
    "a b p",
    r_gdtrix,
    D(LU(1e-2, 1e2), LU(1e-2, 1e2), U(0, 1), n=2000),
    timeout=120,
    note="inv",
)

section("Generalized hyperbolic distribution", 3)


def gh_sampler(r: Random) -> tuple[float, float, float, float]:
    a = 10 ** r.uniform(-2, 2)
    return r.uniform(-10, 10), r.uniform(-5, 5), a, a * r.uniform(-0.99, 0.99)


for f, r in [
    ("genhyperbolic_pdf", r_genhyperbolic_pdf),
    ("genhyperbolic_logpdf", r_genhyperbolic_logpdf),
]:
    add(
        f,
        "x p a b",
        r,
        D(
            sampler=gh_sampler,
            label="x in [-10, 10], p in [-5, 5], a in [1e-2, 100] (log), b in a [-0.99, 0.99]",
            pretty=f"{nobreak('x', IN, U(-10, 10).pretty())}, {nobreak('p', IN, U(-5, 5).pretty())}, "
            + f"{nobreak('a', IN, LU(1e-2, 100).pretty())} (log), {nobreak('b', IN, 'a [' + MINUS + '0.99, 0.99]')}",
        ),
        note="zeros_c" if "log" in f else "",
    )

section("Generalized inverse Gaussian distribution", 3)
for f, r in [
    ("geninvgauss_pdf", r_geninvgauss_pdf),
    ("geninvgauss_logpdf", r_geninvgauss_logpdf),
]:
    add(
        f,
        "x p b",
        r,
        D(LU(1e-3, 1e3), U(-10, 10), LU(1e-3, 1e3)),
        note="zeros_c" if "log" in f else "",
    )

section("Negative binomial distribution", 3)
for f, r in [("nbdtr", r_nbdtr), ("nbdtrc", r_nbdtrc)]:
    add(
        f,
        "k n p",
        r,
        D(Int(0, 1000), Int(1, 1000), U(0, 1), n=5000),
        D(ILU(1, 2 * 10**9), ILU(1, 2 * 10**9), LU(1e-300, 1), n=2000),
        timeout=60,
    )
add(
    "nbdtri",
    "k n y",
    r_nbdtri,
    D(Int(0, 1000), Int(1, 1000), U(0, 1), n=2000),
    timeout=120,
    note="inv",
)
add(
    "nbdtrik",
    "y n p",
    r_nbdtrik,
    D(U(0, 1), LU(1e-2, 1000), U(0, 1), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)
add(
    "nbdtrin",
    "k y p",
    r_nbdtrin,
    D(LU(1e-2, 1000), U(0, 1), U(0, 1), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)

section("Normal distribution", 3)
add("ndtr", "x", r_ndtr, D(U(-10, 10)), D(LU(1e-300, 38, 0)))
add("ndtri", "p", r_ndtri, D(U(0, 1)), D(LU(1e-300, 0.5)))
add("log_ndtr", "x", r_log_ndtr, D(U(-20, 6)), D(LU(1e-300, 1e10, 0)))
for f, r in [("ndtr", c_ndtr), ("log_ndtr", c_log_ndtr)]:
    add(f"{f}#c", "re im", _c(r), CZ, complex=True)
add("ndtri_exp", "y", r_ndtri_exp, D(U(-20, 0)), D(LU(1e-300, 1e300, -1)))
add(
    "nrdtrimn",
    "p std x",
    r_nrdtrimn,
    D(U(0, 1), LU(1e-3, 1e3), U(-10, 10)),
    note="zeros_c",
)
add("nrdtrisd", "mean p x", r_nrdtrisd, D(U(-10, 10), U(0, 1), U(-10, 10)))

section("Poisson distribution", 3)
for f, r in [("pdtr", r_pdtr), ("pdtrc", r_pdtrc)]:
    add(
        f,
        "k m",
        r,
        D(Int(0, 100), U(0, 100), n=5000),
        D(
            sampler=lambda r: (
                max(0, round((m := 10 ** r.uniform(1, 9)) + r.uniform(-6, 6) * m**0.5)),
                m,
            ),
            label="m in [10, 1e9] (log), k in m +- 6 sqrt(m)",
            pretty=f"{nobreak('m', IN, LU(10, 1e9).pretty())} (log), {nobreak('k', IN, f'm {PM} 6{SQRT}m')}",
            n=2000,
        ),
        timeout=60,
    )
add("pdtri", "k y", r_pdtri, D(Int(0, 100), U(0, 1), n=2000), timeout=120, note="inv")
add(
    "pdtrik",
    "p m",
    r_pdtrik,
    D(U(0, 1), LU(1e-2, 1e4), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)

section("Student's t distribution", 3)
add(
    "stdtr",
    "nu t",
    r_stdtr,
    D(LU(1e-3, 1e3), U(-10, 10)),
    D(LU(1e-3, 1e12), LU(1e-300, 1e300, 0)),
    cpp=False,
    timeout=60,
)
add(
    "stdtri",
    "nu p",
    r_stdtri,
    D(LU(1e-3, 1e3), U(0, 1), n=2000),
    D(LU(1e3, 1e20), U(0, 1), n=2000),
    D(LU(1, 1e20), LU(1e-323, 1e-10), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)
add(
    "stdtridf",
    "p t",
    r_stdtridf,
    root_domain(
        "stdtr", r_stdtr, "p t", {"nu": LU(1e-2, 1e3), "t": LU(1e-2, 10, 0)}, n=2000
    ),
    cpp=False,
    timeout=120,
    note="inv root",
)

section("Noncentral t distribution", 3)
add(
    "nctdtr",
    "nu nc t",
    r_nctdtr,
    D(LU(1e-2, 1e3), U(-10, 10), U(-20, 20), n=1000),
    # `t < 0 < nc`, where the reference cancels
    D(LU(1e-2, 1e3), U(0, 10), U(-20, 0), n=500),
    cpp=False,
    timeout=120,
)
add(
    "nctdtrit",
    "nu nc p",
    r_nctdtrit,
    D(LU(1e-2, 1e3), U(-10, 10), U(0, 1), n=500),
    cpp=False,
    timeout=120,
    note="inv",
)
add(
    "nctdtridf",
    "p nc t",
    r_nctdtridf,
    root_domain(
        "nctdtr",
        r_nctdtr,
        "p nc t",
        {"nu": LU(1e-2, 1e3), "nc": U(-10, 10), "t": U(-20, 20)},
    ),
    cpp=False,
    timeout=300,
    note="inv root",
)
add(
    "nctdtrinc",
    "nu p t",
    r_nctdtrinc,
    D(LU(1e-2, 1e3), U(0, 1), U(-20, 20), n=500),
    cpp=False,
    timeout=300,
    note="inv",
)

section("Chi square distribution", 3)
for f, r in [("chdtr", r_chdtr), ("chdtrc", r_chdtrc)]:
    add(
        f,
        "v x",
        r,
        D(LU(1e-3, 1e3), LU(1e-3, 1e3), n=5000),
        D(
            sampler=lambda r: (
                v := 10 ** r.uniform(1, 8),
                max(0.0, v + r.uniform(-6, 6) * (2 * v) ** 0.5),
            ),
            label="v in [10, 1e8] (log), x in v +- 6 sqrt(2v)",
            pretty=f"{nobreak('v', IN, LU(10, 1e8).pretty())} (log), {nobreak('x', IN, f'v {PM} 6{SQRT}(2v)')}",
            n=2000,
        ),
        timeout=60,
    )
add(
    "chdtri",
    "v y",
    r_chdtri,
    D(LU(1e-3, 1e3), U(0, 1), n=2000),
    timeout=120,
    note="inv",
)
add(
    "chdtriv",
    "p x",
    r_chdtriv,
    D(U(0, 1), LU(1e-3, 1e3), n=2000),
    cpp=False,
    timeout=120,
    note="inv",
)

section("Noncentral chi square distribution", 3)


def chndtr_mean_pm_sigma(r: Random) -> tuple[float, float, float]:
    """(x, v, nc) with v in [1e-3, 1e4], nc in [1, 1e4] (log), and x > 0 in mu +- 10 sigma: the
    transition region and its tails."""
    while True:
        v, nc = 10 ** r.uniform(-3, 4), 10 ** r.uniform(0, 4)
        mu, sigma = v + nc, (2 * (v + 2 * nc)) ** 0.5
        x = mu + r.uniform(-10, 10) * sigma
        if x > 0:
            return x, v, nc


add(
    "chndtr",
    "x v nc",
    r_chndtr,
    D(LU(1e-3, 1e3), LU(1e-3, 1e3), LU(1e-3, 1e3), n=2000),
    D(
        sampler=chndtr_mean_pm_sigma,
        label="v in [1e-3, 1e4], nc in [1, 1e4] (log), x in mu +- 10 sigma (x > 0)",
        pretty=f"{nobreak('v', IN, LU(1e-3, 1e4).pretty())}, "
        + f"{nobreak('nc', IN, LU(1, 1e4).pretty())} (log), "
        + nobreak("x", IN, f"{MU} {PM} 10{SIGMA}")
        + " (x > 0)",
        n=1000,
    ),
    cpp=False,
    timeout=120,
    note="ncx2",
)


def chndtridf_args(r: Random) -> tuple[float, float, float]:
    """(x, p, nc) with x in [1e-3, 1e3], nc in [1e-3, 30] (log), and p in [0, exp(-nc / 2)],
    below the CDF at df = 0, so that there is a root."""
    x, nc = 10 ** r.uniform(-3, 3), 10 ** r.uniform(-3, math.log10(30))
    return x, r.uniform(0, 1) * math.exp(-nc / 2), nc


def chndtrinc_args(r: Random) -> tuple[float, float, float]:
    """(x, v, p) with x, v in [1e-3, 1e3] (log), and p in (0, chdtr(v, x)], below the CDF at
    nc = 0, so that there is a finite root."""
    while True:
        x, v = 10 ** r.uniform(-3, 3), 10 ** r.uniform(-3, 3)
        p = r.uniform(0, 1) * float(mp.gammainc(v / 2, 0, x / 2, regularized=True))
        if p > 0:
            return x, v, p


add(
    "chndtri",
    "p v nc",
    r_chndtri,
    D(U(0, 1), LU(1e-3, 1e3), LU(1e-3, 1e3), n=500),
    cpp=False,
    timeout=120,
    note="inv",
)
add(
    "chndtridf",
    "x p nc",
    r_chndtridf,
    D(
        sampler=chndtridf_args,
        label="x in [1e-3, 1e3], nc in [1e-3, 30] (log), p in [0, exp(-nc/2)]",
        pretty=f"{nobreak('x', IN, LU(1e-3, 1e3).pretty())}, "
        + f"{nobreak('nc', IN, LU(1e-3, 30).pretty())} (log), "
        + nobreak("p", IN, "[0, e^(\u2212nc/2)]"),
        n=500,
    ),
    cpp=False,
    timeout=120,
    note="inv",
)
add(
    "chndtrinc",
    "x v p",
    r_chndtrinc,
    D(
        sampler=chndtrinc_args,
        label="x, v in [1e-3, 1e3] (log), p in (0, chdtr(v, x)]",
        pretty=f"{nobreak('x, v', IN, LU(1e-3, 1e3).pretty())} (log), "
        + nobreak("p", IN, "(0, chdtr(v, x)]"),
        n=500,
    ),
    cpp=False,
    timeout=120,
    note="inv",
)

section("Kolmogorov distribution", 3)
for f, r in [
    ("kolmogorov", r_kolmogorov),
    ("kolmogc", r_kolmogc),
    ("kolmogp", r_kolmogp),
]:
    add(f, "x", r, D(U(0, 3)), D(LU(1e-3, 30)))
add("kolmogi", "q", r_kolmogi, D(U(0, 1), n=2000), D(LU(1e-300, 1), n=2000), note="inv")
add(
    "kolmogci",
    "p",
    r_kolmogci,
    D(U(0, 1), n=2000),
    D(LU(1e-300, 1), n=2000),
    note="inv",
)

section("Kolmogorov-Smirnov distribution", 3)
for f, r in [
    ("smirnov", r_smirnov),
    ("smirnovc", r_smirnovc),
    ("smirnovp", r_smirnovp),
]:
    add(
        f,
        "n x",
        r,
        D(ILU(1, 1000), U(0, 1), n=5000),
        D(ILU(1, 1000), LU(1e-300, 1), n=2000),
        timeout=60,
    )
for f, r in [("smirnovi", r_smirnovi), ("smirnovci", r_smirnovci)]:
    add(
        f,
        "n q" if f == "smirnovi" else "n p",
        r,
        D(ILU(1, 1000), U(0, 1), n=1000),
        timeout=120,
        note="inv",
    )

section("Box-Cox transformation", 3)
add(
    "boxcox",
    "x lambda",
    r_boxcox,
    D(LU(1e-3, 1e3), U(-5, 5)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e3, 0)),
)
add(
    "boxcox1p",
    "x lambda",
    r_boxcox1p,
    D(U(-1, 10), U(-5, 5)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e3, 0)),
)
for f, r in [("inv_boxcox", r_inv_boxcox), ("inv_boxcox1p", r_inv_boxcox1p)]:
    add(
        f,
        "y lambda",
        r,
        D(U(-1, 1), U(-1, 1)),
        D(LU(1e-300, 10, 0), LU(1e-300, 0.1, 0)),
    )

section("Sigmoidal functions", 3)
add("logit", "p", r_logit, D(U(0, 1)), D(LU(1e-300, 0.5)), note="zeros_c")
add("expit", "x", r_expit, D(U(-10, 10)), D(LU(1e-300, 750, 0)))
add("log_expit", "x", r_log_expit, D(U(-20, 20)), D(LU(1e-300, 1e300, 0)))

section("Cram\u00e9r-von Mises", 3)
add("cdf_cvm_inf", "x", r_cdf_cvm_inf, D(U(0, 2)), D(LU(1e-3, 100)))

section("Multivariate normal distribution", 3)
add(
    "bivariate_normal_cdf",
    "h k r",
    r_bivariate_normal_cdf,
    D(U(-10, 10), U(-10, 10), U(-1, 1), n=2000),
    timeout=60,
)


def tvn_sampler(r: Random) -> tuple[float, ...]:
    """(h1, h2, h3, r12, r13, r23) with h in [-3, 3], and a positive definite correlation matrix:
    r12, r13 in [-0.9, 0.9], and r23 in r12 r13 +- 0.9 sqrt((1 - r12^2) (1 - r13^2))."""
    r12, r13 = r.uniform(-0.9, 0.9), r.uniform(-0.9, 0.9)
    r23 = r12 * r13 + r.uniform(-0.9, 0.9) * ((1 - r12 * r12) * (1 - r13 * r13)) ** 0.5
    return r.uniform(-3, 3), r.uniform(-3, 3), r.uniform(-3, 3), r12, r13, r23


add(
    "trivariate_normal_cdf",
    "h1 h2 h3 r12 r13 r23",
    r_trivariate_normal_cdf,
    D(
        sampler=tvn_sampler,
        label="h in [-3, 3], positive definite R",
        pretty=f"{nobreak('h\u2081, h\u2082, h\u2083', IN, U(-3, 3).pretty())}, "
        + f"{nobreak('r\u2081\u2082, r\u2081\u2083', IN, U(-0.9, 0.9).pretty())}, "
        + nobreak("r\u2082\u2083", IN, "r\u2081\u2082r\u2081\u2083")
        # on two lines, to keep the column narrow
        + "<br>"
        + nobreak(
            f"{PM} 0.9{SQRT}((1 {MINUS} r\u2081\u2082\u00b2)(1 {MINUS} r\u2081\u2083\u00b2))"
        ),
        n=1000,
    ),
    timeout=120,
    note="tvn abs",
)

section("Miscellaneous", 3)
add(
    "tukeylambdacdf",
    "x lambda",
    r_tukeylambdacdf,
    D(U(-5, 5), U(-2, 2), n=2000),
    timeout=60,
)
add(
    "von_mises_cdf",
    "k x",
    r_von_mises_cdf,
    D(U(0, 10), U(-math.pi, math.pi), n=2000),
    D(LU(1e-3, 1e5), U(-math.pi, math.pi), n=2000),
    timeout=60,
)
add("cosine_cdf", "x", r_cosine_cdf, D(U(-math.pi, math.pi)))
add("cosine_invcdf", "p", r_cosine_invcdf, D(U(0, 1)), D(LU(1e-300, 0.5)))
add(
    "owens_t",
    "h a",
    r_owens_t,
    D(U(-5, 5), U(-5, 5)),
    D(LU(1e-5, 40, 0), LU(1e-5, 1e3, 0)),
)

section("Information Theory functions")
add("entr", "x", r_entr, D(U(0, 2)), D(LU(1e-300, 1e300)), note="zeros_c")
add(
    "rel_entr",
    "x y",
    r_rel_entr,
    D(LU(1e-3, 1e3), LU(1e-3, 1e3)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e300)),
    note="zeros_c",
)
add(
    "kl_div",
    "x y",
    r_kl_div,
    D(LU(1e-3, 1e3), LU(1e-3, 1e3)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e300)),
    note="zeros_c",
)
add("huber", "delta r", r_huber, D(U(0, 5), U(-10, 10)))
add(
    "pseudo_huber",
    "delta r",
    r_pseudo_huber,
    D(U(0, 5), U(-10, 10)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e300, 0)),
)

section("Gamma and related functions")
add("gamma", "x", r_gamma, D(U(-30, 30)), D(U(-170, 171.6)), D(LU(1e-300, 171)))
add("gammaln", "x", r_gammaln, D(U(-30, 30)), D(LU(1e-300, 1e300)), note="zeros")
add("rgamma", "x", r_rgamma, D(U(-30, 30)), D(U(-200, 200)))
add("digamma", "x", r_digamma, D(U(-30, 30)), D(LU(1e-300, 1e300)), note="zeros")
add("digammainv", "y", r_digammainv, D(U(-10, 10)), D(U(-1000, 700)))
add("trigamma", "x", r_trigamma, D(U(-30, 30)), D(LU(1e-300, 1e300)))
add(
    "beta",
    "a b",
    r_beta,
    D(U(0, 20), U(0, 20)),
    D(LU(1e-300, 170), LU(1e-300, 170)),
    D(LU(1e-3, 1e6, 0), LU(1e-3, 1e6)),
    note="neg",
)
add(
    "betaln",
    "a b",
    r_betaln,
    D(U(0, 20), U(0, 20)),
    D(LU(1e-300, 170), LU(1e-300, 170)),
    D(LU(1e-3, 1e6, 0), LU(1e-3, 1e6)),
    note="zeros neg",
)
add(
    "pow_rising",
    "a m",
    r_poch,
    D(U(0, 10), U(-10, 10)),
    D(LU(1e-3, 100), LU(1e-3, 100, 0)),
    D(LU(1e-3, 1e6, 0), LU(1e-3, 100, 0)),
    note="neg",
)
add("loggamma", "x", r_loggamma, D(U(0, 30)), D(LU(1e-300, 1e300)), note="zeros_c")
add("gammasgn", "x", r_gammasgn, D(U(-200, 200)))
add("lgam1p", "x", r_lgam1p, D(U(-0.9, 3)), D(LU(1e-300, 0.5, 0)))
add(
    "polygamma",
    "n x",
    r_polygamma,
    D(Int(1, 10), U(-10, 10)),
    D(Int(1, 50), LU(1e-3, 1e5)),
    cpp=False,
    note="zeros_c",
)
add(
    "multigammaln",
    "a d",
    r_multigammaln,
    D(
        sampler=lambda r: (
            (d := r.randint(1, 20)) / 2 - 0.5 + 10 ** r.uniform(-3, 3),
            d,
        ),
        label="d in {1, ..., 20}, a in (d - 1)/2 + [1e-3, 1e3] (log)",
        pretty=f"{nobreak('d', IN, Int(1, 20).pretty())}, "
        + nobreak("a", IN, f"(d {MINUS} 1)/2 + {LU(1e-3, 1e3).pretty()}")
        + " (log)",
    ),
    cpp=False,
    note="zeros_c",
)
add(
    "pow_falling",
    "x m",
    r_ff,
    D(U(0, 10), U(-10, 10)),
    D(LU(1e-3, 100), LU(1e-3, 100, 0)),
    cpp=False,
    note="neg",
)
for f, r in [
    ("gammainc", r_gammainc),
    ("gammaincc", r_gammaincc),
    ("log_gammainc", r_log_gammainc),
    ("log_gammaincc", r_log_gammaincc),
]:
    add(
        f,
        "a x",
        r,
        D(U(0, 20), U(0, 50)),
        D(LU(1e-300, 1e5), LU(1e-300, 1e5)),
        A_PM_SQRT_A,
        cpp=not f.startswith("log_"),
        timeout=60,
    )
for f, r in [("betainc", r_betainc), ("betaincc", r_betaincc)]:
    add(
        f,
        "a b x",
        r,
        D(U(0, 10), U(0, 10), U(0, 1)),
        D(LU(1e-5, 1e3), LU(1e-5, 1e3), U(0, 1)),
        AB_MU_PM_SIGMA,
        cpp=f == "betainc",
        note="beta_dist",
        timeout=60,
    )
for f, r in [("gammaincinv", r_gammaincinv), ("gammainccinv", r_gammainccinv)]:
    add(
        f,
        "a p" if f == "gammaincinv" else "a q",
        r,
        D(U(0, 20), U(0, 1), n=2000),
        D(LU(1e-3, 1e5), LU(1e-300, 1), n=2000),
        timeout=120,
        note="inv",
    )
for f, r in [("betaincinv", r_betaincinv), ("betainccinv", r_betainccinv)]:
    add(
        f,
        "a b y",
        r,
        D(U(0, 10), U(0, 10), U(0, 1), n=2000),
        D(LU(1e-3, 1e3), LU(1e-3, 1e3), LU(1e-300, 1), n=2000),
        cpp=f == "betaincinv",
        timeout=120,
        note="inv",
    )
for f, r in [
    ("gamma", mp.gamma),
    ("rgamma", mp.rgamma),
    ("loggamma", mp.loggamma),
    ("digamma", mp.digamma),
]:
    add(
        f"{f}#c",
        "re im",
        _c(r),
        CZ,
        complex=True,
        note="zeros_c" if f in ("loggamma", "digamma") else "",
    )
add("trigamma#c", "re im", _c(lambda z: mp.psi(1, z)), CZ, complex=True)

section("Error function and Fresnel integrals")
add("erf", "x", r_erf, D(U(-6, 6)), D(LU(1e-300, 10, 0)))
add("erfc", "x", r_erfc, D(U(-6, 27)), D(LU(1e-300, 27, 0)))
add("erfcx", "x", r_erfcx, D(U(-26, 26)), D(LU(1e-300, 1e300)))
add("erfinv", "y", r_erfinv, D(U(-1, 1)), D(Offset(1, LU(1e-16, 1), -1)))
add("erfcinv", "y", r_erfcinv, D(U(0, 2)), D(LU(1e-300, 1)))
add("erfi", "x", r_erfi, D(U(-6, 6)), D(LU(1e-300, 26.6, 0)))
add("dawsn", "x", r_dawsn, D(U(-10, 10)), D(LU(1e-300, 1e300, 0)))
add(
    "voigt_profile",
    "x sigma gamma",
    r_voigt_profile,
    D(U(-10, 10), U(0, 5), U(0, 5)),
    D(LU(1e-3, 1e3, 0), LU(1e-3, 1e3), LU(1e-3, 1e3)),
)
for f, r in [
    ("erf", mp.erf),
    ("erfc", mp.erfc),
    ("wofz", c_wofz),
    ("erfcx", c_erfcx),
    ("erfi", mp.erfi),
    ("dawsn", c_dawsn),
]:
    add(f"{f}#c", "re im", _c(r), D(U(-10, 10), U(-10, 10), n=5000), complex=True)
add("fresnel.0", "x", r_fresnel_s, D(U(-10, 10)), D(LU(1e-300, 1e300, 0)))
add("fresnel.1", "x", r_fresnel_c, D(U(-10, 10)), D(LU(1e-300, 1e300, 0)))
for i in range(2):
    add(
        f"fresnel.{i}#c",
        "re im",
        lambda re, im, i=i: c_fresnel(re, im)[i],
        D(U(-10, 10), U(-10, 10), n=5000),
        complex=True,
    )
for i, f in enumerate(["plus.0", "plus.1", "minus.0", "minus.1"]):
    add(
        f"modified_fresnel_{f}",
        "x",
        lambda x, i=i: c_modified_fresnel(x)[i],
        D(U(0, 10)),
        D(LU(1e-300, 1e4)),
        complex=True,
    )

add("erf_zeros", "k", _erf_zero, each_k(), complex=True, note="zeros_n")
for i in range(2):
    add(
        f"fresnel_zeros.{i}",
        "k",
        lambda k, i=i: _fresnel_zero(i, k),
        each_k(),
        complex=True,
        note="zeros_k",
    )

section("Legendre functions")


def nm_sampler(r: Random, hi: int = 50) -> tuple[int, int]:
    n = r.randint(0, hi)
    return n, r.randint(-n, n)


def nm_nonneg(r: Random, hi: int) -> tuple[int, int]:
    """`nm_sampler` with m in {0, ..., n}."""
    n = r.randint(0, hi)
    return n, r.randint(0, n)


def legendre_note(f: str) -> str:
    """The notes of `f`; an `_all` row follows the row of its scalar function, with the same
    domains and reference."""
    return "zeros_c all" if f.endswith("_all") else "zeros_c"


NM = f"{nobreak('n', IN, Int(0, 50).pretty())}, {nobreak('m', IN, f'{{{MINUS}n, {ELLIPSIS}, n}}')}"
for f in ["legendre_p", "legendre_p_all"]:
    add(
        f,
        "n x",
        r_legendre,
        D(Int(0, 100), U(-1, 1)),
        D(Int(0, 30), LU(1, 1e3, 0)),
        note=legendre_note(f),
    )
for f, r in [
    ("assoc_legendre_p", r_assoc_legendre_p),
    ("assoc_legendre_p_all", r_assoc_legendre_p),
    ("assoc_legendre_p_norm", r_assoc_legendre_p_norm),
    ("assoc_legendre_p_norm_all", r_assoc_legendre_p_norm),
]:
    add(
        f,
        "n m x",
        r,
        D(
            sampler=lambda r: (*nm_sampler(r), r.uniform(-1, 1)),
            label="n in {0, ..., 50}, m in {-n, ..., n}, x in [-1, 1]",
            pretty=f"{NM}, {nobreak('x', IN, U(-1, 1).pretty())}",
        ),
        note=legendre_note(f),
    )
for f in ["sph_legendre_p", "sph_legendre_p_all"]:
    add(
        f,
        "n m theta",
        r_sph_legendre_p,
        D(
            sampler=lambda r: (*nm_sampler(r), r.uniform(0, math.pi)),
            label="n in {0, ..., 50}, m in {-n, ..., n}, theta in [0, pi]",
            pretty=f"{NM}, {nobreak(GREEK['theta'], IN, U(0, math.pi).pretty())}",
        ),
        note=legendre_note(f),
    )
for f in ["sph_harm_y", "sph_harm_y_all"]:
    add(
        f,
        "n m theta phi",
        c_sph_harm_y,
        D(
            sampler=lambda r: (
                *nm_sampler(r),
                r.uniform(0, math.pi),
                r.uniform(0, 2 * math.pi),
            ),
            label="n in {0, ..., 50}, m in {-n, ..., n}, theta in [0, pi], phi in [0, 2 pi]",
            pretty=f"{NM}, {nobreak(GREEK['theta'], IN, U(0, math.pi).pretty())}, "
            + nobreak(GREEK["phi"], IN, "[0, 2\u03c0]"),
            n=5000,
        ),
        complex=True,
        note=legendre_note(f),
    )
add("pmv", "m v x", r_pmv, D(Int(0, 10), U(-50, 50), U(-1, 1), n=5000), note="zeros_c")
for i, r in enumerate([r_legendre_q, r_legendre_qp]):
    add(
        f"legendre_q_all.{i}",
        "n x",
        r,
        D(Int(0, 100), U(-1, 1)),
        D(Int(0, 30), LU(1, 1e3, 0)),
        note="zeros_c all",
    )
for i, r in enumerate([r_assoc_legendre_q, r_assoc_legendre_qp]):
    add(
        f"assoc_legendre_q_all.{i}",
        "n m x",
        r,
        D(
            sampler=lambda r: (*nm_nonneg(r, 50), r.uniform(-1, 1)),
            label="n in {0, ..., 50}, m in {0, ..., n}, x in [-1, 1]",
            pretty=f"{nobreak('n', IN, Int(0, 50).pretty())}, "
            + f"{nobreak('m', IN, f'{{0, {ELLIPSIS}, n}}')}, "
            + nobreak("x", IN, U(-1, 1).pretty()),
        ),
        D(
            sampler=lambda r: (*nm_nonneg(r, 30), LU(1, 1e3, 0).sample(r)),
            label="n in {0, ..., 30}, m in {0, ..., n}, x in +-[1, 1e3] (log)",
            pretty=f"{nobreak('n', IN, Int(0, 30).pretty())}, "
            + f"{nobreak('m', IN, f'{{0, {ELLIPSIS}, n}}')}, "
            + f"{nobreak('x', IN, LU(1, 1e3, 0).pretty())} (log)",
        ),
        note="zeros_c all",
    )

section("Orthogonal polynomials")


def poly(
    name: str,
    args: str,
    ref: FunctionType,
    dists: list[Dist],
    int_n: int = 50,
    cpp: bool = True,
) -> None:
    """A spec with an integer-degree domain (via the `_l` variant), and a real-degree one."""
    add(
        name,
        args,
        ref,
        D(Int(0, int_n), *dists, call=f"{name}_l"),
        D(U(0, int_n), *dists),
        cpp=cpp,
        note="zeros_c poly",
    )


poly("eval_jacobi", "n alpha beta x", r_jacobi, [U(-1, 10), U(-1, 10), U(-1, 1)])
poly(
    "eval_sh_jacobi",
    "n alpha beta x",
    r_sh_jacobi,
    [U(-1, 10), U(-1, 10), U(0, 1)],
    cpp=False,
)
poly("eval_gengegenbauer", "n p q x", r_gengegenbauer, [U(1, 10), U(0.5, 1), U(0, 1)])
poly("eval_gegenbauer", "n alpha x", r_gegenbauer, [U(-0.5, 10), U(-1, 1)])
poly("eval_chebyshev_t", "n x", r_chebyt, [U(-1, 1)], 100)
poly("eval_chebyshev_u", "n x", r_chebyu, [U(-1, 1)], 100)
poly("eval_chebyshev_c", "n x", r_chebyc, [U(-2, 2)], 100)
poly("eval_chebyshev_s", "n x", r_chebys, [U(-2, 2)], 100)
poly("eval_sh_chebyshev_t", "n x", r_sh_chebyt, [U(0, 1)], 100)
poly("eval_sh_chebyshev_u", "n x", r_sh_chebyu, [U(0, 1)], 100)
poly("eval_legendre", "n x", r_legendre, [U(-1, 1)], 100)
poly("eval_sh_legendre", "n x", r_sh_legendre, [U(0, 1)], 100)
poly("eval_genlaguerre", "n alpha x", r_genlaguerre, [U(-1, 10), U(0, 50)])
poly("eval_laguerre", "n x", r_laguerre, [U(0, 50)])
for f, r in [("eval_hermite_h", r_hermite_h), ("eval_hermite_he", r_hermite_he)]:
    add(
        f,
        "n x",
        r,
        D(Int(0, 50), U(-10, 10), call=f"{f}_l"),
        D(U(0, 50), U(-10, 10), cpp=False),
        note="zeros_c poly",
    )


def hahn_domain(name: str, pretty: str, x: Callable[[Random, int], float]) -> Domain:
    """N in {1, ..., 1000}, n in {0, ..., N}, alpha and beta in [-1, 100], and x(rng, N)."""

    def sample(r: Random) -> tuple[float, ...]:
        big_n = r.randint(1, 1000)
        a, b = r.uniform(-1, 100), r.uniform(-1, 100)
        return float(r.randint(0, big_n)), float(big_n), a, b, x(r, big_n)

    return D(
        sampler=sample,
        label=f"N in {{1, ..., 1000}}, n in {{0, ..., N}}, alpha, beta in [-1, 100], x {name}",
        pretty=f"{nobreak('N', IN, Int(1, 1000).pretty())}, "
        + f"{nobreak('n', IN, f'{{0, {ELLIPSIS}, N}}')}, "
        + f"{nobreak(GREEK['alpha'] + ',', GREEK['beta'], IN, U(-1, 100).pretty())}, "
        + nobreak("x", IN, pretty),
        n=2000,
        cpp=False,
    )


add(
    "eval_hahn",
    "n N alpha beta x",
    r_hahn_q,
    hahn_domain(
        "in {0, ..., N}", f"{{0, {ELLIPSIS}, N}}", lambda r, m: float(r.randint(0, m))
    ),
    hahn_domain("in [0, N]", "[0, N]", lambda r, m: r.uniform(0, m)),
    hahn_domain(
        "in [-N, 0] or [N, 2N]",
        f"[{MINUS}N, 0] \u222a [N, 2N]",
        lambda r, m: r.uniform(-m, 0) if r.random() < 0.5 else r.uniform(m, 2 * m),
    ),
    note="zeros_c",
)

NK100 = Int(1, 100)
GG = D(
    sampler=lambda r: (
        n := r.randint(1, 100),
        r.randint(1, n),
        (q := r.uniform(0.5, 10)) + r.uniform(-1, 10),
        q,
    ),
    label="n in {1, ..., 100}, k in {1, ..., n}, p - q in [-1, 10], q in [0.5, 10]",
    pretty=f"{nobreak('n', IN, NK100.pretty())}, "
    + f"{nobreak('k', IN, f'{{1, {ELLIPSIS}, n}}')}, "
    + f"{nobreak(f'p {MINUS} q', IN, U(-1, 10).pretty())}, "
    + nobreak("q", IN, U(0.5, 10).pretty()),
    n=1000,
)
AB = {"alpha": U(-1, 10), "beta": U(-1, 10)}
for f, args, dom in [
    ("jacobi", "alpha beta", nk_domain(NK100, AB)),
    ("sh_jacobi", "alpha beta", nk_domain(NK100, AB)),
    ("gengegenbauer", "p q", GG),
    ("legendre", "", nk_domain(ILU(1, 1000), {})),
    ("sh_legendre", "", nk_domain(ILU(1, 1000), {})),
    *(
        (f"{s}chebyshev_{t}", "", nk_domain(ILU(1, 1000), {}, 2000))
        for s in ["", "sh_"]
        for t in ["t", "u", "c", "s"]
        if not (s and t in "cs")
    ),
    ("gegenbauer", "alpha", nk_domain(NK100, {"alpha": U(-0.5, 10)})),
    ("genlaguerre", "alpha", nk_domain(NK100, {"alpha": U(-1, 10)})),
    ("laguerre", "", nk_domain(ILU(1, 1000), {})),
    ("hermite_h", "", nk_domain(ILU(1, 1000), {})),
    ("hermite_he", "", nk_domain(ILU(1, 1000), {})),
]:
    for i in range(2):
        add(
            f"roots_{f}.{i}",
            f"n k {args}".strip(),
            lambda *a, f=f, i=i: r_gauss(f, i, *a),
            dom,
            cpp=False,
            note="gauss",
        )

section("Hypergeometric functions")
add("hyp0f0", "z", r_hyp0f0, D(U(-700, 700)), cpp=False)
add("hyp1f0", "a z", r_hyp1f0, D(U(-10, 10), U(-10, 1)), cpp=False)
add("hyp0f1", "b z", r_hyp0f1, D(LU(1e-3, 1e6, 0), LU(1e-3, 1e6, 0)), note="osc_hyp0f1")
add(
    "hyp0f0#c",
    "re im",
    _c(mp.exp),
    D(U(-700, 700), U(-700, 700), n=5000),
    complex=True,
    cpp=False,
)
add(
    "hyp1f0#c",
    "a re im",
    c_hyp1f0,
    D(U(-10, 10), U(-10, 10), U(-10, 10), n=5000),
    complex=True,
    cpp=False,
)
add(
    "hyp0f1#c",
    "b re im",
    c_hyp0f1,
    D(LU(1e-3, 100, 0), U(-50, 50), U(-50, 50), n=5000),
    D(LU(100, 1e4, 0), U(-50, 50), U(-50, 50), n=2000),
    complex=True,
    note="osc_hyp0f1_c",
)
add(
    "hyp1f1",
    "a b x",
    r_hyp1f1,
    D(U(-20, 20), U(-20, 20), U(-50, 50), n=5000),
    D(LU(1e-3, 1e3, 0), LU(1e-3, 1e3, 0), LU(1e-3, 1e3, 0), n=2000),
    timeout=60,
    note="zeros_c",
)
add(
    "hypu",
    "a b x",
    r_hypu,
    D(U(-10, 10), U(-10, 10), U(0, 50), n=2000),
    D(LU(1e-3, 100, 0), LU(1e-3, 100, 0), LU(1e-3, 1e3), n=2000),
    timeout=60,
    note="zeros_c",
)
add(
    "hyp2f1",
    "a b c x",
    r_hyp2f1,
    D(U(-10, 10), U(-10, 10), U(-10, 10), U(-1, 1), n=5000),
    D(LU(1e-3, 100, 0), LU(1e-3, 100, 0), LU(1e-3, 100, 0), U(-10, 1), n=2000),
    D(LU(1e-3, 1e3, 0), LU(1e-3, 1e3, 0), LU(1e-3, 1e3, 0), U(-10, 1), n=1000),
    timeout=60,
    note="zeros_c",
)
add(
    "hyp3f2",
    "a1 a2 a3 b1 b2 z",
    r_hyppfq,
    D(*[U(-10, 10)] * 5, U(-1, 1), n=5000),
    D(*[LU(1e-3, 100, 0)] * 5, U(-10, 1), n=2000),
    D(Int(-30, -1), *[U(-10, 10)] * 4, U(-10, 10), n=2000),
    cpp=False,
    timeout=60,
    note="none",
)
add(
    "hyp2f2",
    "a1 a2 b1 b2 z",
    r_hyp2f2,
    D(*[U(-10, 10)] * 4, U(-30, 30), n=5000),
    D(*[U(-10, 10)] * 4, LU(1e-3, 1e3, 0), n=2000),
    D(*[LU(1e-3, 100, 0)] * 4, LU(1e-3, 1e3, 0), n=2000),
    D(Int(-30, -1), *[U(-10, 10)] * 3, U(-30, 30), n=1000),
    cpp=False,
    timeout=60,
    note="none",
)
add(
    "hyp1f2",
    "a b1 b2 z",
    r_hyp1f2,
    D(*[U(-10, 10)] * 3, U(-225, 225), n=5000),
    D(*[U(-10, 10)] * 3, LU(1e-3, 1e6, 0), n=2000),
    D(*[LU(1e-3, 100, 0)] * 3, LU(1e-3, 1e6, 0), n=2000),
    D(Int(-30, -1), *[U(-10, 10)] * 2, U(-225, 225), n=1000),
    cpp=False,
    timeout=60,
    note="none",
)
for name, args in [
    ("hyppfq#4f3", "a1 a2 a3 a4 b1 b2 b3 z"),
    ("hyppfq#5f4", "a1 a2 a3 a4 a5 b1 b2 b3 b4 z"),
]:
    params = len(args.split()) - 1
    add(
        name,
        args,
        r_hyppfq,
        D(*[U(-10, 10)] * params, U(-1, 1), n=2000),
        D(*[LU(1e-3, 100, 0)] * params, U(-10, 1), n=2000),
        cpp=False,
        timeout=60,
        note="none",
    )
add(
    "hyppfq#2f0",
    "a1 a2 z",
    r_hyppfq_arb("hyppfq#2f0"),
    D(U(-10, 10), U(-10, 10), LU(1e-3, 1e6, -1), n=2000, call="hyppfq2"),
    D(LU(1e-3, 100, 0), LU(1e-3, 100, 0), LU(1e-3, 1e6, -1), n=2000, call="hyppfq2"),
    cpp=False,
    timeout=60,
    note="none",
)
for p, args, zmax in [
    (3, "a1 a2 a3 b1 b2 b3 z", 1e3),
    (2, "a1 a2 b1 b2 b3 z", 1e6),
    (1, "a b1 b2 b3 z", 1e8),
    (0, "b1 b2 z", 1e8),
    (0, "b1 b2 b3 z", 1e8),
]:
    params, call = len(args.split()) - 1, f"hyppfq{p}"
    name = f"hyppfq#{p}f{params - p}"
    add(
        name,
        args,
        r_hyppfq_arb(name),
        D(*[U(-10, 10)] * params, LU(1e-3, zmax, 0), n=2000, call=call),
        D(*[LU(1e-3, 100, 0)] * params, LU(1e-3, zmax, 0), n=2000, call=call),
        cpp=False,
        timeout=60,
        note="none",
    )
add(
    "hyp1f1#c",
    "a b re im",
    c_hyp1f1,
    D(U(-10, 10), U(-10, 10), U(-20, 20), U(-20, 20), n=2000),
    complex=True,
    timeout=60,
)
add(
    "hyp2f1#c",
    "a b c re im",
    c_hyp2f1,
    D(U(-5, 5), U(-5, 5), U(-5, 5), U(-2, 2), U(-2, 2), n=2000),
    complex=True,
    timeout=60,
)

section("Parabolic cylinder functions")
for f, r in [("pbdv", r_pbdv), ("pbvv", r_pbvv)]:
    for seq in ["", "_seq"]:
        for i in range(2):
            add(
                f"{f}{seq}.{i}",
                "v x",
                lambda v, x, r=r, i=i: r(v, x)[i],
                # the sequences have non-negative orders
                D(U(-10 if not seq else 0, 10), U(-10, 10), n=2000),
                D(LU(20, 300, 0 if not seq else 1), U(-60, 60), n=2000),
                note="zeros_c" + (" seq" if seq else ""),
            )
for i in range(2):
    add(
        f"pbwa.{i}",
        "a x",
        lambda a, x, i=i: r_pbwa(a, x)[i],
        D(U(-5, 5), U(-5, 5), n=2000),
        note="zeros_c",
    )
for i in range(2):
    add(
        f"pbdn_seq.{i}",
        "n re im",
        lambda n, re, im, i=i: c_pbdn(n, re, im)[i],
        D(Int(0, 20), U(-10, 10), U(-10, 10), n=2000),
        D(Int(-20, -1), U(-10, 10), U(-10, 10), n=2000),
        D(Int(0, 100), U(-30, 30), U(-30, 30), n=1000),
        complex=True,
        note="seq",
    )

section("Mathieu functions")


def mathieu_coef_domain(odd: int, q: Dist) -> Domain:
    """(m, q, k) with k the index of one of the first floor(m / 2) + 7 coefficients, which the
    output has for every q."""
    m = Int(odd, 30)

    def sampler(r: Random) -> tuple[float, float, int]:
        order = m.sample(r)
        return order, q.sample(r), r.randint(0, int(order) // 2 + 6)

    log = " (log)" if _is_log(q) else ""
    return D(
        sampler=sampler,
        label=f"m in {m.desc()}, q in {q.desc()}, k in {{0, ..., floor(m / 2) + 6}}",
        pretty=f"{nobreak('m', IN, m.pretty())}, {nobreak('q', IN, q.pretty())}{log}, "
        + nobreak("k", IN, f"{{0, {ELLIPSIS}, \u230am/2\u230b + 6}}"),
        n=2000,
    )


def mathieu_mod_domain(odd: int, hi: float, osc: bool) -> Domain:
    """(m, q, x) with q in [0, hi], where the modified Mathieu functions of order m oscillate,
    2 q cosh(2 x) > m^2 (about the characteristic value), or where they do not."""
    m, q, x = Int(odd, 30), U(0, hi), U(0, 3)

    def sampler(r: Random) -> tuple[float, float, float]:
        while True:
            a = m.sample(r), q.sample(r), x.sample(r)
            if (2 * a[1] * math.cosh(2 * a[2]) > a[0] ** 2) == osc:
                return a

    cond = f"2q cosh 2x {'>' if osc else '<'} m\u00b2"
    return D(
        sampler=sampler,
        label=f"m in {m.desc()}, q in {q.desc()}, x in {x.desc()}, {cond}",
        pretty=f"{nobreak('m', IN, m.pretty())}, {nobreak('q', IN, q.pretty())}, "
        + f"{nobreak('x', IN, x.pretty())}, {nobreak(cond)}",
        n=2000,
    )


for odd, f, r in [(0, "mathieu_a", r_mathieu_a), (1, "mathieu_b", r_mathieu_b)]:
    add(
        f,
        "m q",
        r,
        D(Int(odd, 30), U(0, 10)),
        D(Int(odd, 30), U(0, 100)),
        D(Int(odd, 30), LU(1e-300, 1, 0)),
        D(Int(odd, 30), LU(100, 1e4, 0)),
        note="mathieu mathieu_cv zeros_c",
    )
for odd, f, r in [
    (0, "mathieu_cem", [r_mathieu_cem, r_mathieu_cemp]),
    (1, "mathieu_sem", [r_mathieu_sem, r_mathieu_semp]),
]:
    for i, ref in enumerate(r):
        add(
            f"{f}.{i}",
            "m q x",
            ref,
            D(Int(odd, 30), U(0, 10), U(-360, 360), n=2000),
            D(Int(odd, 30), U(0, 100), U(-360, 360), n=2000),
            D(Int(odd, 30), LU(1e-300, 1, 0), U(-360, 360), n=2000),
            D(Int(odd, 30), LU(100, 2000, 0), U(-360, 360), n=2000),
            note="mathieu_osc",
        )
for odd, f, r in [
    (0, "mathieu_ce", [r_mathieu_ce, r_mathieu_cep]),
    (1, "mathieu_se", [r_mathieu_se, r_mathieu_sep]),
]:
    for i, ref in enumerate(r):
        add(
            f"{f}.{i}",
            "m q x",
            ref,
            D(Int(odd, 30), U(0, 10), U(-2 * math.pi, 2 * math.pi), n=2000),
            D(Int(odd, 30), U(0, 100), U(-2 * math.pi, 2 * math.pi), n=2000),
            cpp=False,
            note="mathieu_osc",
        )
for odd, f in enumerate(["mathieu_modcem", "mathieu_modsem"]):
    for kind in (1, 2):
        for i in range(2):
            add(
                f"{f}{kind}.{i}",
                "m q x",
                lambda m, q, x, odd=odd, kind=kind, i=i: _mathieu_radial(
                    odd, kind, m, q, x, bool(i)
                ),
                mathieu_mod_domain(odd, 10, True),
                mathieu_mod_domain(odd, 10, False),
                mathieu_mod_domain(odd, 100, True),
                D(Int(odd, 30), LU(1e-10, 1e3), U(0, 3), n=2000),
                note="mathieu_mod",
            )
for odd, f, r in [
    (0, "mathieu_even_coef", r_mathieu_even_coef),
    (1, "mathieu_odd_coef", r_mathieu_odd_coef),
]:
    add(
        f,
        "m q k",
        r,
        mathieu_coef_domain(odd, U(0, 10)),
        mathieu_coef_domain(odd, U(0, 100)),
        mathieu_coef_domain(odd, LU(1e-7, 1e3)),
        mathieu_coef_domain(odd, LU(1e-20, 1e-7)),
        note="mathieu_coef zeros_c",
    )

section("Spheroidal wave functions")
SW_N = 20  # the largest n; must match `SEGV_N` of the dump programs


def sw_domain(kd: int, c: Dist, x: Dist | None = None, cv: bool = False) -> Domain:
    """(m, n, c, x) with m <= n <= SW_N: without x if it is None, and with the characteristic
    value rounded to f64 after c if `cv`."""

    def sampler(r: Random) -> list[float]:
        n, m = nm_nonneg(r, SW_N)
        args = [m, n, c.sample(r)]
        if cv:
            with mp.workdps(30):
                args.append(float(_sw_cv(m, n, kd * mp.mpf(args[2]) ** 2)[0]))
        return args if x is None else [*args, x.sample(r)]

    dists = [("c", c)] if x is None else [("c", c), ("x", x)]
    return D(
        sampler=sampler,
        label=f"n in {{0, ..., {SW_N}}}, m in {{0, ..., n}}, "
        + ", ".join(f"{a} in {d.desc()}" for a, d in dists),
        pretty=f"{nobreak('n', IN, Int(0, SW_N).pretty())}, "
        + f"{nobreak('m', IN, f'{{0, {ELLIPSIS}, n}}')}, "
        + ", ".join(
            nobreak(a, IN, d.pretty()) + (" (log)" if _is_log(d) else "")
            for a, d in dists
        ),
        n=2000,
    )


def sw_ref(r: FunctionType, kd: int, i: int, cv: bool) -> FunctionType:
    """Output i of `r(m, n, c, x, kd)` as a function of (m, n, c, x), or of (m, n, c, cv, x)."""
    if cv:
        return lambda m, n, c, cv, x: r(m, n, c, x, kd)[i]
    return lambda m, n, c, x: r(m, n, c, x, kd)[i]


SW_C = (U(0, 10), U(10, 50))
SW_KD = (("prolate", 1), ("oblate", -1))
for seq in ("", "_seq"):
    for f, kd in SW_KD:
        add(
            f"{f}_segv{seq}",
            "m n c",
            lambda m, n, c, kd=kd: _sw_cv(int(m), int(n), kd * c * c)[0],
            *(sw_domain(kd, c) for c in SW_C),
            note=" ".join(k for k, on in [("sw_seq", seq), ("zeros_c", kd < 0)] if on),
        )
for f, kd in SW_KD:
    x_rad = Offset(1, LU(1e-3, 10)) if kd == 1 else LU(1e-3, 10)
    for g, r, x, note in [
        ("aswfa", _sw_angular, U(-1, 1), "sw_ang"),
        ("radial1", _sw_rad1, x_rad, "sw_rad"),
        ("radial2", _sw_rad2, x_rad, "sw_rad"),
    ]:
        for cv in (False, True):
            for i in range(2):
                add(
                    f"{f}_{g}{'' if cv else '_nocv'}.{i}",
                    "m n c cv x" if cv else "m n c x",
                    sw_ref(r, kd, i, cv),
                    *(sw_domain(kd, c, x, cv) for c in SW_C),
                    note=f"{note} sw_cv" if cv else note,
                )

section("Kelvin functions")
for i, f in enumerate(["ber", "bei", "ker", "kei", "berp", "beip", "kerp", "keip"]):
    pos = f.startswith("k")
    add(
        f,
        "x",
        lambda x, i=i: _kelvin(i, x),
        D(U(0, 20) if pos else U(-20, 20)),
        D(LU(1e-300, 1e3, 1 if pos else 0)),
        note="kelvin",
    )

for i, f in enumerate(["ber", "bei", "ker", "kei", "berp", "beip", "kerp", "keip"]):
    add(
        f"{f}_zeros",
        "k",
        lambda k, i=i: _kelvin_zero(i, k),
        each_k(200),
        note="zeros_k",
    )

section("Combinatorics")
add("perm", "n k", r_ff, D(U(0, 100), Int(0, 30)), cpp=False)
add(
    "stirling2",
    "n k",
    r_stirling2,
    D(
        sampler=lambda r: (n := r.randint(0, 150), r.randint(0, n)),
        label="n in {0, ..., 150}, k in {0, ..., n}",
        pretty=f"{nobreak('n', IN, Int(0, 150).pretty())}, {nobreak('k', IN, f'{{0, {ELLIPSIS}, n}}')}",
    ),
)

section("Exponential integrals")
add("expi", "x", r_expi, D(U(-20, 20)), D(LU(1e-300, 700, 0)))
add("expn", "n x", r_expn, D(Int(0, 10), U(0, 50)), D(Int(0, 1000), LU(1e-300, 700)))
add("exp1", "x", r_exp1, D(U(0, 20)), D(LU(1e-300, 750)))
add("scaled_exp1", "x", r_scaled_exp1, D(U(0, 20)), D(LU(1e-300, 1e300)))
for f, r in [("exp1", mp.e1), ("expi", mp.ei)]:
    add(f"{f}#c", "re im", _c(r), CZ, complex=True)

section("Zeta functions")
add(
    "zeta",
    "x q",
    r_hurwitz_zeta,
    D(U(1, 20), U(0, 10)),
    D(Offset(1, LU(1e-15, 1e3)), LU(1e-300, 1e300)),
)
add(
    "riemann_zeta",
    "x",
    r_riemann_zeta,
    D(U(-20, 20)),
    D(LU(1e-300, 1e3, 0)),
    note="zeros",
)
add(
    "riemann_zeta#c",
    "re im",
    _c(mp.zeta),
    CZ,
    D(U(-20, 20), LU(20, 1e3, 0), n=5000),
    complex=True,
)
add("zetac", "x", r_zetac, D(U(-20, 20)), D(LU(1e-300, 1e3, 0)), note="zeros")

section("Harmonic numbers", 3)
add(
    "gen_harmonic",
    "n a",
    r_gen_harmonic,
    D(Int(1, 1000), U(-5, 5)),
    D(ILU(1, 10**15), Offset(1, LU(1e-2, 20))),
)


def ngh_sampler(r: Random) -> tuple[int, int, int, float]:
    n = r.randint(1, 1000)
    j = r.randint(1, n)
    return j, r.randint(j, n), n, r.uniform(-5, 5)


add(
    "normalized_gen_harmonic",
    "j k n a",
    r_normalized_gen_harmonic,
    D(
        sampler=ngh_sampler,
        label="n in {1, ..., 1000}, 1 <= j <= k <= n, a in [-5, 5]",
        pretty=f"{nobreak('n', IN, Int(1, 1000).pretty())}, 1 \u2264 j \u2264 k \u2264 n, "
        + nobreak("a", IN, U(-5, 5).pretty()),
        n=2000,
    ),
)

section("Lambert W and related functions")
add("lambertw", "x", r_lambertw, D(U(INV_E, 10)), D(LU(1e-300, 1e300)))
add("wrightomega", "x", r_wrightomega, D(U(-20, 20)), D(LU(1e-300, 1e300, 0)))
add("lambertw#c", "re im", _c(c_lambertw), CZ, complex=True)

section("Other special functions")
add("agm", "a b", r_agm, D(U(0, 10), U(0, 10)), D(LU(1e-300, 1e300), LU(1e-300, 1e300)))
add(
    "binom",
    "n k",
    r_binom,
    D(
        sampler=lambda r: (n := r.uniform(0, 100), r.uniform(0, n)),
        label="n in [0, 100], k in [0, n]",
        pretty=f"{nobreak('n', IN, '[0, 100]')}, {nobreak('k', IN, '[0, n]')}",
    ),
    D(
        sampler=lambda r: (n := r.randint(0, 1000), r.randint(0, n)),
        label="n in {0, ..., 1000}, k in {0, ..., n}",
        pretty=f"{nobreak('n', IN, f'{{0, {ELLIPSIS}, 1000}}')}, "
        + nobreak("k", IN, f"{{0, {ELLIPSIS}, n}}"),
    ),
)
add("diric", "n x", r_diric, D(Int(1, 20), U(-10, 10)), cpp=False, note="zeros_c")
add("sici.0", "x", r_si, D(U(-30, 30)), D(LU(1e-300, 1e300, 0)))
add("sici.1", "x", r_ci, D(U(-30, 30)), D(LU(1e-300, 1e300, 0)), note="zeros_c")
add("shichi.0", "x", r_shi, D(U(-30, 30)), D(LU(1e-300, 1e3, 0)))
add("shichi.1", "x", r_chi, D(U(-30, 30)), D(LU(1e-300, 1e3, 0)), note="zeros_c")
add("spence", "x", r_spence, D(U(0, 2)), D(LU(1e-300, 1e300)))
add("softplus", "x", r_softplus, D(U(-20, 20)), D(LU(1e-300, 1e300, 0)), cpp=False)
add("log1mexp", "x", r_log1mexp, D(U(-20, 0)), D(LU(1e-300, 1e300, -1)))
for f, r in [("bernoulli", r_bernoulli), ("euler", r_euler)]:
    add(f, "n", r, each_k(ZEROS_N - 1, 0), note="numbers")

section("Convenience functions")
add("cbrt", "x", r_cbrt, D(LU(1e-300, 1e300, 0)), std=True)
add("exp10", "x", r_exp10, D(U(-1, 1)), D(U(-307, 308)))
add("exp2", "x", r_exp2, D(U(-1, 1)), D(U(-1022, 1023)), std=True)
add("radian", "d m s", r_radian, D(U(0, 360), U(0, 60), U(0, 60)))
for f, r in [
    ("sindg", r_sindg),
    ("cosdg", r_cosdg),
    ("tandg", r_tandg),
    ("cotdg", r_cotdg),
]:
    add(f, "x", r, D(U(-360, 360)), D(LU(1e-300, 1e14, 0)))
add("cospi", "x", r_cospi, D(U(-10, 10)), D(LU(1e-300, 1e15, 0)))
add("sinpi", "x", r_sinpi, D(U(-10, 10)), D(LU(1e-300, 1e15, 0)))
for f, r in [("cospi", mp.cospi), ("sinpi", mp.sinpi)]:
    add(
        f"{f}#c",
        "re im",
        _c(r),
        CZ,
        D(LU(1e-300, 1e15, 0), LU(1e-300, 100, 0), n=5000),
        complex=True,
        note="zeros_c",
    )
add(
    "log1p",
    "x",
    r_log1p,
    D(U(-1, 1)),
    D(LU(1e-300, 1e300)),
    D(LU(1e-300, 1, -1)),
    std=True,
)
add("log1pmx", "x", r_log1pmx, D(U(-1, 2)), D(LU(1e-300, 1, 0)), D(LU(1, 1e300)))
add("expm1", "x", r_expm1, D(U(-2, 2)), D(LU(1e-300, 709, 0)), std=True)
for f, r in [("log1p", mp.log1p), ("expm1", mp.expm1)]:
    add(
        f"{f}#c",
        "re im",
        _c(r),
        CZ,
        D(LU(1e-300, 1, 0), LU(1e-300, 1, 0), n=5000),
        complex=True,
    )
add(
    "cosm1",
    "x",
    r_cosm1,
    D(U(-math.pi, math.pi)),
    D(LU(1e-300, 1e3, 0)),
    note="zeros_c",
)
add(
    "powm1",
    "x y",
    r_powm1,
    D(LU(1e-3, 1e3), U(-5, 5)),
    D(Offset(1, LU(1e-15, 0.1, 0)), LU(1e-3, 1e3, 0)),
    D(LU(1e-300, 1e300), LU(1e-300, 1e-3, 0)),
    D(U(-10, 0), Int(-50, 50)),
    cpp=False,
)
add(
    "xlogy",
    "x y",
    r_xlogy,
    D(LU(1e-300, 1e300, 0), LU(1e-300, 1e300)),
    D(U(-10, 10), Offset(1, LU(1e-16, 1e-3, 0))),
)
add(
    "xlog1py",
    "x y",
    r_xlog1py,
    D(LU(1e-300, 1e300, 0), LU(1e-300, 1e300)),
    D(U(-10, 10), U(-1, 1)),
)
add(
    "xlogy#c",
    "xr xi yr yi",
    c_xlogy,
    D(U(-10, 10), U(-10, 10), U(-10, 10), U(-10, 10), n=5000),
    D(
        U(-10, 10),
        U(-10, 10),
        Offset(1, LU(1e-16, 1e-3, 0)),
        LU(1e-16, 1e-3, 0),
        n=5000,
    ),
    complex=True,
)
add(
    "xlog1py#c",
    "xr xi yr yi",
    c_xlog1py,
    D(U(-10, 10), U(-10, 10), U(-10, 10), U(-10, 10), n=5000),
    D(U(-10, 10), U(-10, 10), LU(1e-300, 1, 0), LU(1e-300, 1, 0), n=5000),
    complex=True,
)
for f, r in [("logaddexp", r_logaddexp), ("logaddexp2", r_logaddexp2)]:
    add(
        f,
        "x y",
        r,
        D(U(-100, 100), U(-100, 100)),
        D(LU(1e-300, 1e300, 0), LU(1e-300, 1e300, 0)),
        cpp=False,
        note="zeros_c",
    )
add("exprel", "x", r_exprel, D(U(-2, 2)), D(LU(1e-300, 709, 0)))
add(
    "sinc",
    "x",
    r_sinc,
    D(U(-10, 10)),
    D(LU(1e-300, 1e15, 0)),
    cpp=False,
    note="zeros_c",
)

section("Functions of arrays")


def elements(r: Random, dist: Dist, n: int, inf: tuple[float, ...]) -> list[float]:
    """`n` samples of `dist`, each replaced by one of `inf` with probability 1/4, but with at
    least one finite."""
    xs = [
        r.choice(inf) if inf and r.random() < 0.25 else dist.sample(r) for _ in range(n)
    ]
    if not any(map(math.isfinite, xs)):
        xs[r.randrange(n)] = dist.sample(r)
    return xs


def array(name: str, dist: Dist, inf: str = "", length: str = "\u207f") -> str:
    """The pretty domain of an array of `length` elements of `dist`, or of the set `inf`."""
    s = dist.pretty() + (f" \u222a {{{inf}}}" if inf else "")
    if inf or not s.startswith(("[", f"{PM}[")):
        s = f"({s})"
    return nobreak(name, IN, s + length) + (" (log)" if _is_log(dist) else "")


def xs_domain(dist: Dist, inf: tuple[float, ...] = (), index: bool = False) -> Domain:
    """An array x of n in {1, ..., 20} elements (see `elements`), preceded by a random index j
    into it if `index`."""

    def sampler(r: Random) -> list[float]:
        n = r.randint(1, 20)
        xs = elements(r, dist, n, inf)
        return [r.randrange(n), *xs] if index else xs

    sign, pretty = ("+-", PM) if math.inf in inf else ("-", MINUS)
    return D(
        sampler=sampler,
        label=f"x in {dist.desc()}^n"
        + (f" or {sign}inf (1/4)" if inf else "")
        + ", n in {1, ..., 20}"
        + (", j in {0, ..., n - 1}" if index else ""),
        pretty=f"{array('x', dist, f'{pretty}\u221e' if inf else '')}, "
        + nobreak("n", IN, Int(1, 20).pretty())
        + (f", {nobreak('j', IN, f'{{0, {ELLIPSIS}, n {MINUS} 1}}')}" if index else ""),
    )


INFS: Final = (-math.inf, math.inf)
add(
    "logsumexp",
    "x",
    r_logsumexp,
    xs_domain(U(-50, 50)),
    xs_domain(LU(1e-300, 1e300, 0)),
    xs_domain(U(-50, 50), INFS),
    cpp=False,
    note="arr",
)
for f, r in [("softmax", r_softmax), ("log_softmax", r_log_softmax)]:
    add(
        f,
        "j x",
        r,
        xs_domain(U(-50, 50), index=True),
        xs_domain(LU(1e-300, 1e300, 0), index=True),
        xs_domain(U(-50, 50), INFS[:1], index=True),
        cpp=False,
        note="arr elem",
    )


def poly_domain(z: Dist, c: Dist) -> Domain:
    """z, and the coefficients c of a polynomial of degree n in {0, ..., 20}."""
    return D(
        sampler=lambda r: (
            z.sample(r),
            z.sample(r),
            *(c.sample(r) for _ in range(r.randint(0, 20) + 1)),
        ),
        label=f"re, im in {z.desc()}, c in {c.desc()}^(n + 1), n in {{0, ..., 20}}",
        pretty=f"{nobreak(GREEK['re'] + ', ' + GREEK['im'], IN, z.pretty())}"
        + (" (log)" if _is_log(z) else "")
        + f", {array('c', c, length='\u207f\u207a\u00b9')}, "
        + nobreak("n", IN, Int(0, 20).pretty()),
    )


add(
    "cevalpoly",
    "re im c",
    c_cevalpoly,
    poly_domain(U(-2, 2), U(-1, 1)),
    poly_domain(LU(1e-3, 1e3, 0), LU(1e-3, 1e3, 0)),
    complex=True,
    note="zeros_c",
)


def ps_domain(dist: Dist) -> Domain:
    """k in {0, ..., n}, and the probabilities p of n in {1, ..., 50} trials."""

    def sampler(r: Random) -> list[float]:
        n = r.randint(1, 50)
        return [r.randint(0, n), *(dist.sample(r) for _ in range(n))]

    return D(
        sampler=sampler,
        label=f"k in {{0, ..., n}}, p in {dist.desc()}^n, n in {{1, ..., 50}}",
        pretty=f"{nobreak('k', IN, f'{{0, {ELLIPSIS}, n}}')}, {array('p', dist)}, "
        + nobreak("n", IN, Int(1, 50).pretty()),
        n=2000,
    )


for f, r in [
    ("poisson_binom_cdf", r_poisson_binom_cdf),
    ("poisson_binom_pmf_all", r_poisson_binom_pmf),
    ("poisson_binom_cdf_all", r_poisson_binom_cdf),
]:
    add(
        f,
        "k p",
        r,
        ps_domain(U(0, 1)),
        ps_domain(LU(1e-10, 0.1)),
        ps_domain(Offset(1, LU(1e-10, 0.1), -1)),
        cpp=f.endswith("_all"),
        note="elem" if f.endswith("_all") else "",
    )


def wilcoxon_domain(lo: int, hi: int, n: int = 10000) -> Domain:
    """n in {lo, ..., hi}, and k in {0, ..., n (n + 1) / 2}, the support."""
    return D(
        sampler=lambda r: (m := r.randint(lo, hi), r.randint(0, m * (m + 1) // 2)),
        label=f"n in {{{lo}, ..., {hi}}}, k in {{0, ..., n (n + 1) / 2}}",
        pretty=f"{nobreak('n', IN, Int(lo, hi).pretty())}, "
        + nobreak("k", IN, f"{{0, {ELLIPSIS}, n(n + 1)/2}}"),
        n=n,
    )


for f, r in [
    ("wilcoxon_pmf_all", r_wilcoxon_pmf),
    ("wilcoxon_cdf_all", r_wilcoxon_cdf),
    ("wilcoxon_pmf", r_wilcoxon_pmf),
    ("wilcoxon_cdf", r_wilcoxon_cdf),
    ("wilcoxon_sf", r_wilcoxon_sf),
]:
    add(
        f,
        "n k",
        r,
        wilcoxon_domain(0, 50),
        wilcoxon_domain(51, 300, 2000),
        note={"wilcoxon_sf": "abs"}.get(f, "elem" if f.endswith("_all") else ""),
    )

NOTES = {
    "zeros_c": "Ill-conditioned near the zeros.",
    "zeros": "Ill-conditioned near the zeros, where the largest errors are.",
    "osc": r"For $`x > \lvert v\rvert`$: $`A = \sqrt{J_v(x)^2 + Y_v(x)^2}`$.",
    "neg": (
        r"Ill-conditioned for negative non-integer $`a`$ with large $`\lvert a\rvert`$, and near the "
        r"poles of $`\Gamma(a + b)`$ (or $`\Gamma(a + m)`$), where the largest errors are."
    ),
    "osc_hyp0f1": (
        r"For $`z < 0`$ with $`x = 2\sqrt{-z} > \lvert b - 1\rvert`$: "
        r"$`A = \lvert\Gamma(b)\rvert (x/2)^{1-b} \sqrt{J_{b-1}(x)^2 + Y_{b-1}(x)^2}`$."
    ),
    "beta_dist": r"$`\mu = a / (a + b)`$, $`\sigma^2 = \mu(1 - \mu) / (a + b + 1)`$.",
    "ncx2": r"$`\mu = v + \mathrm{nc}`$, $`\sigma^2 = 2(v + 2\,\mathrm{nc})`$.",
    "root": (
        r"Reference: the root of $`F(x) = p`$, with $`F`$ the forward reference, which must be the "
        r"only sign change on a log-uniform grid of 225 points in $`[10^{-100}, 10^{12}]`$; the "
        "inputs with more than one are left out."
    ),
    "seq": r"The element of order $`v`$ (or $`n`$) of the sequence.",
    "none": (
        "inf (k) counts the `None` results, where it cannot compute the result to a relative "
        "error of about $`10^{-13}`$."
    ),
}

NOTES |= {
    "inv": (
        "Reference: the mpmath root of the forward reference. Ill-conditioned where the forward "
        "function is flat, e.g. in its tails."
    ),
    "osc_sph": (
        r"For $`x > n + 1/2`$: "
        r"$`A = \sqrt{\pi/(2x)}\,\sqrt{J_{n+1/2}(x)^2 + Y_{n+1/2}(x)^2}`$."
    ),
    "osc_c": r"$`A = \sqrt{\lvert J_v(z)\rvert^2 + \lvert Y_v(z)\rvert^2}`$.",
    "osc_cz": (
        r"For $`\lvert z\rvert > \lvert v\rvert`$: "
        r"$`A = s\sqrt{\lvert J_v(z)\rvert^2 + \lvert Y_v(z)\rvert^2}`$, with "
        r"$`s = e^{-\lvert\mathrm{Im}\,z\rvert}`$ for `bessel_je` and `bessel_ye`, else 1."
    ),
    "osc_hyp0f1_c": (
        r"For $`4\lvert z\rvert > (b - 1)^2`$, with $`x = 2\sqrt{-z}`$: "
        r"$`A = \lvert\Gamma(b)\rvert \lvert x/2\rvert^{1-b} "
        r"\sqrt{\lvert J_{b-1}(x)\rvert^2 + \lvert Y_{b-1}(x)\rvert^2}`$."
    ),
    "kelvin": (
        r"$`A`$ = the modulus of the complex function, "
        r"e.g. $`\lvert\mathrm{ber}\,x + i\,\mathrm{bei}\,x\rvert`$ for `ber` and `bei`."
    ),
    "airy": (
        r"For $`x < 0`$: $`A = \sqrt{\mathrm{Ai}(x)^2 + \mathrm{Bi}(x)^2}`$, "
        r"or $`\sqrt{\mathrm{Ai}'(x)^2 + \mathrm{Bi}'(x)^2}`$ for the derivatives."
    ),
    "binom_dist": r"$`\sigma^2 = np(1 - p)`$.",
    "tvn": "With `epsi` = $`10^{-14}`$, its lower bound.",
    "abs": (
        r"Absolute error $`\sim 10^{-16}`$ at best, so the relative error is large for tiny "
        "probabilities."
    ),
    "all": (
        "The last element: degree $`n`$ (and order $`m`$) of the output up to degree $`n`$ (and "
        r"order $`\lvert m\rvert`$)."
    ),
    "riccati": (
        r"Order $`n`$ of the output for the orders $`0, \ldots, 100`$ (`N = 101`). For "
        r"$`x > n + 1/2`$: "
        r"$`A = \sqrt{\pi x/2}\,\sqrt{J_{n+1/2}(x)^2 + Y_{n+1/2}(x)^2}`$."
    ),
    "lmbda": (
        r"The last order, $`v`$. For $`x > u`$, with $`u = v`$, or $`v + 1`$ for the derivative: "
        r"$`A = \Gamma(v + 1) (2/x)^v \sqrt{J_u(x)^2 + Y_u(x)^2}`$."
    ),
    "arr": (
        r"In the rows with $`\infty`$, each element is $`\pm\infty`$ with probability 1/4 (but "
        r"at least one is finite), or only $`-\infty`$ for `softmax` and `log_softmax`, which can "
        r"be NaN for $`+\infty`$, like in SciPy."
    ),
    "elem": r"Element $`j`$ or $`k`$ of the output.",
    "lame": (
        r"For $`s < k`$: $`A = \max\lvert E^p_n\rvert`$ over 129 points of $`[0, h]`$ or "
        r"$`[h, k]`$, whichever has $`s`$."
    ),
    "zeros_n": (
        r"Element $`k`$ of the outputs for `N = 500` (counting from 1): the zero $`k`$, or the "
        "value there."
    ),
    "zeros_k": r"The last element: zero $`k`$ of the first $`k`$ (`nt = k`).",
    "jnjnp": (
        r"Element $`k`$ of the outputs for `N` zeros: the zero $`x_k`$, its order $`n_k`$, its "
        r"serial number $`m_k`$, and whether it is a zero of $`J_{n_k}'`$ (1) or $`J_{n_k}`$ (0). "
        r"The zeros of $`J_1`$ and $`J_0' = -J_1`$ come in either order, by rounding; they are "
        r"compared as if $`J_1`$ is first."
    ),
    "y_zeros": (
        r"The real zeros (`complex = false`), or the complex ones with "
        r"$`\mathrm{Im}\, z > 0`$ (`complex = true`)."
    ),
    "gauss": (
        r"Node (`.0`) or weight (`.1`) $`k`$ of the $`n`$-point rule, with the nodes in ascending "
        "order."
    ),
    "numbers": r"Element $`n`$ of the output for `N = 500` (counting from 0).",
    "poly": (
        r"Rows with $`n \in \{\ldots\}`$ use an integer degree (`i32`, or `u32` for Hermite), "
        "the others an `f64` degree."
    ),
    "mathieu": (
        "Reference of this section: the eigenvalues and eigenvectors of the recurrence relations "
        "[DLMF 28.4.5-8](https://dlmf.nist.gov/28.4#ii), the Fourier series "
        "[28.4.1-4](https://dlmf.nist.gov/28.4#i), and the series of Bessel function products "
        "[28.24.2-5](https://dlmf.nist.gov/28.24), in mpmath."
    ),
    "mathieu_cv": (
        "For some $`m`$ and $`q`$, the result is the characteristic value of another order of "
        "the same parity, which the other functions of this section then use. For $`m = 0`$ and "
        r"$`\lvert q\rvert \lesssim 10^{-162}`$ it is NaN."
    ),
    "mathieu_osc": (
        r"$`A = \sqrt{y^2 + y'^2 / w}`$ for $`y = \mathrm{ce}_m`$ or $`\mathrm{se}_m`$, and "
        r"$`\sqrt{w y^2 + y'^2}`$ for $`y'`$, with "
        r"$`w = \max(\lvert a - 2 q \cos 2x\rvert, 1)`$ and $`a`$ the characteristic value."
    ),
    "mathieu_mod": (
        r"For $`2 q \cosh 2x > a`$, with $`a`$ the characteristic value: "
        r"$`A = \sqrt{M^{(1)}(x)^2 + M^{(2)}(x)^2}`$ for the first and second kind $`M^{(j)}`$ "
        r"(or their derivatives). For $`q \le 10^{-7}`$, the results are infinite for "
        r"$`m \ge 4`$ (`modcem`) or $`m \ge 5`$ (`modsem`)."
    ),
    "mathieu_coef": (
        r"Coefficient $`k`$ of the output. For $`q \le 10^{-7}`$, only the largest coefficient "
        "and the two next to it are nonzero."
    ),
    "sw_seq": r"Mode $`n`$ of the output for the modes $`m, \ldots, 20`$.",
    "sw_cv": "With the characteristic value `cv` rounded to the nearest f64.",
    "sw_ang": (
        r"For $`q = \lambda \mp c^2 x^2 - m^2 / (1 - x^2) > 0`$ (prolate: $`-`$, oblate: $`+`$), "
        r"with the characteristic value $`\lambda`$: "
        r"$`A = \sqrt{S^2 + (1 - x^2) S'^2 / \hat q}`$ for the angular function $`S`$, and "
        r"$`\sqrt{\hat q / (1 - x^2)}\,A`$ for $`S'`$, with "
        r"$`\hat q = \max(q, ((1 - x^2) q'^2)^{1/3})`$."
    ),
    "sw_rad": (
        r"For $`c^2 x^2 > \lambda \pm m^2 / (x^2 \mp 1)`$ (prolate: upper signs, oblate: lower): "
        r"$`A = \sqrt{R_1^2 + R_2^2}`$ for the radial functions of the first and second kind "
        r"$`R_1`$ and $`R_2`$, and $`\sqrt{R_1'^2 + R_2'^2}`$ for their derivatives."
    ),
}

NOT_COVERED = """\
## Not covered

- Exact functions: `comb`, `comb_rep`, `factorial`, `multifactorial` (and `_checked`) on unsigned
  integers, where an error is a wrong integer, not a rounding error, and `round`, which is exact
- Complex input: `sph_bessel_*`, `elliprc`, `elliprd`, `elliprf`, `elliprg`, `elliprj`, `sici`,
  `shichi`, `spence`, `wrightomega`, `*legendre_p*`, `*legendre_q_all`, `eval_*`, and `zeta`
  (NaN for non-real `z`, unless `q = 1`)
- Other: `cdf_cvm` (an approximation itself), `lanczos_sum_expg_scaled` (internal), the
  `extended_*_error` metrics, `kelvin`, and `kelvin_zeros` (their parts are above)
"""

BY_NAME: Final = {s.name: s for s in SPECS}

# the same functions in Arb (python-flint), on balls, to check the references (`--arb-check`)
ARB: Final[dict[str, Callable[..., flint.arb | flint.acb]]] = {
    "gamma": lambda x: x.gamma(),
    "rgamma": lambda x: x.rgamma(),
    "digamma": lambda x: x.digamma(),
    "polygamma": lambda n, x: flint.acb(x).polygamma(n),
    "erfinv": lambda y: y.erfinv(),
    "erfcinv": lambda y: y.erfcinv(),
    "bessel_j0": lambda x: x.bessel_j(0),
    "bessel_j1": lambda x: x.bessel_j(1),
    "bessel_y0": lambda x: x.bessel_y(0),
    "bessel_y1": lambda x: x.bessel_y(1),
    "bessel_i0": lambda x: x.bessel_i(0),
    "bessel_i1": lambda x: x.bessel_i(1),
    "bessel_k0": lambda x: x.bessel_k(0),
    "bessel_k1": lambda x: x.bessel_k(1),
    "bessel_j": lambda v, x: x.bessel_j(v),
    "bessel_y": lambda v, x: x.bessel_y(v),
    "bessel_i": lambda v, x: x.bessel_i(v),
    "bessel_k": lambda v, x: x.bessel_k(v),
    "zeta": lambda x, q: x.zeta(q),
    "riemann_zeta": lambda x: x.zeta(),
    "elliprc": flint.acb.elliptic_rc,
    "elliprd": flint.acb.elliptic_rd,
    "elliprf": flint.acb.elliptic_rf,
    "elliprg": flint.acb.elliptic_rg,
    "elliprj": flint.acb.elliptic_rj,
    "hyp0f1": lambda b, z: z.hypgeom_0f1(b),
    "hyp1f1": lambda a, b, x: x.hypgeom_1f1(a, b),
    "hypu": lambda a, b, x: x.hypgeom_u(a, b),
    "hyp2f1": lambda a, b, c, x: x.hypgeom_2f1(a, b, c),
    "hyp3f2": lambda a1, a2, a3, b1, b2, z: z.hypgeom([a1, a2, a3], [b1, b2]),
    "hyp2f2": lambda a1, a2, b1, b2, z: z.hypgeom([a1, a2], [b1, b2]),
    "hyp1f2": lambda a, b1, b2, z: z.hypgeom([a], [b1, b2]),
    # the Borel sum of the divergent series, for z < 0 (DLMF 13.6.21)
    "hyppfq#2f0": lambda a1, a2, z: (
        (-1 / z) ** a1 * (-1 / z).hypgeom_u(a1, a1 - a2 + 1)
    ),
    **{
        f"hyppfq#{p}f{q}": lambda *v, p=p: v[-1].hypgeom(list(v[:p]), list(v[p:-1]))
        for p, q in [(3, 3), (2, 3), (1, 3), (0, 2), (0, 3)]
    },
}
