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
a direct Euler-Maclaurin sum for the Hurwitz zeta function, quadrature for `expn`, the power
series for `gammainc` with large `a`, and a continued fraction for `betainc`. The references of
inverse functions find the root of the reference of the function itself (see `_invert`).
"""

import math
from dataclasses import dataclass
from fractions import Fraction

import mpmath as mp

DBL_MAX = 1.7976931348623157e308

# ---------------------------------------------------------------- sampling

# Unicode math for the table, escaped to keep the source ASCII
IN, MINUS, PM, ELLIPSIS, SQRT = "\u2208", "\u2212", "\u00b1", "\u2026", "\u221a"
HALF, MU, SIGMA, GAMMA, EPSILON, PSI = "\u00bd", "\u03bc", "\u03c3", "\u0393", "\u03b5", "\u03c8"
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
}


def nobreak(*parts):
    """Join the parts of one argument's domain with no-break spaces, and keep a leading sign
    attached to its interval, so that table cells only wrap between the arguments, or after
    a comma within the interval (the last part)."""
    *head, tail = (p.replace(" ", "\u00a0") for p in parts)
    s = "\u00a0".join([*head, tail.replace(",\u00a0", ", ")])
    return s.replace(f"{PM}[", f"{PM}\u2060[").replace(f"{MINUS}[", f"{MINUS}\u2060[")


def wrap(s):
    """Put each argument's domain on its own line, i.e. break at the ", " between them (the ones
    outside brackets; within an argument's domain, `nobreak` leaves only those in intervals)."""
    parts, depth, start = [], 0, 0
    for i, c in enumerate(s):
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == "," and depth == 0 and s[i + 1 : i + 2] == " ":
            parts.append(s[start:i])
            start = i + 2
    return "<br>".join([*parts, s[start:]])


@dataclass(frozen=True)
class U:
    """Uniform on [lo, hi]."""

    lo: float
    hi: float

    def sample(self, rng):
        return rng.uniform(self.lo, self.hi)

    def desc(self):
        return f"[{_g(self.lo)}, {_g(self.hi)}]"

    def pretty(self):
        return f"[{_u(self.lo)}, {_u(self.hi)}]"


@dataclass(frozen=True)
class LU:
    """Log-uniform on [lo, hi] (0 < lo < hi), optionally with a random sign."""

    lo: float
    hi: float
    sign: int = 1  # 1, -1, or 0 for a random sign

    def sample(self, rng):
        x = 10 ** rng.uniform(math.log10(self.lo), math.log10(self.hi))
        s = self.sign or rng.choice((-1, 1))
        return s * min(x, self.hi)

    def desc(self):
        pre = {1: "", -1: "-", 0: "+-"}[self.sign]
        return f"{pre}[{_g(self.lo)}, {_g(self.hi)}] (log)"

    def pretty(self):
        pre = {1: "", -1: MINUS, 0: PM}[self.sign]
        return f"{pre}[{_u(self.lo)}, {_u(self.hi)}]"


@dataclass(frozen=True)
class Int:
    """Uniform integer in [lo, hi]."""

    lo: int
    hi: int

    def sample(self, rng):
        return float(rng.randint(self.lo, self.hi))

    def desc(self):
        return f"{{{self.lo}, ..., {self.hi}}}"

    def pretty(self):
        return f"{{{self.lo}, {ELLIPSIS}, {self.hi}}}"


@dataclass(frozen=True)
class ILU:
    """Log-uniform integer in [lo, hi] (1 <= lo < hi)."""

    lo: int
    hi: int

    def sample(self, rng):
        return float(min(int(LU(self.lo, self.hi + 1).sample(rng)), self.hi))

    def desc(self):
        return f"{{{_g(self.lo)}, ..., {_g(self.hi)}}} (log)"

    def pretty(self):
        return f"{{{_u(self.lo)}, {ELLIPSIS}, {_u(self.hi)}}}"


@dataclass(frozen=True)
class Offset:
    """c + sign * t (rounded to f64), where t is sampled from `dist`."""

    c: float
    dist: object
    sign: int = 1

    def sample(self, rng):
        return self.c + self.sign * self.dist.sample(rng)

    def desc(self):
        return f"{_g(self.c)} {'+' if self.sign > 0 else '-'} {self.dist.desc()}"

    def pretty(self):
        return f"{_u(self.c)} {'+' if self.sign > 0 else MINUS} {self.dist.pretty()}"


def _g(x):
    """Short repr of a domain bound."""
    if x == int(x) and abs(x) < 1e6:
        return str(int(x))
    s = f"{x:.4g}"
    return s.replace("e+", "e").replace("e-0", "e-").replace("e0", "e")


def _u(x):
    """Unicode repr of a domain bound, e.g. 10^-3 with a superscript exponent."""
    sign = MINUS if x < 0 else ""
    a = abs(x)
    named = {math.pi: "\u03c0", math.pi / 2: "\u03c0/2", -INV_E: "1/e"}
    if a in named:
        return sign + named[a]
    e = round(math.log10(a)) if a else 0
    if (e <= -3 or e >= 4) and a == float(f"1e{e}"):
        return f"{sign}10{str(e).translate(SUP)}"
    return sign + _g(a)


def _is_log(d):
    return isinstance(d, (LU, ILU)) or isinstance(d, Offset) and _is_log(d.dist)


class Domain:
    """Independent per-argument distributions, or a custom `sampler(rng) -> args`.

    `call` overrides the function name passed to the dump programs (e.g. for an integer-degree
    variant), and `cpp=False` marks a variant that xsf does not have.
    """

    def __init__(self, *dists, n=10000, sampler=None, label="", pretty="", call="", cpp=True):
        self.dists, self.n, self.sampler = dists, n, sampler
        self.label, self.pretty_label = label, pretty
        self.call, self.cpp = call, cpp

    def desc(self, argnames):
        """Plain description; also seeds the sampling."""
        if self.label:
            return self.label
        return ", ".join(f"{a} in {d.desc()}" for a, d in zip(argnames, self.dists))

    def pretty(self, argnames):
        """Unicode description for the table; consecutive equal distributions are merged."""
        if self.pretty_label:
            return wrap(self.pretty_label)
        groups = []
        for a, d in zip(argnames, self.dists):
            if groups and groups[-1][1] == d:
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

    def sample(self, rng):
        if self.sampler is not None:
            return tuple(float(a) for a in self.sampler(rng))
        return tuple(d.sample(rng) for d in self.dists)


D = Domain


@dataclass
class Spec:
    name: str  # a suffix "#..." distinguishes specs of the same function, e.g. for complex input
    args: tuple
    ref: object
    domains: list
    section: str
    std: bool = False  # has a Rust std equivalent
    cpp: bool = True  # has an xsf (C++) equivalent
    complex: bool = False  # complex-valued
    dps: int = 50
    timeout: int = 30
    note: str = ""  # space-separated keys into NOTES
    rust = True

    @property
    def func(self):
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
    """P(a, x) via the power series (DLMF 8.11.4); mpmath does not converge for a > 1e5."""
    with wdps(xd(a, x) * 2):
        t = s = mp.mpf(1)
        k = 0
        while t > mp.eps * s:
            k += 1
            t *= x / (a + k)
            s += t
        return mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1)) * s


def r_gammainc(a, x):
    if a > 1e5:
        return _gammainc_series(a, x)
    with wdps(xd(a, x) * 2):
        return mp.gammainc(a, 0, x, regularized=True)


def r_gammaincc(a, x):
    with wdps(xd(a, x) * 2):
        if a > 1e5:
            return 1 - _gammainc_series(a, x)
        return mp.gammainc(a, x, mp.inf, regularized=True)


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


def r_hyp0f1(b, z):
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
        return d * (d - 1) / 4 * mp.log(mp.pi) + mp.fsum(mp.loggamma(a - mp.mpf(j) / 2) for j in range(d))


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
        if prev is not None and abs(p - prev) <= abs(p) * mp.mpf(10) ** (-mp.mp.dps - 5):
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
        s = mp.fsum(c**k * mp.binomial(n, k) * f(v - n + 2 * k, x) for k in range(n + 1))
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
        return a**v / ((lmbda + v + 1) * mp.gamma(v + 1)) * mp.hyp1f2(c, v + 1, c + 1, -a * a)


def r_itstruve0(x):
    # int_0^x H_0(t) dt = x^2/pi 2F3(1, 1; 3/2, 3/2, 2; -x^2/4)
    with wdps(xd(x) + 10):
        return x * x / mp.pi * mp.hyper([1, 1], [1.5, 1.5, 2], -x * x / 4)


def r_it2struve0(x):
    # int_x^inf H_0(t) / t dt = pi/2 - 2x/pi 2F3(1/2, 1; 3/2, 3/2, 3/2; -x^2/4)
    with wdps(xd(x) + 10):
        return mp.pi / 2 - 2 * x / mp.pi * mp.hyper([0.5, 1], [1.5, 1.5, 1.5], -x * x / 4)


def r_itmodstruve0(x):
    # int_0^x L_0(t) dt = x^2/pi 2F3(1, 1; 3/2, 3/2, 2; x^2/4)
    with wdps(xd(x) + 10):
        return x * x / mp.pi * mp.hyper([1, 1], [1.5, 1.5, 2], x * x / 4)


# ---------------------------------------------------------------- Kelvin and Airy functions


def _kelvin(kind, x):
    """ber, bei, ker, kei (kind 0 to 3), or their derivatives (kind 4 to 7), with the modulus of
    the complex function (e.g. |ber + i bei|) as scale. ber and bei are even."""
    s2 = mp.sqrt(2)
    # the sums of the order-1 functions cancel for small x, by up to a factor x^2
    with wdps(xd(x) + 2 * max(0, int(-mp.log10(abs(x))))):
        a = abs(x)
        if kind < 2:
            f = mp.ber(0, a), mp.bei(0, a)
        elif kind < 4:
            f = mp.ker(0, a), mp.kei(0, a)
        elif kind < 6:  # ber'(x) = (ber_1(x) + bei_1(x)) / sqrt(2), DLMF 10.63.1
            b1, c1 = mp.ber(1, a), mp.bei(1, a)
            f = mp.sign(x) * (b1 + c1) / s2, mp.sign(x) * (c1 - b1) / s2
        else:
            b1, c1 = mp.ker(1, a), mp.kei(1, a)
            f = (b1 + c1) / s2, (c1 - b1) / s2
        return f[kind % 2], mp.hypot(*f)


def _airy(x, scaled=False):
    """Ai, Ai', Bi, Bi', with the scale sqrt(Ai^2 + Bi^2) (or of the derivatives) for x < 0."""
    with wdps(xd(x) + 10):
        ai, aip = mp.airyai(x), mp.airyai(x, derivative=1)
        bi, bip = mp.airybi(x), mp.airybi(x, derivative=1)
        if scaled and x > 0:
            z = 2 * x * mp.sqrt(x) / 3
            ai, aip, bi, bip = ai * mp.exp(z), aip * mp.exp(z), bi * mp.exp(-z), bip * mp.exp(-z)
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
            return c1 * t * mp.hyp1f2(third, 2 * third, 4 * third, z), c2 * t * t / 2 * mp.hyp1f2(
                2 * third, 4 * third, 5 * third, z
            )

        f, g = fg(x)
        fn, gn = fg(-x)
        ai = mp.airyai(x, derivative=-1) if x <= 100 else third - mp.airyai(x) / mp.sqrt(x)
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


def r_pbwa(a, x):
    with wdps(20):
        return mp.pcfw(a, x), mp.diff(lambda t: mp.pcfw(a, t), x)


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


def _gammaincc_cf(a, x):
    """Q(a, x) for x > a, via the modified Lentz continued fraction for Gamma(a, x)."""
    tiny = mp.mpf(10) ** (-2 * mp.mp.dps)
    b = x + 1 - a
    c, d = 1 / tiny, 1 / b
    h, i = d, 0
    while True:
        i += 1
        an, b = -i * (i - a), b + 2
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1 / d
        h *= d * c
        if abs(d * c - 1) < mp.eps:
            return mp.exp(-x + a * mp.log(x) - mp.loggamma(a)) * h


def _gammainc_pq(a, x, upper):
    """P(a, x) or Q(a, x), also where mpmath does not converge."""
    try:
        return r_gammaincc(a, x) if upper else r_gammainc(a, x)
    except (mp.libmp.NoConvergence, ValueError):
        with wdps(xd(a, x) * 2):
            q = _gammaincc_cf(a, x)
            return q if upper else 1 - q


def r_gammaincinv(a, p):
    return _invert_pos(lambda x: _gammainc_pq(a, x, False), p, a)


def r_gammainccinv(a, q):
    return _invert_pos(lambda x: _gammainc_pq(a, x, True), q, a, decreasing=True)


def r_betaincinv(a, b, y):
    return _invert_unit(lambda x, z: _betainc(a, b, x, z), y, a / (a + b))


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
        if prev is not None and mp.isfinite(y) and abs(y - prev) <= abs(y) * mp.mpf(10) ** (-mp.mp.dps - 5):
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
        sn, cn, dn = (mp.ellipfun(k, u, m=m) for k in ("sn", "cn", "dn"))
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


def r_hyp2f1(a, b, c, x):
    with wdps(20):
        return mp.hyp2f1(a, b, c, x, maxterms=10**6)


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
    return _invert_unit(lambda p, q: _betainc(n - k, k + 1, q, p), y, (k + 1) / (n + 1), True)


def r_nbdtr(k, n, p):
    return _betainc(n, k + 1, p, 1 - p)


def r_nbdtrc(k, n, p):
    return _betainc(k + 1, n, 1 - p, p)


def r_nbdtri(k, n, y):
    # p such that nbdtr(k, n, p) = y
    return _invert_unit(lambda p, q: _betainc(n, k + 1, p, q), y, n / (n + k + 1))


def r_pdtr(k, m):
    return _gammainc_pq(k + 1, m, True)


def r_pdtrc(k, m):
    return _gammainc_pq(k + 1, m, False)


def r_pdtri(k, y):
    # m such that pdtr(k, m) = y
    return _invert_pos(lambda m: _gammainc_pq(k + 1, m, True), y, k + 1, decreasing=True)


def r_gdtr(a, b, x):
    return _gammainc_pq(b, a * x, False)


def r_gdtrc(a, b, x):
    return _gammainc_pq(b, a * x, True)


def r_gdtria(p, b, x):
    # a such that gdtr(a, b, x) = p
    return _invert_pos(lambda a: _gammainc_pq(b, a * x, False), p, b / x)


def r_gdtrib(a, p, x):
    # b such that gdtr(a, b, x) = p
    return _invert_pos(lambda b: _gammainc_pq(b, a * x, False), p, a * x, decreasing=True)


def r_gdtrix(a, b, p):
    # x such that gdtr(a, b, x) = p
    return _invert_pos(lambda x: _gammainc_pq(b, a * x, False), p, b / a)


def r_chdtr(v, x):
    return _gammainc_pq(v / 2, x / 2, False)


def r_chdtrc(v, x):
    return _gammainc_pq(v / 2, x / 2, True)


def r_chdtri(v, y):
    # x such that chdtrc(v, x) = y
    return _invert_pos(lambda x: _gammainc_pq(v / 2, x / 2, True), y, v, decreasing=True)


def r_fdtri(a, b, y):
    # x such that fdtr(a, b, x) = y
    return _invert_pos(lambda x: r_fdtr(a, b, x), y, 1)


def _stdtr_tail(v, t):
    """P(T <= -|t|) = I_x(v/2, 1/2) / 2, with x = v / (v + t^2)."""
    return _betainc(v / 2, mp.mpf(1) / 2, v / (v + t * t), t * t / (v + t * t)) / 2


def r_stdtr(v, t):
    with wdps(xd(v, t)):
        tail = _stdtr_tail(v, t)
        return tail if t < 0 else 1 - tail


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
        f = lambda t: mp.exp(k * (mp.cos(t) - c))  # noqa: E731
        w = 1 / mp.sqrt(k)
        pts = [-mp.pi] + [t for t in (-10 * w, -w, 0, w, 10 * w) if -mp.pi < t < x] + [x]
        return mp.quad(f, pts) * mp.exp(k * (c - 1)) / (2 * mp.pi * mp.besseli(0, k) * mp.exp(-k))


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
            t = mp.gamma(k + 0.5) / (mp.gamma(0.5) * mp.factorial(k)) * mp.sqrt(4 * k + 1)
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
        mp.binomial(n, j) * (x + mp.mpf(j) / n) ** (j - 1) * (mp.mpf(n - j) / n - x) ** (n - j)
        for j in js
    ]


def _smirnov_jmax(n, x):
    """floor(n (1 - x)), exactly."""
    return math.floor(n * (1 - Fraction(float(x))))


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
            s = mp.nsum(lambda k: (-1) ** k * 4 * k * k * x * mp.exp(-2 * k * k * x * x), [1, mp.inf])
        return 2 * s, 1 - 2 * s if not deriv else -2 * s
    c = mp.sqrt(2 * mp.pi)

    def a(k):
        return (2 * k - 1) ** 2 * mp.pi**2 / 8

    s = c / x * mp.nsum(lambda k: mp.exp(-a(k) / x**2), [1, mp.inf])
    if deriv:
        s = mp.nsum(lambda k: mp.exp(-a(k) / x**2) * (2 * a(k) / x**4 - 1 / x**2), [1, mp.inf]) * c
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
    s, t, k, peaked = mp.mpf(0), mp.mpf(1), 0, False
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
    terms = [(-1) ** k * mp.binomial(n, k) * f(v - n + 2 * k, mp.mpc(re, im)) for k in range(n + 1)]
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


# ---------------------------------------------------------------- specs

SPECS = []
SECTIONS = []  # section titles, in order of appearance
_section = None


def section(title):
    global _section
    _section = title
    SECTIONS.append(title)


def add(name, args, ref, *domains, **kw):
    """Add a spec for the function `name` with the space-separated argument names `args`."""
    SPECS.append(Spec(name, tuple(args.split()), ref, list(domains), _section, **kw))


INV_E = float(-1 / mp.e)  # the smallest f64 > -1/e
if INV_E < -1 / mp.e:
    INV_E = math.nextafter(INV_E, 0)


def a_pm_sqrt_a(r):
    """(a, x) with a in [10, 1e8] (log), and x in a +- 6 sqrt(a): the transition region."""
    a = 10 ** r.uniform(1, 8)
    return a, max(0.0, a + r.uniform(-6, 6) * a**0.5)


A_PM_SQRT_A = D(
    sampler=a_pm_sqrt_a,
    label="a in [10, 1e8] (log), x in a +- 6 sqrt(a)",
    pretty=f"{nobreak('a', IN, LU(10, 1e8).pretty())} (log), "
    + nobreak("x", IN, f"a {PM} 6{SQRT}a"),
)


def ab_mu_pm_sigma(r):
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


def v_near_half_int(r):
    """(v, x) with v in -k - 1/2 +- [1e-12, 1e-3] (log) for k in {0, ..., 20}, and x in [0, 100]."""
    t = 10 ** r.uniform(-12, -3)
    return -r.randint(0, 20) - 0.5 + r.choice((-1, 1)) * t, r.uniform(0, 100)


V_NEAR_HALF_INT = D(
    sampler=v_near_half_int,
    label="v in -k - 1/2 +- [1e-12, 1e-3] (log), k in {0, ..., 20}, x in [0, 100]",
    pretty=nobreak("v", IN, f"{MINUS}k {MINUS} {HALF} {PM} {LU(1e-12, 1e-3).pretty()}")
    + f" (log), {nobreak('k', IN, Int(0, 20).pretty())}, {nobreak('x', IN, U(0, 100).pretty())}",
)

section("Gamma and related functions")
add("gamma", "x", r_gamma, D(U(-30, 30)), D(U(-170, 171.6)), D(LU(1e-300, 171)))
add("gammaln", "x", r_gammaln, D(U(-30, 30)), D(LU(1e-300, 1e300)), note="zeros")
add("rgamma", "x", r_rgamma, D(U(-30, 30)), D(U(-200, 200)))
add("digamma", "x", r_digamma, D(U(-30, 30)), D(LU(1e-300, 1e300)), note="zeros")
add("digammainv", "y", r_digammainv, D(U(-10, 10)), D(U(-1000, 700)))
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
        sampler=lambda r: ((d := r.randint(1, 20)) / 2 - 0.5 + 10 ** r.uniform(-3, 3), d),
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
add("perm", "n k", r_ff, D(U(0, 100), Int(0, 30)), cpp=False)

section("Error function and normal distribution")
add("erf", "x", r_erf, D(U(-6, 6)), D(LU(1e-300, 10, 0)))
add("erfc", "x", r_erfc, D(U(-6, 27)), D(LU(1e-300, 27, 0)))
add("erfcx", "x", r_erfcx, D(U(-26, 26)), D(LU(1e-300, 1e300)))
add("erfinv", "y", r_erfinv, D(U(-1, 1)), D(Offset(1, LU(1e-16, 1), -1)))
add("erfcinv", "y", r_erfcinv, D(U(0, 2)), D(LU(1e-300, 1)))
add("ndtr", "x", r_ndtr, D(U(-10, 10)), D(LU(1e-300, 38, 0)))
add("ndtri", "p", r_ndtri, D(U(0, 1)), D(LU(1e-300, 0.5)))
add("log_ndtr", "x", r_log_ndtr, D(U(-20, 6)), D(LU(1e-300, 1e10, 0)))

add("erfi", "x", r_erfi, D(U(-6, 6)), D(LU(1e-300, 26.6, 0)))
add("dawsn", "x", r_dawsn, D(U(-10, 10)), D(LU(1e-300, 1e300, 0)))
add(
    "voigt_profile",
    "x sigma gamma",
    r_voigt_profile,
    D(U(-10, 10), U(0, 5), U(0, 5)),
    D(LU(1e-3, 1e3, 0), LU(1e-3, 1e3), LU(1e-3, 1e3)),
)
add("ndtri_exp", "y", r_ndtri_exp, D(U(-20, 0)), D(LU(1e-300, 1e300, -1)))
add("owens_t", "h a", r_owens_t, D(U(-5, 5), U(-5, 5)), D(LU(1e-5, 40, 0), LU(1e-5, 1e3, 0)))
add(
    "bivariate_normal_cdf",
    "h k r",
    r_bivariate_normal_cdf,
    D(U(-10, 10), U(-10, 10), U(-1, 1), n=2000),
    timeout=60,
)


def tvn_sampler(r):
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
        + nobreak(
            "r\u2082\u2083",
            IN,
            f"r\u2081\u2082r\u2081\u2083 {PM} 0.9{SQRT}((1 {MINUS} r\u2081\u2082\u00b2)(1 {MINUS} r\u2081\u2083\u00b2))",
        ),
        n=1000,
    ),
    timeout=120,
    note="tvn abs",
)
add(
    "nrdtrimn",
    "p std x",
    r_nrdtrimn,
    D(U(0, 1), LU(1e-3, 1e3), U(-10, 10)),
    note="zeros_c",
)
add("nrdtrisd", "mean p x", r_nrdtrisd, D(U(-10, 10), U(0, 1), U(-10, 10)))

section("Elementary functions")
add("expm1", "x", r_expm1, D(U(-2, 2)), D(LU(1e-300, 709, 0)), std=True)
add("log1p", "x", r_log1p, D(U(-1, 1)), D(LU(1e-300, 1e300)), D(LU(1e-300, 1, -1)), std=True)
add("cbrt", "x", r_cbrt, D(LU(1e-300, 1e300, 0)), std=True)
add("exp2", "x", r_exp2, D(U(-1, 1)), D(U(-1022, 1023)), std=True)
add("exp10", "x", r_exp10, D(U(-1, 1)), D(U(-307, 308)))
add("sinpi", "x", r_sinpi, D(U(-10, 10)), D(LU(1e-300, 1e15, 0)))
add("cospi", "x", r_cospi, D(U(-10, 10)), D(LU(1e-300, 1e15, 0)))
add("cosine_cdf", "x", r_cosine_cdf, D(U(-math.pi, math.pi)))
add("cosine_invcdf", "p", r_cosine_invcdf, D(U(0, 1)), D(LU(1e-300, 0.5)))

add("cosm1", "x", r_cosm1, D(U(-math.pi, math.pi)), D(LU(1e-300, 1e3, 0)), note="zeros_c")
for f, r in [("sindg", r_sindg), ("cosdg", r_cosdg), ("tandg", r_tandg), ("cotdg", r_cotdg)]:
    add(f, "x", r, D(U(-360, 360)), D(LU(1e-300, 1e14, 0)))
add("radian", "d m s", r_radian, D(U(0, 360), U(0, 60), U(0, 60)))
add("sinc", "x", r_sinc, D(U(-10, 10)), D(LU(1e-300, 1e15, 0)), cpp=False, note="zeros_c")
add("diric", "n x", r_diric, D(Int(1, 20), U(-10, 10)), cpp=False, note="zeros_c")
add("agm", "a b", r_agm, D(U(0, 10), U(0, 10)), D(LU(1e-300, 1e300), LU(1e-300, 1e300)))

section("Exponential and logarithmic functions")
add("expit", "x", r_expit, D(U(-10, 10)), D(LU(1e-300, 750, 0)))
add("logit", "p", r_logit, D(U(0, 1)), D(LU(1e-300, 0.5)), note="zeros_c")
add("log_expit", "x", r_log_expit, D(U(-20, 20)), D(LU(1e-300, 1e300, 0)))
add("exprel", "x", r_exprel, D(U(-2, 2)), D(LU(1e-300, 709, 0)))
add("log1mexp", "x", r_log1mexp, D(U(-20, 0)), D(LU(1e-300, 1e300, -1)))
add("log1pmx", "x", r_log1pmx, D(U(-1, 2)), D(LU(1e-300, 1, 0)), D(LU(1, 1e300)))
add(
    "xlogy",
    "x y",
    r_xlogy,
    D(LU(1e-300, 1e300, 0), LU(1e-300, 1e300)),
    D(U(-10, 10), Offset(1, LU(1e-16, 1e-3, 0))),
)
add("xlog1py", "x y", r_xlog1py, D(LU(1e-300, 1e300, 0), LU(1e-300, 1e300)), D(U(-10, 10), U(-1, 1)))
add("softplus", "x", r_softplus, D(U(-20, 20)), D(LU(1e-300, 1e300, 0)), cpp=False)
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

section("Box-Cox transformations and convex analysis")
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
    add(f, "y lambda", r, D(U(-1, 1), U(-1, 1)), D(LU(1e-300, 10, 0), LU(1e-300, 0.1, 0)))
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

section("Bessel and Struve functions")
for f, r in [("j0", r_j0), ("j1", r_j1), ("y0", r_y0), ("y1", r_y1)]:
    add(f"bessel_{f}", "x", r, D(U(0, 30)), D(LU(1e-300, 1e8)), note="osc")
for f, r in [("i0", r_i0), ("i1", r_i1)]:
    add(f"bessel_{f}", "x", r, D(U(-30, 30)), D(LU(1e-300, 700)))
for f, r in [("k0", r_k0), ("k1", r_k1)]:
    add(f"bessel_{f}", "x", r, D(U(0, 30)), D(LU(1e-300, 700)))
add("bessel_j", "v x", r_jv, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 1e4)), note="osc")
add("bessel_y", "v x", r_yv, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 1e4)), note="osc")
add("bessel_i", "v x", r_iv, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 700)))
add("bessel_k", "v x", r_kv, D(U(0, 20), U(0, 50)), D(U(0, 100), LU(1e-300, 700)))
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

section("Scaled Bessel functions, derivatives, and ratios")
for f, r in [("i0e", r_i0e), ("i1e", r_i1e)]:
    add(f"bessel_{f}", "x", r, D(U(-30, 30)), D(LU(1e-300, 1e300, 0)))
for f, r in [("k0e", r_k0e), ("k1e", r_k1e)]:
    add(f"bessel_{f}", "x", r, D(U(0, 30)), D(LU(1e-300, 1e300)))
add("bessel_je", "v x", r_jv, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 1e8)), note="osc")
add("bessel_ye", "v x", r_yv, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 1e8)), note="osc")
add("bessel_ie", "v x", r_ie, D(U(0, 20), U(0, 50)), D(U(-100, 100), LU(1e-300, 1e8)))
add("bessel_ke", "v x", r_ke, D(U(0, 20), U(0, 50)), D(U(0, 100), LU(1e-300, 1e8)))
for f, r in [("j", r_jvp), ("y", r_yvp), ("i", r_ivp), ("k", r_kvp)]:
    add(
        f"bessel_{f}_prime",
        "v x n",
        r,
        D(U(0, 20), U(0, 50), Int(1, 3)),
        cpp=False,
        note="osc" if f in "jy" else "",
    )
add("iv_ratio", "v x", r_iv_ratio, D(U(0.5, 20), U(0, 50)), D(LU(0.5, 1e4), LU(1e-300, 1e4)))
add("iv_ratio_c", "v x", r_iv_ratio_c, D(U(0.5, 20), U(0, 50)), D(LU(0.5, 1e4), LU(1e-300, 1e4)))
add("iv_ratioinv", "v r", r_iv_ratioinv, D(U(0.5, 20), U(0, 1), n=2000), note="inv")

section("Spherical Bessel functions")
for f, r, rp in [("j", r_sph_j, r_sph_jp), ("y", r_sph_y, r_sph_yp)]:
    for g, ref in [(f"sph_bessel_{f}", r), (f"sph_bessel_{f}_prime", rp)]:
        add(g, "n x", ref, D(Int(0, 50), U(0, 50)), D(Int(0, 100), LU(1e-300, 1e4)), note="osc_sph")
for f, r, rp in [("i", r_sph_i, r_sph_ip), ("k", r_sph_k, r_sph_kp)]:
    for g, ref in [(f"sph_bessel_{f}", r), (f"sph_bessel_{f}_prime", rp)]:
        add(g, "n x", ref, D(Int(0, 50), U(0, 50)), D(Int(0, 100), LU(1e-300, 700)))

section("Integrals of Bessel and Struve functions")
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
    add(f, "x", r, D(U(0, 30)), D(LU(1e-300, hi)), note="zeros_c" if f[-1] == "1" else "")
add("besselpoly", "a lambda v", r_besselpoly, D(U(0, 10), U(0, 10), U(0, 10)), note="zeros_c")
add("itstruve0", "x", r_itstruve0, D(U(0, 30)), D(LU(1e-300, 1e4)))
add("it2struve0", "x", r_it2struve0, D(U(0, 30)), D(LU(1e-300, 1e4)))
add("itmodstruve0", "x", r_itmodstruve0, D(U(0, 30)), D(LU(1e-300, 700)))

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

section("Airy functions")
for i, f in enumerate(["Ai", "Ai'", "Bi", "Bi'"]):
    add(f"airy[{i}]", "x", lambda x, i=i: _airy(x)[i], D(U(-20, 20)), D(LU(1e-300, 1e6, 0)), note="airy")
for i in range(4):
    dom = D(U(0, 20)) if i < 2 else D(U(-20, 20))
    dom2 = D(LU(1e-300, 1e6)) if i < 2 else D(LU(1e-300, 1e6, 0))
    add(f"airy_scaled[{i}]", "x", lambda x, i=i: _airy(x, True)[i], dom, dom2, note="airy")
for i in range(4):
    add(
        f"airy_integrals[{i}]",
        "x",
        lambda x, i=i: _airy_int(x)[i],
        D(U(0, 30)),
        D(LU(1e-300, 1e3)),
        note="zeros_c" if i == 3 else "",
    )

section("Parabolic cylinder functions")
for f, r in [("pbdv", r_pbdv), ("pbvv", r_pbvv)]:
    large = [D(LU(20, 300, 0), U(-60, 60), n=2000)] if f == "pbdv" else []
    for i in range(2):
        add(f"{f}.{i}", "v x", lambda v, x, r=r, i=i: r(v, x)[i], D(U(-10, 10), U(-10, 10), n=2000), *large, note="zeros_c")
for i in range(2):
    add(f"pbwa.{i}", "a x", lambda a, x, i=i: r_pbwa(a, x)[i], D(U(-5, 5), U(-5, 5), n=2000), note="zeros_c")

section("Zeta functions")
add(
    "zeta",
    "x q",
    r_hurwitz_zeta,
    D(U(1, 20), U(0, 10)),
    D(Offset(1, LU(1e-15, 1e3)), LU(1e-300, 1e300)),
)
add("riemann_zeta", "x", r_riemann_zeta, D(U(-20, 20)), D(LU(1e-300, 1e3, 0)), note="zeros")
add("zetac", "x", r_zetac, D(U(-20, 20)), D(LU(1e-300, 1e3, 0)), note="zeros")

add(
    "gen_harmonic",
    "n a",
    r_gen_harmonic,
    D(Int(1, 1000), U(-5, 5)),
    D(ILU(1, 10**15), Offset(1, LU(1e-2, 20))),
)


def ngh_sampler(r):
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

section("Exponential integrals")
add("expi", "x", r_expi, D(U(-20, 20)), D(LU(1e-300, 700, 0)))
add("expn", "n x", r_expn, D(Int(0, 10), U(0, 50)), D(Int(0, 1000), LU(1e-300, 700)))

add("exp1", "x", r_exp1, D(U(0, 20)), D(LU(1e-300, 750)))
add("scaled_exp1", "x", r_scaled_exp1, D(U(0, 20)), D(LU(1e-300, 1e300)))

section("Sine, cosine, and Fresnel integrals")
add("sici.0", "x", r_si, D(U(-30, 30)), D(LU(1e-300, 1e300, 0)))
add("sici.1", "x", r_ci, D(U(-30, 30)), D(LU(1e-300, 1e300, 0)), note="zeros_c")
add("shichi.0", "x", r_shi, D(U(-30, 30)), D(LU(1e-300, 1e3, 0)))
add("shichi.1", "x", r_chi, D(U(-30, 30)), D(LU(1e-300, 1e3, 0)), note="zeros_c")
add("fresnel.0", "x", r_fresnel_s, D(U(-10, 10)), D(LU(1e-300, 1e300, 0)))
add("fresnel.1", "x", r_fresnel_c, D(U(-10, 10)), D(LU(1e-300, 1e300, 0)))

section("Incomplete gamma and beta functions")
for f, r in [("gammainc", r_gammainc), ("gammaincc", r_gammaincc)]:
    add(
        f,
        "a x",
        r,
        D(U(0, 20), U(0, 50)),
        D(LU(1e-300, 1e5), LU(1e-300, 1e5)),
        A_PM_SQRT_A,
        timeout=60,
    )
add(
    "betainc",
    "a b x",
    r_betainc,
    D(U(0, 10), U(0, 10), U(0, 1)),
    D(LU(1e-5, 1e3), LU(1e-5, 1e3), U(0, 1)),
    AB_MU_PM_SIGMA,
    note="beta_dist",
    timeout=60,
)
for f, r in [("fdtr", r_fdtr), ("fdtrc", r_fdtrc)]:
    add(f, "a b x", r, D(LU(1e-3, 1e6), LU(1e-3, 1e6), LU(1e-30, 1e30)), timeout=60)

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
add(
    "betaincinv",
    "a b y",
    r_betaincinv,
    D(U(0, 10), U(0, 10), U(0, 1), n=2000),
    D(LU(1e-3, 1e3), LU(1e-3, 1e3), LU(1e-300, 1), n=2000),
    timeout=120,
    note="inv",
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
add("elliprc", "x y", r_elliprc, D(LU(1e-3, 1e3), LU(1e-3, 1e3)), D(LU(1e-300, 1e300), LU(1e-300, 1e300)))
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

section("Hypergeometric functions")
add("hyp0f0", "z", r_hyp0f0, D(U(-700, 700)), cpp=False)
add("hyp1f0", "a z", r_hyp1f0, D(U(-10, 10), U(-10, 1)), cpp=False)
add("hyp0f1", "b z", r_hyp0f1, D(LU(1e-3, 1e6, 0), LU(1e-3, 1e6, 0)), note="osc_hyp0f1")
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
    timeout=60,
    note="zeros_c",
)

section("Orthogonal polynomials")


def poly(name, args, ref, dists, int_n=50, cpp=True):
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
poly("eval_sh_jacobi", "n alpha beta x", r_sh_jacobi, [U(-1, 10), U(-1, 10), U(0, 1)], cpp=False)
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

section("Legendre functions")


def nm_sampler(r, hi=50):
    n = r.randint(0, hi)
    return n, r.randint(-n, n)


NM = f"{nobreak('n', IN, Int(0, 50).pretty())}, {nobreak('m', IN, f'{{{MINUS}n, {ELLIPSIS}, n}}')}"
add("legendre_p", "n x", r_legendre, D(Int(0, 100), U(-1, 1)), D(Int(0, 30), LU(1, 1e3, 0)), note="zeros_c")
for f, r in [("assoc_legendre_p", r_assoc_legendre_p), ("assoc_legendre_p_norm", r_assoc_legendre_p_norm)]:
    add(
        f,
        "n m x",
        r,
        D(
            sampler=lambda r: (*nm_sampler(r), r.uniform(-1, 1)),
            label="n in {0, ..., 50}, m in {-n, ..., n}, x in [-1, 1]",
            pretty=f"{NM}, {nobreak('x', IN, U(-1, 1).pretty())}",
        ),
        note="zeros_c",
    )
add(
    "sph_legendre_p",
    "n m theta",
    r_sph_legendre_p,
    D(
        sampler=lambda r: (*nm_sampler(r), r.uniform(0, math.pi)),
        label="n in {0, ..., 50}, m in {-n, ..., n}, theta in [0, pi]",
        pretty=f"{NM}, {nobreak(GREEK['theta'], IN, U(0, math.pi).pretty())}",
    ),
    note="zeros_c",
)
add("pmv", "m v x", r_pmv, D(Int(0, 10), U(-50, 50), U(-1, 1), n=5000), note="zeros_c")

section("Discrete distributions")


def kn_sampler(r, n_hi=1000):
    n = int(10 ** r.uniform(0, math.log10(n_hi)))
    return r.randint(0, n), n


def k_near_np(r):
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
        sampler=lambda r: (*(lambda k, n: (min(k, n - 1), n))(*kn_sampler(r)), r.uniform(0, 1)),
        label="n in {1, ..., 1000} (log), k in {0, ..., n - 1}, y in [0, 1]",
        pretty=f"{nobreak('n', IN, '{1, ' + ELLIPSIS + ', 1000}')} (log), "
        + f"{nobreak('k', IN, '{0, ' + ELLIPSIS + f', n {MINUS} 1}}')}, "
        + nobreak("y", IN, U(0, 1).pretty()),
        n=2000,
    ),
    timeout=120,
    note="inv",
)
for f, r in [("nbdtr", r_nbdtr), ("nbdtrc", r_nbdtrc)]:
    add(
        f,
        "k n p",
        r,
        D(Int(0, 1000), Int(1, 1000), U(0, 1), n=5000),
        D(ILU(1, 2 * 10**9), ILU(1, 2 * 10**9), LU(1e-300, 1), n=2000),
        timeout=60,
    )
add("nbdtri", "k n y", r_nbdtri, D(Int(0, 1000), Int(1, 1000), U(0, 1), n=2000), timeout=120, note="inv")
for f, r in [("pdtr", r_pdtr), ("pdtrc", r_pdtrc)]:
    add(
        f,
        "k m",
        r,
        D(Int(0, 100), U(0, 100), n=5000),
        D(
            sampler=lambda r: (lambda m: (max(0, round(m + r.uniform(-6, 6) * m**0.5)), m))(10 ** r.uniform(1, 9)),
            label="m in [10, 1e9] (log), k in m +- 6 sqrt(m)",
            pretty=f"{nobreak('m', IN, LU(10, 1e9).pretty())} (log), {nobreak('k', IN, f'm {PM} 6{SQRT}m')}",
            n=2000,
        ),
        timeout=60,
    )
add("pdtri", "k y", r_pdtri, D(Int(0, 100), U(0, 1), n=2000), timeout=120, note="inv")

section("Continuous distributions")


def gamma_trans(r):
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
add("gdtria", "p b x", r_gdtria, D(U(0, 1), LU(1e-2, 1e2), LU(1e-2, 1e2), n=2000), timeout=120, note="inv")
add("gdtrib", "a p x", r_gdtrib, D(LU(1e-2, 1e2), U(0, 1), LU(1e-2, 1e2), n=2000), timeout=120, note="inv")
add("gdtrix", "a b p", r_gdtrix, D(LU(1e-2, 1e2), LU(1e-2, 1e2), U(0, 1), n=2000), timeout=120, note="inv")
for f, r in [("chdtr", r_chdtr), ("chdtrc", r_chdtrc)]:
    add(
        f,
        "v x",
        r,
        D(LU(1e-3, 1e3), LU(1e-3, 1e3), n=5000),
        D(
            sampler=lambda r: (lambda v: (v, max(0.0, v + r.uniform(-6, 6) * (2 * v) ** 0.5)))(10 ** r.uniform(1, 8)),
            label="v in [10, 1e8] (log), x in v +- 6 sqrt(2v)",
            pretty=f"{nobreak('v', IN, LU(10, 1e8).pretty())} (log), {nobreak('x', IN, f'v {PM} 6{SQRT}(2v)')}",
            n=2000,
        ),
        timeout=60,
    )
add("chdtri", "v y", r_chdtri, D(LU(1e-3, 1e3), U(0, 1), n=2000), timeout=120, note="inv")
add("fdtri", "a b y", r_fdtri, D(LU(1e-3, 1e3), LU(1e-3, 1e3), U(0, 1), n=2000), timeout=120, note="inv")
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
    cpp=False,
    timeout=120,
    note="inv",
)
add("tukeylambdacdf", "x lambda", r_tukeylambdacdf, D(U(-5, 5), U(-2, 2), n=2000), timeout=60)
add(
    "von_mises_cdf",
    "k x",
    r_von_mises_cdf,
    D(U(0, 10), U(-math.pi, math.pi), n=2000),
    D(LU(1e-3, 1e5), U(-math.pi, math.pi), n=2000),
    timeout=60,
)
for f, r in [("geninvgauss_pdf", r_geninvgauss_pdf), ("geninvgauss_logpdf", r_geninvgauss_logpdf)]:
    add(f, "x p b", r, D(LU(1e-3, 1e3), U(-10, 10), LU(1e-3, 1e3)), note="zeros_c" if "log" in f else "")


def gh_sampler(r):
    a = 10 ** r.uniform(-2, 2)
    return r.uniform(-10, 10), r.uniform(-5, 5), a, a * r.uniform(-0.99, 0.99)


for f, r in [("genhyperbolic_pdf", r_genhyperbolic_pdf), ("genhyperbolic_logpdf", r_genhyperbolic_logpdf)]:
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
add("cdf_cvm_inf", "x", r_cdf_cvm_inf, D(U(0, 2)), D(LU(1e-3, 100)))

section("Kolmogorov-Smirnov distributions")
for f, r in [("smirnov", r_smirnov), ("smirnovc", r_smirnovc), ("smirnovp", r_smirnovp)]:
    add(f, "n x", r, D(ILU(1, 1000), U(0, 1), n=5000), D(ILU(1, 1000), LU(1e-300, 1), n=2000), timeout=60)
for f, r in [("smirnovi", r_smirnovi), ("smirnovci", r_smirnovci)]:
    add(f, "n q" if f == "smirnovi" else "n p", r, D(ILU(1, 1000), U(0, 1), n=1000), timeout=120, note="inv")
for f, r in [("kolmogorov", r_kolmogorov), ("kolmogc", r_kolmogc), ("kolmogp", r_kolmogp)]:
    add(f, "x", r, D(U(0, 3)), D(LU(1e-3, 30)))
add("kolmogi", "q", r_kolmogi, D(U(0, 1), n=2000), D(LU(1e-300, 1), n=2000), note="inv")
add("kolmogci", "p", r_kolmogci, D(U(0, 1), n=2000), D(LU(1e-300, 1), n=2000), note="inv")

section("Other functions")
add("lambertw", "x", r_lambertw, D(U(INV_E, 10)), D(LU(1e-300, 1e300)))
add("spence", "x", r_spence, D(U(0, 2)), D(LU(1e-300, 1e300)))
add("wrightomega", "x", r_wrightomega, D(U(-20, 20)), D(LU(1e-300, 1e300, 0)))
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
add("wright_bessel", "a b x", r_wright_bessel, D(U(0, 5), U(0, 5), U(0, 100), n=2000), timeout=60)
add(
    "log_wright_bessel",
    "a b x",
    r_log_wright_bessel,
    D(U(0, 5), U(0, 5), LU(1e-3, 1e3), n=2000),
    timeout=60,
    note="zeros_c",
)

section("Complex functions")
CZ = D(U(-20, 20), U(-20, 20), n=5000)
for f, r in [("gamma", mp.gamma), ("rgamma", mp.rgamma), ("loggamma", mp.loggamma), ("digamma", mp.digamma)]:
    add(f"{f}#c", "re im", _c(r), CZ, complex=True, note="zeros_c" if f in ("loggamma", "digamma") else "")
for f, r in [("erf", mp.erf), ("erfc", mp.erfc), ("wofz", c_wofz)]:
    add(f"{f}#c", "re im", _c(r), D(U(-10, 10), U(-10, 10), n=5000), complex=True)
for f, r in [("exp1", mp.e1), ("expi", mp.ei), ("lambertw", c_lambertw)]:
    add(f"{f}#c", "re im", _c(r), CZ, complex=True)
add("hyp1f1#c", "a b re im", c_hyp1f1, D(U(-10, 10), U(-10, 10), U(-20, 20), U(-20, 20), n=2000), complex=True, timeout=60)
add("hyp2f1#c", "a b c re im", c_hyp2f1, D(U(-5, 5), U(-5, 5), U(-5, 5), U(-2, 2), U(-2, 2), n=2000), complex=True, timeout=60)
for f, r in [
    ("bessel_j", c_besselj),
    ("bessel_y", c_bessely),
    ("bessel_i", c_besseli),
    ("bessel_k", c_besselk),
    ("hankel_1", c_hankel1),
    ("hankel_2", c_hankel2),
]:
    add(f"{f}#c", "v re im", r, D(U(-20, 20), U(-50, 50), U(-50, 50), n=5000), complex=True, note="osc_c" if f in ("bessel_j", "bessel_y") else "")
for f, r in [("hankel_1e", c_hankel1e), ("hankel_2e", c_hankel2e)]:
    add(f, "v re im", r, D(U(-20, 20), U(-50, 50), U(-50, 50), n=5000), complex=True)
for f, r in [("hankel_1_prime", c_hankel1p), ("hankel_2_prime", c_hankel2p)]:
    add(f, "v re im n", r, D(U(-20, 20), U(-50, 50), U(-50, 50), Int(1, 3), n=5000), complex=True, cpp=False)
for i, f in enumerate(["plus.0", "plus.1", "minus.0", "minus.1"]):
    add(
        f"modified_fresnel_{f}",
        "x",
        lambda x, i=i: c_modified_fresnel(x)[i],
        D(U(0, 10)),
        D(LU(1e-300, 1e4)),
        complex=True,
    )
add(
    "sph_harm_y",
    "n m theta phi",
    c_sph_harm_y,
    D(
        sampler=lambda r: (*nm_sampler(r), r.uniform(0, math.pi), r.uniform(0, 2 * math.pi)),
        label="n in {0, ..., 50}, m in {-n, ..., n}, theta in [0, pi], phi in [0, 2 pi]",
        pretty=f"{NM}, {nobreak(GREEK['theta'], IN, U(0, math.pi).pretty())}, "
        + nobreak(GREEK["phi"], IN, "[0, 2\u03c0]"),
        n=5000,
    ),
    complex=True,
    note="zeros_c",
)

NOTES = {
    "zeros_c": "Ill-conditioned near the zeros.",
    "zeros": "Ill-conditioned near the zeros, where the largest errors are.",
    "osc": f"For x > |v|: A = {SQRT}(J\u1d65(x)\u00b2 + Y\u1d65(x)\u00b2).",
    "neg": (
        f"Ill-conditioned for negative non-integer a with large |a|, and near the poles of "
        f"{GAMMA}(a + b) (or {GAMMA}(a + m)), where the largest errors are."
    ),
    "osc_hyp0f1": (
        f"For z < 0 with x = 2{SQRT}{MINUS}z > |b {MINUS} 1|: A = |{GAMMA}(b)| (x/2)^(1 {MINUS} b) "
        f"{SQRT}(J_(b{MINUS}1)(x)\u00b2 + Y_(b{MINUS}1)(x)\u00b2)."
    ),
    "beta_dist": f"{MU} = a / (a + b), {SIGMA}\u00b2 = {MU}(1 {MINUS} {MU}) / (a + b + 1).",
}

NOTES |= {
    "inv": (
        "Reference: the mpmath root of the forward reference. Ill-conditioned where the forward "
        "function is flat, e.g. in its tails."
    ),
    "osc_sph": (
        f"For x > n + {HALF}: A = {SQRT}(\u03c0/(2x)) {SQRT}(J_(n+{HALF})(x)\u00b2 + "
        f"Y_(n+{HALF})(x)\u00b2)."
    ),
    "osc_c": f"A = {SQRT}(|J\u1d65(z)|\u00b2 + |Y\u1d65(z)|\u00b2).",
    "kelvin": "A = the modulus of the complex function, e.g. |ber x + i bei x| for `ber` and `bei`.",
    "airy": (
        f"For x < 0: A = {SQRT}(Ai(x)\u00b2 + Bi(x)\u00b2), or {SQRT}(Ai\u2032(x)\u00b2 + "
        "Bi\u2032(x)\u00b2) for the derivatives."
    ),
    "binom_dist": f"{SIGMA}\u00b2 = np(1 {MINUS} p).",
    "tvn": "With `epsi` = 10\u207b\u00b9\u2074, its lower bound.",
    "abs": "Absolute error ~10\u207b\u00b9\u2076 at best, so the relative error is large for tiny probabilities.",
    "poly": (
        f"Rows with n \u2208 {{{ELLIPSIS}}} use an integer degree (`i32`, or `u32` for Hermite), "
        "the others an `f64` degree."
    ),
}

NOT_COVERED = """\
## Not covered

- No mpmath reference: Mathieu (`mathieu_*`) and spheroidal wave functions (`prolate_*`,
  `oblate_*`)
- Zeros and coefficients: `*_zeros`, `bernoulli`, `euler`, `mathieu_*_coef`
- All orders or degrees: `*_all`, `riccati_j`, `riccati_y`, `jahnke_emden_lambda`
- Array input: `logsumexp`, `softmax`, `log_softmax`, `cevalpoly`, `poisson_binom_*`, `wilcoxon_*`
- Exact integer functions: `comb`, `comb_rep`, `factorial`, `multifactorial` (and `_checked`),
  `round`
- Other: `cdf_cvm` (an approximation itself), `lanczos_sum_expg_scaled` (internal), the
  `extended_*_error` metrics, `kelvin` (its parts are above), and complex variants not in the
  last table
"""

BY_NAME = {s.name: s for s in SPECS}
