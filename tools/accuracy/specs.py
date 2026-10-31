"""Functions, domains, and high-precision mpmath references for the accuracy table.

Each `Spec` names an xsf-rust function (the same name is used by the Rust and C++ dump
programs), its argument names, a reference function, and one or more `Domain`s to sample.

A reference function takes `mpf` arguments (the exact values of the f64 inputs) and returns
either the exact function value, or a `(value, scale)` tuple. With a scale, the error is
`min(|y - f| / |f|, |y - f| / scale)`; this is used for oscillatory functions, where the
relative error is meaningless near the zeros (see `NOTES`).

Where mpmath is known to be inaccurate or slow, the references use alternative formulas, e.g.
a direct Euler-Maclaurin sum for the Hurwitz zeta function, quadrature for `expn`, the power
series for `gammainc` with large `a`, and a continued fraction for `betainc`.
"""

import math
from dataclasses import dataclass

import mpmath as mp

DBL_MAX = 1.7976931348623157e308

# ---------------------------------------------------------------- sampling

# Unicode math for the table, escaped to keep the source ASCII
IN, MINUS, PM, ELLIPSIS, SQRT = "\u2208", "\u2212", "\u00b1", "\u2026", "\u221a"
HALF, MU, SIGMA, GAMMA, EPSILON, PSI = "\u00bd", "\u03bc", "\u03c3", "\u0393", "\u03b5", "\u03c8"
SUP = str.maketrans(
    "-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079"
)
GREEK = {"phi": "\u03c6"}


def nobreak(*parts):
    """Join the parts of one argument's domain with no-break spaces, and keep a leading sign
    attached to its interval, so that table cells only wrap between the arguments, or after
    a comma within the interval (the last part)."""
    *head, tail = (p.replace(" ", "\u00a0") for p in parts)
    s = "\u00a0".join([*head, tail.replace(",\u00a0", ", ")])
    return s.replace(f"{PM}[", f"{PM}\u2060[").replace(f"{MINUS}[", f"{MINUS}\u2060[")


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
    return isinstance(d, LU) or isinstance(d, Offset) and _is_log(d.dist)


class Domain:
    """Independent per-argument distributions, or a custom `sampler(rng) -> args`."""

    def __init__(self, *dists, n=10000, sampler=None, label="", pretty=""):
        self.dists, self.n, self.sampler = dists, n, sampler
        self.label, self.pretty_label = label, pretty

    def desc(self, argnames):
        """Plain description; also seeds the sampling."""
        if self.label:
            return self.label
        return ", ".join(f"{a} in {d.desc()}" for a, d in zip(argnames, self.dists))

    def pretty(self, argnames):
        """Unicode description for the table; consecutive equal distributions are merged."""
        if self.pretty_label:
            return self.pretty_label
        groups = []
        for a, d in zip(argnames, self.dists):
            if groups and groups[-1][1] == d:
                groups[-1][0].append(a)
            else:
                groups.append(([a], d))
        return ", ".join(
            nobreak(", ".join(GREEK.get(a, a) for a in names), IN, d.pretty())
            + (" (log)" if _is_log(d) else "")
            for names, d in groups
        )

    def sample(self, rng):
        if self.sampler is not None:
            return tuple(float(a) for a in self.sampler(rng))
        return tuple(d.sample(rng) for d in self.dists)


D = Domain


@dataclass
class Spec:
    name: str
    args: tuple
    ref: object
    domains: list
    section: str
    std: bool = False  # has a Rust std equivalent
    dps: int = 50
    timeout: int = 30
    note: str = ""  # space-separated keys into NOTES


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

section("Error function and normal distribution")
add("erf", "x", r_erf, D(U(-6, 6)), D(LU(1e-300, 10, 0)))
add("erfc", "x", r_erfc, D(U(-6, 27)), D(LU(1e-300, 27, 0)))
add("erfcx", "x", r_erfcx, D(U(-26, 26)), D(LU(1e-300, 1e300)))
add("erfinv", "y", r_erfinv, D(U(-1, 1)), D(Offset(1, LU(1e-16, 1), -1)))
add("erfcinv", "y", r_erfcinv, D(U(0, 2)), D(LU(1e-300, 1)))
add("ndtr", "x", r_ndtr, D(U(-10, 10)), D(LU(1e-300, 38, 0)))
add("ndtri", "p", r_ndtri, D(U(0, 1)), D(LU(1e-300, 0.5)))
add("log_ndtr", "x", r_log_ndtr, D(U(-20, 6)), D(LU(1e-300, 1e10, 0)))

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

section("Exponential integrals")
add("expi", "x", r_expi, D(U(-20, 20)), D(LU(1e-300, 700, 0)))
add("expn", "n x", r_expn, D(Int(0, 10), U(0, 50)), D(Int(0, 1000), LU(1e-300, 700)))

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

section("Elliptic integrals")
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

section("Other functions")
add("hyp0f1", "b z", r_hyp0f1, D(LU(1e-3, 1e6, 0), LU(1e-3, 1e6, 0)), note="osc_hyp0f1")
add("lambertw", "x", r_lambertw, D(U(INV_E, 10)), D(LU(1e-300, 1e300)))
add("spence", "x", r_spence, D(U(0, 2)), D(LU(1e-300, 1e300)))

NOTES = {
    "zeros": (
        "The relative error is ill-conditioned near the zeros of the function, which is where "
        "the largest errors come from."
    ),
    "osc": (
        "In the oscillatory region x > |v|, the error is the absolute error divided by the "
        f"local amplitude {SQRT}(J\u1d65(x)\u00b2 + Y\u1d65(x)\u00b2), if that is smaller than "
        "the relative error."
    ),
    "neg": (
        f"For negative non-integer a, these are ill-conditioned for large |a|, and near the poles "
        f"of {GAMMA}(a + b) (or {GAMMA}(a + m)): a relative change {EPSILON} of a changes "
        f"{GAMMA}(a + b) / {GAMMA}(a) by a factor of about 1 + {EPSILON} a ({PSI}(a + b) "
        f"{MINUS} {PSI}(a)). This is where the largest errors come from."
    ),
    "osc_hyp0f1": (
        f"In the oscillatory region z < 0 with x = 2{SQRT}{MINUS}z > |b {MINUS} 1|, where the "
        f"function is {GAMMA}(b) (x/2)^(1 {MINUS} b) J_(b{MINUS}1)(x), the error is the absolute "
        f"error divided by the local amplitude |{GAMMA}(b)| (x/2)^(1 {MINUS} b) "
        f"{SQRT}(J_(b{MINUS}1)(x)\u00b2 + Y_(b{MINUS}1)(x)\u00b2), if that is smaller than the "
        "relative error."
    ),
    "beta_dist": (
        f"{MU} = a / (a + b) and {SIGMA}\u00b2 = {MU}(1 {MINUS} {MU}) / (a + b + 1) are the mean "
        "and variance of the beta distribution."
    ),
}

BY_NAME = {s.name: s for s in SPECS}
