"""SciPy (or NumPy) versions of the functions that xsf lacks, for the second column of the table.

Each is called on the same f64 inputs, so that the column shows what a SciPy user gets. Only
where SciPy has the same function, with the same definition and argument convention; so not:

| xsf-rust                     | SciPy                                        |
| ---------------------------- | -------------------------------------------- |
| `eval_hermite_h(e)`, f64 `n` | integer `n` only                             |
| `eval_sh_jacobi`             | the generalized Gegenbauer polynomial        |
| `mathieu_ce`, `mathieu_se`   | `mathieu_cem` and `mathieu_sem` take degrees |
| `perm`                       | 0 for `k > n`, not the falling factorial     |
"""

import warnings
from collections.abc import Callable, Sequence
from typing import Final

import numpy as np
import scipy
import scipy.special as sp
import scipy.stats as st

type Args = tuple[float, ...]
type Values = Sequence[object] | np.ndarray

# keys of the footnotes in `NOTES`: what SciPy uses
BOOST: Final = "scipy_boost"
CDFLIB: Final = "scipy_cdflib"
XSF: Final = "scipy_xsf"
OWN: Final = "scipy"  # its own Python or Cython code
NUMPY: Final = "numpy"

NOTES: Final = {
    BOOST: "The second column is SciPy, which uses Boost.Math here.",
    CDFLIB: "The second column is SciPy, which uses its cdflib here.",
    XSF: "The second column is SciPy, which wraps xsf here.",
    OWN: "The second column is SciPy.",
    NUMPY: "The second column is NumPy.",
}
VERSIONS: Final = {"scipy": scipy.__version__, "numpy": np.__version__}


def cols(
    f: Callable[..., np.ndarray], ints: Sequence[int] = ()
) -> Callable[[list[Args]], Values]:
    """`f` of the argument columns, with those at `ints` as integers."""

    def g(args: list[Args]) -> Values:
        xs = list(np.array(args).T)
        for i in ints:
            xs[i] = xs[i].astype(np.int64)
        return f(*xs)

    return g


def each(f: Callable[..., object]) -> Callable[[list[Args]], Values]:
    """`f` of each input, for arguments that are arrays or must be scalars."""
    return lambda args: [f(*a) for a in args]


def deriv(f: Callable[..., object]) -> Callable[[list[Args]], Values]:
    """A derivative `f(v, z, n)` of a Bessel function of real `z`, with a scalar `n`."""
    return each(lambda v, x, n: f(v, x, int(n)))


def deriv_c(f: Callable[..., object]) -> Callable[[list[Args]], Values]:
    """The same, of complex `z` = `re` + i `im`."""
    return each(lambda v, re, im, n: f(v, complex(re, im), int(n)))


# by spec name: the SciPy (or NumPy) version, and the key of its footnote
FUNCS: Final[dict[str, tuple[Callable[[list[Args]], Values], str]]] = {
    "ellip_harm": (cols(sp.ellip_harm, ints=(2, 3)), OWN),
    "ellip_harm_2": (cols(sp.ellip_harm_2, ints=(2, 3)), OWN),
    "ellip_normal": (cols(sp.ellip_normal, ints=(2, 3)), OWN),
    "bessel_j_prime": (deriv(sp.jvp), OWN),
    "bessel_y_prime": (deriv(sp.yvp), OWN),
    "bessel_i_prime": (deriv(sp.ivp), OWN),
    "bessel_k_prime": (deriv(sp.kvp), OWN),
    "bessel_j_prime#c": (deriv_c(sp.jvp), OWN),
    "bessel_y_prime#c": (deriv_c(sp.yvp), OWN),
    "bessel_i_prime#c": (deriv_c(sp.ivp), OWN),
    "bessel_k_prime#c": (deriv_c(sp.kvp), OWN),
    "hankel_1_prime": (deriv_c(sp.h1vp), OWN),
    "hankel_2_prime": (deriv_c(sp.h2vp), OWN),
    "bdtrik": (cols(sp.bdtrik), BOOST),
    "bdtrin": (cols(sp.bdtrin), BOOST),
    "btdtria": (cols(sp.btdtria), BOOST),
    "btdtrib": (cols(sp.btdtrib), BOOST),
    "fdtridfd": (cols(sp.fdtridfd), CDFLIB),
    "ncfdtr": (cols(sp.ncfdtr), BOOST),
    "ncfdtri": (cols(sp.ncfdtri), BOOST),
    "ncfdtridfd": (cols(sp.ncfdtridfd), CDFLIB),
    "ncfdtridfn": (cols(sp.ncfdtridfn), CDFLIB),
    "ncfdtrinc": (cols(sp.ncfdtrinc), CDFLIB),
    "nbdtrik": (cols(sp.nbdtrik), BOOST),
    "nbdtrin": (cols(sp.nbdtrin), BOOST),
    "pdtrik": (cols(sp.pdtrik), BOOST),
    "stdtr": (cols(sp.stdtr), BOOST),
    "stdtri": (cols(sp.stdtrit), BOOST),
    "stdtridf": (cols(sp.stdtridf), CDFLIB),
    "nctdtr": (cols(sp.nctdtr), BOOST),
    "nctdtrit": (cols(sp.nctdtrit), BOOST),
    "nctdtridf": (cols(sp.nctdtridf), CDFLIB),
    "nctdtrinc": (cols(sp.nctdtrinc), CDFLIB),
    "chdtriv": (cols(sp.chdtriv), BOOST),
    "chndtr": (cols(sp.chndtr), BOOST),
    "chndtri": (cols(sp.chndtrix), BOOST),
    "chndtridf": (cols(sp.chndtridf), BOOST),
    "chndtrinc": (cols(sp.chndtrinc), BOOST),
    "polygamma": (cols(sp.polygamma), OWN),
    # `d` must be a scalar
    "multigammaln": (each(lambda a, d: sp.multigammaln(a, int(d))), OWN),
    "betaincc": (cols(sp.betaincc), BOOST),
    "betainccinv": (cols(sp.betainccinv), BOOST),
    "diric": (cols(lambda n, x: sp.diric(x, n)), OWN),
    "softplus": (cols(sp.softplus), OWN),
    "powm1": (cols(sp.powm1), BOOST),
    "logaddexp": (cols(np.logaddexp), NUMPY),
    "logaddexp2": (cols(np.logaddexp2), NUMPY),
    "sinc": (cols(np.sinc), NUMPY),
    "logsumexp": (each(lambda *x: sp.logsumexp(x)), OWN),
    "softmax": (each(lambda j, *x: sp.softmax(x)[int(j)]), OWN),
    "log_softmax": (each(lambda j, *x: sp.log_softmax(x)[int(j)]), OWN),
    "poisson_binom_cdf": (each(lambda k, *p: st.poisson_binom.cdf(k, p)), XSF),
}


def evaluate(name: str, args: list[Args]) -> list[complex]:
    """The SciPy (or NumPy) results of the spec `name` (a key of `FUNCS`)."""
    # like the other implementations, without warnings (e.g. of `quad` in `ellip_harm_2`)
    with np.errstate(all="ignore"), warnings.catch_warnings(action="ignore"):
        ys = FUNCS[name][0](args)
    return [
        complex(y) if isinstance(y, complex) else float(y)
        for y in np.asarray(ys).tolist()
    ]
