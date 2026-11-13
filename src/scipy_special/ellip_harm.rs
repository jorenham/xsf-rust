//! Translated from scipy/special/_ellip_harm.pxd

use super::tridiag;
use crate::utils::pow2;
use crate::xsf::cephes::dd_real::{frexp_f64, ldexp_f64};

/// `x^2 - c`, without the rounding error of `x^2`, which matters if `x^2 ~ c`
fn sq_sub(x: f64, c: f64) -> f64 {
    // Veltkamp's split `x = hi + lo` (as in Dekker's `two_prod`), where `hi` and `lo` have at
    // most 26 bits, so that `hi^2`, `hi lo`, and `lo^2` are exact
    const SPLITTER: f64 = pow2(27) + 1.0;
    let x2 = x * x;
    if x2.is_infinite() {
        // and `SPLITTER x` might overflow
        return x2;
    }
    let t = SPLITTER * x;
    let hi = t - (t - x);
    let lo = x - hi;
    ((hi * hi - c) + 2.0 * hi * lo) + lo * lo
}

/// The type of a Lamé function, i.e. its factor `psi` in front of the polynomial (DLMF 29.12.1)
#[derive(Clone, Copy)]
enum Kind {
    K,
    L,
    M,
    N,
}

/// Lamé function of the first kind `E_n^p`, as `psi(s)` times a polynomial in `s^2 - h2`
pub(crate) struct Lame {
    h2: f64,
    k2: f64,
    odd: bool,
    kind: Kind,
    /// of ascending powers of `s^2 - h2`
    coefs: Vec<f64>,
}

impl Lame {
    /// `lame_coefficients`, or `None` for an invalid `p`
    ///
    /// Unlike in SciPy, the eigenpair of the symmetric tridiagonal matrix is computed without
    /// LAPACK's `dstevr`, see [`tridiag`].
    pub(crate) fn new(h2: f64, k2: f64, n: u32, p: u32) -> Option<Self> {
        let (n, p) = (n as usize, p as usize);
        let odd = n % 2 == 1;

        let r = n / 2;
        let alpha = h2;
        let beta = k2 - h2;
        let gamma = alpha - beta;

        // the type, the order `tp` within that type, and the number of coefficients
        let (kind, tp, size) = if p == 0 {
            return None;
        } else if p - 1 < r + 1 {
            (Kind::K, p, r + 1)
        } else if p - 1 < (n - r) + (r + 1) {
            (Kind::L, p - (r + 1), n - r)
        } else if p - 1 < (n - r) + (n - r) + (r + 1) {
            (Kind::M, p - (n - r) - (r + 1), n - r)
        } else if p - 1 < 2 * n + 1 {
            (Kind::N, p - (n - r) - (n - r) - (r + 1), r)
        } else {
            return None;
        };

        // the tridiagonal matrix with diagonal `d`, and `g` above and `f` below it
        let mut g = vec![0.0; size];
        let mut f = vec![0.0; size];
        let mut d = vec![0.0; size];
        let r = r as f64;
        for i in 0..size {
            let j = i as f64;
            (g[i], f[i], d[i]) = match (kind, odd) {
                (Kind::K, true) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 1.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0)) + 2.0) * (2.0 * ((j + 1.0) + r) + 1.0),
                    ((2.0 * r + 1.0) * (2.0 * r + 2.0) - 4.0 * j * j) * alpha
                        + (2.0 * j + 1.0) * (2.0 * j + 1.0) * beta,
                ),
                (Kind::K, false) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 1.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0)) + 2.0) * (2.0 * (r + (j + 1.0)) - 1.0),
                    2.0 * r * (2.0 * r + 1.0) * alpha - 4.0 * j * j * gamma,
                ),
                (Kind::L, true) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 3.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0)) + 2.0) * (2.0 * ((j + 1.0) + r) + 1.0),
                    (2.0 * r + 1.0) * (2.0 * r + 2.0) * alpha
                        - (2.0 * j + 1.0) * (2.0 * j + 1.0) * gamma,
                ),
                (Kind::L, false) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 3.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0))) * (2.0 * (r + (j + 1.0)) + 1.0),
                    (2.0 * r * (2.0 * r + 1.0) - (2.0 * j + 1.0) * (2.0 * j + 1.0)) * alpha
                        + (2.0 * j + 2.0) * (2.0 * j + 2.0) * beta,
                ),
                (Kind::M, true) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 1.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0)) + 2.0) * (2.0 * ((j + 1.0) + r) + 1.0),
                    ((2.0 * r + 1.0) * (2.0 * r + 2.0) - (2.0 * j + 1.0) * (2.0 * j + 1.0)) * alpha
                        + 4.0 * j * j * beta,
                ),
                (Kind::M, false) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 1.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0))) * (2.0 * (r + (j + 1.0)) + 1.0),
                    2.0 * r * (2.0 * r + 1.0) * alpha - (2.0 * j + 1.0) * (2.0 * j + 1.0) * gamma,
                ),
                (Kind::N, true) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 3.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0))) * (2.0 * ((j + 1.0) + r) + 3.0),
                    (2.0 * r + 1.0) * (2.0 * r + 2.0) * alpha
                        - (2.0 * j + 2.0) * (2.0 * j + 2.0) * gamma,
                ),
                (Kind::N, false) => (
                    -(2.0 * j + 2.0) * (2.0 * j + 3.0) * beta,
                    -alpha * (2.0 * (r - (j + 1.0))) * (2.0 * (r + (j + 1.0)) + 1.0),
                    2.0 * r * (2.0 * r + 1.0) * alpha - (2.0 * j + 2.0) * (2.0 * j + 2.0) * alpha
                        + (2.0 * j + 1.0) * (2.0 * j + 1.0) * beta,
                ),
            };
        }

        // the diagonal similarity transform `ss` that makes it symmetric
        let mut ss = vec![1.0; size];
        for i in 1..size {
            ss[i] = (g[i - 1] / f[i - 1]).sqrt() * ss[i - 1];
        }
        let dd: Vec<f64> = (0..size - 1).map(|i| g[i] * ss[i] / ss[i + 1]).collect();

        // the eigenvector of the `tp`-th smallest eigenvalue
        let eigv = tridiag::eigvec(&d, &dd, tridiag::eigvals(&d, &dd)[tp - 1]);
        let eigv: Vec<f64> = eigv.iter().zip(&ss).map(|(v, s)| v / s).collect();

        // the coefficients of the polynomial in `s^2 - h2`, monic in `s^2`, so
        // `v_j (-h2)^(size - 1 - j) / v_last`, with the exponents kept apart
        //
        // Unlike in SciPy, not of the one in `1 - s^2 / h2`, whose normalization by
        // `(-h2)^(size - 1)` underflows for small `h2` and large `n`.
        let (m_h, e_h) = frexp_f64(-h2);
        let (m_last, e_last) = frexp_f64(eigv[size - 1]);
        // `(-h2)^(size - 1 - j) / v_last = m 2^e`
        let (mut m, mut e) = (1.0 / m_last, -e_last);
        let mut coefs = vec![0.0; size];
        for (c, &v) in coefs.iter_mut().zip(&eigv).rev() {
            let (m_v, e_v) = frexp_f64(v);
            *c = ldexp_f64(m_v * m, e_v + e);
            let (m_next, e_next) = frexp_f64(m * m_h);
            (m, e) = (m_next, e + e_next + e_h);
        }
        Some(Self {
            h2,
            k2,
            odd,
            kind,
            coefs,
        })
    }

    /// `ellip_harm_eval`
    ///
    /// Unlike in SciPy, `s^2 - h2` and `s^2 - k2` are computed with [`sq_sub`], and their square
    /// roots are taken separately, and without `s^2` if that overflows, so that they can't
    /// overflow needlessly.
    pub(crate) fn eval(&self, s: f64, signm: f64, signn: f64) -> f64 {
        let (h2, k2) = (self.h2, self.k2);
        let (s2_h2, s2_k2) = (sq_sub(s, h2), sq_sub(s, k2));
        // `sqrt(|s^2 - c|)`
        let root = |s2_c: f64, c: f64| {
            if s2_c.is_finite() {
                s2_c.abs().sqrt()
            } else {
                s.abs() * (1.0 - c / s / s).abs().sqrt()
            }
        };
        // `pow(s, n - 2 r)` and `pow(s, 1 - n + 2 r)`
        let (s_odd, s_even) = if self.odd { (s, 1.0) } else { (1.0, s) };
        let psi = match self.kind {
            Kind::K => s_odd,
            Kind::L => s_even * signm * root(s2_h2, h2),
            Kind::M => s_even * signn * root(s2_k2, k2),
            Kind::N => s_odd * signm * signn * root(s2_h2, h2) * root(s2_k2, k2),
        };

        let (&last, rest) = self.coefs.split_last().unwrap();
        let pp = rest.iter().rev().fold(last, |pp, &c| pp * s2_h2 + c);
        pp * psi
    }
}

/// Ellipsoidal harmonic function $E^p_n(s)$
///
/// Also known as the Lamé function of the first kind, a solution of the Lamé equation
///
/// $$
/// (s^2 - h^2)(s^2 - k^2)E^{\prime\prime}(s) + s(2s^2 - h^2 - k^2)E^\prime(s) +
/// (a - n(n+1) s^2)E(s) = 0
/// $$
///
/// for an eigenvalue $a$.
///
/// Corresponds to [`scipy.special.ellip_harm`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellip_harm.html
///
/// # Arguments
/// - `h2`: $h^2$
/// - `k2`: $k^2$, should be larger than `h2`
/// - `n`: degree
/// - `p`: order, `1 <= p <= 2n + 1`
/// - `s`: coordinate
/// - `signm`, `signn`: the sign, `+1.0` or `-1.0`, of the prefactor of the functions of type L
///   (`signm`), M (`signn`), and N (`signm * signn`), see DLMF 29.12.1
///
/// Returns NaN if `p` or a sign is invalid.
///
/// # See also
/// - [DLMF 29.12](https://dlmf.nist.gov/29.12)
#[must_use]
pub fn ellip_harm(h2: f64, k2: f64, n: u32, p: u32, s: f64, signm: f64, signn: f64) -> f64 {
    if signm.abs() != 1.0 || signn.abs() != 1.0 {
        return f64::NAN;
    }
    Lame::new(h2, k2, n, p).map_or(f64::NAN, |e| e.eval(s, signm, signn))
}
