//! Hahn polynomials (DLMF 18.19-18.26)
//!
//! With `a1 = alpha + 1` and `b1 = beta + 1`, which don't cancel for alpha and beta near -1.

use crate::utils::{neg_one_pow, pow2};
use crate::xsf::Real;
use crate::xsf::cephes::dd_real::DoubleDouble as Dd;

/// The `|q|` above which [`recurrence`] rescales, so that its products stay finite
const BIG: f64 = pow2(512);

/// `q_n` of `q_{k+1} = (u q_k - v q_{k-1}) / w`, where `(u, v, w) = coef(k)`
///
/// The values are rescaled by powers of 2, so that the result only overflows if `q_n` does.
fn recurrence<T: Real>(
    n: u32,
    mut q0: T,
    mut q1: T,
    coef: impl Fn(f64) -> (T, T, T),
    cancellation_tol: f64,
) -> (T, bool) {
    if n == 0 {
        return (q0, false);
    }
    let big = T::from_f64(BIG);
    let mut scale = 0;
    let (mut p, mut q, mut difference) = (q0, q0, q0);
    for k in 1..n {
        let (u, v, w) = coef(f64::from(k));
        (p, q) = (u * q1, v * q0);
        difference = p - q;
        (q0, q1) = (q1, difference / w);
        if q1.abs() > big {
            (q0, q1, scale) = (q0 / big, q1 / big, scale + 1);
        }
    }
    let cancels = cancellation_tol > 0.0 && {
        let mass = p.abs().to_f64() + q.abs().to_f64();
        mass.is_normal() && difference.abs().to_f64() <= cancellation_tol * mass
    };
    ((0..scale).fold(q1, |q, _| q * big), cancels)
}

/// `(A_k, C_k)` of the recurrence in `n` (DLMF 18.22.3), for `k >= 1`
fn hahn_a_c(k: f64, a1: f64, b1: f64, big_n: f64) -> (f64, f64) {
    let (ab2, j) = (a1 + b1, k - 1.0);
    let t = 2.0 * j + ab2; // 2k + alpha + beta
    let a = (j + ab2) * (k + a1) * (big_n - k) / ((t + 1.0) * (t + 2.0));
    let c = k * (j + ab2 + big_n) * (j + b1) / (t * (t + 1.0));
    (a, c)
}

/// `Q_n(x)` by the recurrence in `n` (DLMF 18.22.2-3)
fn hahn_q_n(n: u32, a1: f64, b1: f64, big_n: f64, x: f64) -> f64 {
    let q1 = 1.0 - (a1 + b1) * x / (a1 * big_n);
    recurrence(
        n,
        1.0,
        q1,
        |k| {
            let (a, c) = hahn_a_c(k, a1, b1, big_n);
            (a + c - x, c, a)
        },
        0.0,
    )
    .0
}

/// `Q_n(x)` for integer `x` by the difference equation in `x` (DLMF 18.22.10-11), from
/// `Q_n(0) = q0`
fn hahn_q_x<T: Real>(
    n: u32,
    a1: T,
    b1: T,
    big_n: f64,
    x: u32,
    q0: T,
    cancellation_tol: f64,
) -> (T, bool) {
    let n = T::from_f64(f64::from(n));
    let big_n = T::from_f64(big_n);
    let one = T::from_f64(1.0);
    let lambda = n * (n - one + a1 + b1);
    let q1 = q0 * (one - lambda / (a1 * big_n));
    recurrence(
        x,
        q0,
        q1,
        |k| {
            let k = T::from_f64(k);
            let (a, c) = ((k + a1) * (k - big_n), k * (k - b1 - big_n));
            (lambda + a + c, c, a)
        },
        cancellation_tol,
    )
}

/// Repeat the same difference equation with the original shape parameters.
#[cold]
fn hahn_q_x_retry(
    n: u32,
    big_n: f64,
    x: u32,
    alpha: f64,
    beta: f64,
    reflect: bool,
    original: f64,
) -> f64 {
    let (a1, b1) = (Dd::new(alpha) + Dd::ONE, Dd::new(beta) + Dd::ONE);
    let (a1, b1, q0) = if reflect {
        (b1, a1, hahn_q_big_n(n, a1, b1))
    } else {
        (a1, b1, Dd::ONE)
    };
    if !q0.hi.is_normal() {
        return original;
    }
    // Keep the original zero when the replay still cancels near double-double precision.
    let tol = if original == 0.0 { pow2(-96) } else { 0.0 };
    let (q, cancels) = hahn_q_x(n, a1, b1, big_n, x, q0, tol);
    if q.is_finite() && !cancels && (q.hi != 0.0 || original != 0.0) {
        q.hi
    } else {
        original
    }
}

/// `t0 Q_n(x)` for `x < 1` by its 3F2 sum (DLMF 18.20.5), whose terms after the first all have the
/// sign of `-t0 x`
fn hahn_q_sum(n: u32, a1: f64, b1: f64, big_n: f64, x: f64, mut t: f64) -> f64 {
    // For subnormal x, the terms after the first (which are proportional to x, since `k - x = k`
    // for `k >= 1`) are scaled by 2^52, so that they stay normal.
    if x != 0.0 && x.abs() < f64::MIN_POSITIVE {
        let scale = pow2(52);
        let s = hahn_q_sum(n, a1, b1, big_n, x * scale, t);
        return t + (s - t) / scale;
    }
    let n = f64::from(n);
    let mut s = t;
    let mut k = 0.0;
    while k < n {
        // `k - x` last, so that a tiny x doesn't underflow
        t *= (n - k) / (k + 1.0) * ((n + k - 1.0 + a1 + b1) / (k + a1)) / (big_n - k) * (k - x);
        s += t;
        k += 1.0;
    }
    s
}

/// `Q_n(N) = (-1)^n (beta + 1)_n / (alpha + 1)_n`, the factor of the reflection
/// `Q_n(x; alpha, beta, N) = Q_n(N) Q_n(N - x; beta, alpha, N)` (by the symmetry of the weight in
/// DLMF Table 18.19.1)
fn hahn_q_big_n<T: Real>(n: u32, a1: T, b1: T) -> T {
    let sign = T::from_f64(neg_one_pow(n));
    (0..n).fold(sign, |q, j| {
        let j = T::from_f64(f64::from(j));
        q * ((b1 + j) / (a1 + j))
    })
}

/// Evaluate the Hahn polynomial $Q_n(x; \alpha, \beta, N)$ at a point.
///
/// Not in SciPy.
///
/// # Definition
///
/// $$
/// Q_n(x; \alpha, \beta, N) =
///   \hyp{3}{2}{-n,\enspace n+\alpha+\beta+1,\enspace -x}{\alpha+1,\enspace -N}{1}
/// $$
///
/// where $n = 0, 1, \ldots, N$, see [DLMF 18.20.5]. For $\alpha, \beta > -1$ these are
/// orthogonal on $x \in \\{0, 1, \ldots, N\\}$, see [DLMF Table 18.19.1]. The result is NaN if
/// `n > N`, or if `alpha` or `beta` isn't greater than -1.
///
/// As $N \to \infty$, they tend to the Jacobi polynomials ([DLMF 18.21.5]):
///
/// $$
/// \lim_{N \to \infty} Q_n(Nx; \alpha, \beta, N) =
///   \frac{P_n^{(\alpha, \beta)}(1 - 2x)}{P_n^{(\alpha, \beta)}(1)}
/// $$
///
/// so for $\alpha = \beta = 0$ to the Legendre polynomials $P_n(1 - 2x)$.
///
/// [DLMF 18.20.5]: https://dlmf.nist.gov/18.20.E5
/// [DLMF Table 18.19.1]: https://dlmf.nist.gov/18.19.T1
/// [DLMF 18.21.5]: https://dlmf.nist.gov/18.21.E5
///
/// # See also
/// - [`eval_jacobi`](crate::eval_jacobi): Jacobi polynomials $P_n^{(\alpha, \beta)}$
/// - [`legendre_p`](crate::legendre_p): Legendre polynomials of the first kind $P_n$
/// - [`hyp3f2`](crate::hyp3f2): Generalized hypergeometric function $_3F_2$
#[must_use]
#[allow(non_snake_case)]
pub fn eval_hahn(n: u32, N: u32, alpha: f64, beta: f64, x: f64) -> f64 {
    if !(alpha > -1.0 && beta > -1.0 && n <= N) {
        return f64::NAN;
    }
    let (a1, b1, nf, big_nf) = (alpha + 1.0, beta + 1.0, f64::from(n), f64::from(N));
    // x < 1 and x > N - 1 (reflected): the sum, whose terms after the first have one sign; in
    // between, it can cancel catastrophically
    if x < 1.0 {
        hahn_q_sum(n, a1, b1, big_nf, x, 1.0)
    } else if x > big_nf - 1.0 {
        hahn_q_sum(n, b1, a1, big_nf, big_nf - x, hahn_q_big_n(n, a1, b1))
    } else if x.fract() == 0.0 && nf > x.min(big_nf - x) {
        // At integer x, Q_n(x) is a minimal solution of the recurrence in n for large n and x near
        // 0 or N. By the duality DLMF 18.21.1, the difference equation in x is stable there, if
        // it starts on the same side of the center A + C of the oscillatory region in x (where
        // the characteristic roots of the recurrence in n at n - 1 are complex).
        let (a, c) = hahn_a_c(nf - 1.0, a1, b1, big_nf);
        let (x, a, b, q0, reflect) = if x <= a + c {
            (x as u32, a1, b1, 1.0, false)
        } else {
            ((big_nf - x) as u32, b1, a1, hahn_q_big_n(n, a1, b1), true)
        };
        // Bound the work of the higher-precision endpoint product and difference equation.
        let work = x.saturating_add(if reflect { n } else { 0 });
        // A relative error of order 1e-13 can arise from rounding the final products alone.
        let tol = if x > 1 && work <= 128 { 1e-3 } else { 0.0 };
        let (q, retry) = hahn_q_x(n, a, b, big_nf, x, q0, tol);
        if retry {
            hahn_q_x_retry(n, big_nf, x, alpha, beta, reflect, q)
        } else {
            q
        }
    } else {
        hahn_q_n(n, a1, b1, big_nf, x)
    }
}
