//! Hahn polynomials (DLMF 18.19-18.26)
//!
//! With `a1 = alpha + 1` and `b1 = beta + 1`, which don't cancel for alpha and beta near -1.

use crate::utils::pow2;

/// The `|q|` above which [`recurrence`] rescales, so that its products stay finite
const BIG: f64 = pow2(512);

/// `q_n` of `q_{k+1} = (u q_k - v q_{k-1}) / w`, where `(u, v, w) = coef(k)`
///
/// The values are rescaled by powers of 2, so that the result only overflows if `q_n` does.
fn recurrence(n: u32, mut q0: f64, mut q1: f64, coef: impl Fn(f64) -> (f64, f64, f64)) -> f64 {
    if n == 0 {
        return q0;
    }
    let mut scale = 0;
    for k in 1..n {
        let (u, v, w) = coef(f64::from(k));
        (q0, q1) = (q1, (u * q1 - v * q0) / w);
        if q1.abs() > BIG {
            (q0, q1, scale) = (q0 / BIG, q1 / BIG, scale + 1);
        }
    }
    (0..scale).fold(q1, |q, _| q * BIG)
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
    recurrence(n, 1.0, q1, |k| {
        let (a, c) = hahn_a_c(k, a1, b1, big_n);
        (a + c - x, c, a)
    })
}

/// `Q_n(x)` for integer `x` by the difference equation in `x` (DLMF 18.22.10-11), from
/// `Q_n(0) = q0`
fn hahn_q_x(n: u32, a1: f64, b1: f64, big_n: f64, x: u32, q0: f64) -> f64 {
    let n = f64::from(n);
    let lambda = n * (n - 1.0 + a1 + b1);
    let q1 = q0 * (1.0 - lambda / (a1 * big_n));
    recurrence(x, q0, q1, |k| {
        let (a, c) = ((k + a1) * (k - big_n), k * (k - b1 - big_n));
        (lambda + a + c, c, a)
    })
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
fn hahn_q_big_n(n: u32, a1: f64, b1: f64) -> f64 {
    let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
    (0..n).fold(sign, |q, j| q * ((b1 + f64::from(j)) / (a1 + f64::from(j))))
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
/// where $n = 0, 1, \ldots, N$, see DLMF 18.20.5 [^DLMF]. For $\alpha, \beta > -1$ these are
/// orthogonal on $x \in \\{0, 1, \ldots, N\\}$, see DLMF Table 18.19.1. The result is NaN if
/// `n > N`, or if `alpha` or `beta` isn't greater than -1.
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/18.20.E5>
///
/// # Algorithm
///
/// - $1 \le x \le N - 1$: the three-term recurrence in $n$ (DLMF 18.22.2), since the sum
///   cancels catastrophically. For integer $x$ and $n > \min(x, N - x)$, where that recurrence
///   can be unstable, the difference equation in $x$ (DLMF 18.22.10) from $Q_n(0) = 1$ or
///   $Q_n(N) = (-1)^n (\beta+1)_n / (\alpha+1)_n$, whichever is on the same side of
///   $A_{n-1} + C_{n-1}$ (DLMF 18.22.3) as $x$.
/// - $x < 1$, and $x > N - 1$ (reflected to $N - x < 1$): the sum, whose terms after the first
///   then all have the same sign.
///
/// # See also
/// - [`eval_jacobi`](crate::eval_jacobi): Evaluate Jacobi polynomials $P_n^{(\alpha, \beta)}$,
///   the limit $Q_n(Nx; \alpha, \beta, N) \to P_n^{(\alpha, \beta)}(1 - 2x) /
///   P_n^{(\alpha, \beta)}(1)$ as $N \to \infty$ (DLMF 18.21.5)
/// - [`legendre_p`](crate::legendre_p): Legendre polynomial of the first kind $P_n$, the limit
///   $Q_n(Nx; 0, 0, N) \to P_n(1 - 2x)$ as $N \to \infty$
#[must_use]
#[allow(non_snake_case)]
pub fn hahn_q(n: u32, N: u32, alpha: f64, beta: f64, x: f64) -> f64 {
    if !(alpha > -1.0 && beta > -1.0 && n <= N) {
        return f64::NAN;
    }
    let (a1, b1, nf, big_nf) = (alpha + 1.0, beta + 1.0, f64::from(n), f64::from(N));
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
        if x <= a + c {
            hahn_q_x(n, a1, b1, big_nf, x as u32, 1.0)
        } else {
            let q_big_n = hahn_q_big_n(n, a1, b1);
            hahn_q_x(n, b1, a1, big_nf, (big_nf - x) as u32, q_big_n)
        }
    } else {
        hahn_q_n(n, a1, b1, big_nf, x)
    }
}
