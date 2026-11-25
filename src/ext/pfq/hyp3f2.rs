//! Real generalized hypergeometric function ₃F₂ and its special cases.

use super::{
    Dd, Estimate, OrBetter, PFQ_DD_MAXITER, PFQ_MODERATE, PFQ_TOL, SeriesReal, first_zero,
    pfq_estimate, pfq_real, series, sort_descending,
};
use crate::ext::{Scaled, hyp2f1_poly_dd, mul_gamma_ratio};
use crate::utils::EPSILON_SQ;
use crate::xsf::cephes::dd_real::{frexp_f64, two_sum};

/// Generalized hypergeometric function $_3F_2$ for real arguments
///
/// This function is not available in SciPy.
///
/// # Definition
///
/// For $\abs z < 1$,
///
/// $$
/// \hyp 3 2 {a_1\enspace a_2\enspace a_3} {b_1\enspace b_2} z
/// = \sum\_{n=0}^\infty
///     {\rpow{a_1}{n} \ \rpow{a_2}{n} \ \rpow{a_3}{n} \over \rpow{b_1}{n} \ \rpow{b_2}{n}}
///     {z^n \over n!}
/// ,
/// $$
///
/// with the rising factorial $\rpow \square n$ ([`pow_rising`](crate::pow_rising)) [^1], and
/// the analytic continuation for $z \le -1$ [^2]. A non-positive integer $a_i = -n$ makes it a
/// polynomial of degree at most $n$, defined for all $z$. An $a_i$ equal to a $b_j$ cancels it,
/// also if it is a non-positive integer.
///
/// It is:
///
/// - $+\infty$ if a $b_j$ is a non-positive integer $-m$ and $z \ne 0$, unless it is a
///   polynomial of degree $\lt m$;
/// - $\pm\infty$ at $z = 1$ for finite parameters if $b_1 + b_2 \le a_1 + a_2 + a_3$, the sign
///   of its limit from below, unless it is a polynomial or that can't be determined (then
///   NaN);
/// - $\pm\infty$ where it detects that it overflows, for $0 < z < 1$;
/// - NaN for $z > 1$, where it is complex, unless it is a polynomial or cancels to
///   $(1 - z)^{-a_i}$ with an integer $a_i$ or to a $_2F_1$ that is a rational function of $z$,
///   and if an argument is NaN;
/// - for an infinite argument, the limit as in [`hyppfq`](crate::hyppfq).
///
/// # Returns NaN
///
/// Where its internal error estimate says that none of its digits are correct:
///
/// - for many $1 - 10^{-5} \lesssim z \le 1$, also where an $a_i$ cancels a $b_j$;
/// - for many $z \lesssim -10^{30}$;
/// - if the terms cancel too much (e.g. for large parameters), overflow, or underflow;
/// - for some infinite arguments, as for [`hyppfq`](crate::hyppfq).
///
/// [^1]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.2>
///
/// [^2]: NIST Digital Library of Mathematical Functions <https://dlmf.nist.gov/16.8.E8>
///
/// # See also
/// - [`hyp2f1`](crate::hyp2f1): Gauss' hypergeometric function, $\hyp 2 1 {a_1\enspace a_2} b z$
/// - [`hyp1f1`](crate::hyp1f1): Kummer's confluent hypergeometric function, $\hyp 1 1 a b z$
#[must_use]
pub fn hyp3f2(a1: f64, a2: f64, a3: f64, b1: f64, b2: f64, z: f64) -> f64 {
    Estimate::or_nan(estimate(a1, a2, a3, b1, b2, z))
}

/// `hyp3f2` and an estimate of its relative error
pub(super) fn estimate(a1: f64, a2: f64, a3: f64, b1: f64, b2: f64, z: f64) -> Option<Estimate> {
    if ![a1, a2, a3, b1, b2, z].into_iter().all(f64::is_finite) {
        return pfq_estimate(&[a1, a2, a3], &[b1, b2], z);
    }
    let (mut a, mut b) = ([a1, a2, a3], [b1, b2]);
    if (0.99..=1.0).contains(&z) && a == [1.0; 3] && b == [2.0; 2] {
        // Li₂(z) / z by the series coefficients. The complement is exact, and the
        // Spence reflection has no large cancellation in this interval.
        let y = crate::xsf::cephes::spence(1.0 - z) / z;
        return Some(Estimate::new(y, 8.0 * f64::EPSILON));
    }
    // a common parameter cancels: 2F1, or 1F0 if that has one as well
    let mut pairs = (0..3).flat_map(|i| (0..2).map(move |j| (i, j)));
    let Some((i, j)) = pairs.find(|&(i, j)| a[i] == b[j]) else {
        return clausen(a, b, z).or_else(|| pfq_real(&mut a, &mut b, z));
    };
    let (mut a, mut b) = ([a[(i + 1) % 3], a[(i + 2) % 3]], [b[1 - j]]);
    let sum = match a.iter().position(|&v| v == b[0]) {
        // the closed form of `pfq_estimate` is also real for z > 1 if `a` is an integer
        Some(i) => return pfq_estimate(&[a[1 - i]], &[], z),
        None if z == 1.0 => gauss_sum(a[0], a[1], b[0]).or_better(|| pfq_real(&mut a, &mut b, z)),
        // `pfq_real` is NaN there, as complex, but this 2F1 can be real
        None if z > 1.0 && first_zero(&a).is_infinite() => return pfq_estimate(&a, &b, z),
        None => pfq_real(&mut a, &mut b, z),
    };
    sum.or_better(|| {
        // Retry an inaccurate finite reduced polynomial with the scaled evaluator
        // on the unit interval.
        if !(0.0..=1.0).contains(&z)
            || !sum?.trusted()?.is_finite()
            || !(1.0..=f64::from(PFQ_DD_MAXITER)).contains(&first_zero(&a))
        {
            return None;
        }
        let y = hyp2f1_poly_dd(a[0], a[1], b[0], z)?;
        if y.m == 0.0 && y.err == 0.0 {
            Some(Estimate::exact(0.0))
        } else {
            Estimate::scaled_certain(y).filter(|v| v.is_accurate())
        }
    })
}

// Clausen's square followed by the Gauss quadratic transformation (DLMF 16.12.2, 15.8.18).
fn clausen(mut a: [f64; 3], b: [f64; 2], z: f64) -> Option<Estimate> {
    if !(0.99..=1.0).contains(&z) || a.iter().any(|&k| k <= 0.0) {
        return None;
    }
    sort_descending(&mut a);
    let [hi, mid, lo] = a;
    let twice = 2.0 * mid;
    let c = if b[0] == twice {
        b[1]
    } else if b[1] == twice {
        b[0]
    } else {
        return None;
    };
    if Dd::new(hi) + lo != twice || Dd::new(mid) + 0.5 != c {
        return None;
    }

    let x = 0.5 * (1.0 - (1.0 - z).sqrt());
    let f = pfq_real(&mut [hi, lo], &mut [c], x)?;
    if !(f.y.is_finite() && f.is_accurate()) {
        return None;
    }
    // Euler's positive integral bounds d log(F)/dx by hi/(1−x).
    // Account for the rounding of x and the square along with the 2F1 error.
    let err = 3.0 * f.err + (4.0 * hi + 2.0) * f64::EPSILON;
    let (m, e) = frexp_f64(f.y);
    Estimate::scaled_certain(Scaled::new(m * m, 2 * e, err)?).filter(|v| v.is_accurate())
}

/// Gauss's sum `2F1(a, b; c; 1)` (DLMF 15.4.20), if it converges and is finite
fn gauss_sum(a: f64, b: f64, c: f64) -> Option<Estimate> {
    let c = Dd::new(c);
    let excess = c - a - b;
    if !(excess.hi > 0.0) {
        return None;
    }
    let y = mul_gamma_ratio([c, excess], [c - a, c - b], || Scaled::new(1.0, 0, 0.0))?;
    // its error bound can't tell if it overflows
    Some(Estimate::scaled(y, None)).filter(|sum| sum.y.is_finite())
}

pub(super) fn unit_upper_at_one(a: &[f64], b: &[f64]) -> Option<Estimate> {
    let (&[1.0, u, v] | &[u, 1.0, v] | &[u, v, 1.0], &[n, w]) = (a, b) else {
        return None;
    };
    for (n, w) in [(n, w), (w, n)] {
        if let Some(value) = shift_at_one(u, v, n, w).or_else(|| shift_at_one(v, u, n, w)) {
            return Some(value);
        }
    }
    thomae_at_one(u.max(v), u.min(v), n.min(w), n.max(w))
}

/// A positive finite sum for `3F2(1,u,v; n,w; 1)` with a small integer shift `d = v - n > 0`.
fn shift_at_one(u: f64, v: f64, n: f64, w: f64) -> Option<Estimate> {
    let (d, low) = two_sum(v, -n);
    // At most nine roundings per update; eight epsilons per term bound their accumulation.
    let err = 8.0 * f64::EPSILON * (d + 1.0);
    let degree =
        (low == 0.0 && d > 0.0 && d.fract() == 0.0 && err <= PFQ_TOL).then_some(d as u32)?;
    let (gap, gap_lo) = two_sum(w, -u);
    // Subtract the integer first, exactly when it is close to the high part (Sterbenz's lemma).
    let gap_minus = |j| (gap - j) + gap_lo;
    (u >= f64::MIN_POSITIVE && n >= f64::MIN_POSITIVE && gap_minus(d + 1.0) > 0.0).then_some(())?;

    // Expand (n+k)_d/(n)_d in binom(k,j) and sum each moment by Gauss (DLMF 15.4.20).
    let mut term = (w - 1.0) / gap_minus(1.0);
    let mut sum = term;
    for j in 0..degree {
        let j = f64::from(j);
        let num = (d - j) * (u + j);
        let den = (n + j) * gap_minus(2.0 + j);
        let ratio = num / den;
        term *= ratio;
        (num.is_normal() && den.is_normal() && ratio.is_normal() && term.is_normal())
            .then_some(())?;
        sum += term;
    }
    sum.mul_add(err, sum)
        .is_finite()
        .then(|| Estimate::new(sum, err))
}

/// Thomae's transformation with a unit upper parameter (DLMF 15.8.1 and 16.5.2).
fn thomae_at_one(u: f64, v: f64, n: f64, w: f64) -> Option<Estimate> {
    const MAX_TERMS: u32 = 512;
    if n < 2.0 || v <= n || w <= u || !PFQ_MODERATE.contains(&w) {
        return None;
    }
    let a = two_sum(n, -u);
    let b = two_sum(n, -v);
    let shift = |(hi, lo), k| (hi + k) + lo;
    // Keep exact integer shifts on the finite-sum path, and bound work before the terms settle.
    let unsuitable =
        |v: (f64, f64)| (v.1 == 0.0 && v.0.fract() == 0.0) || shift(v, f64::from(MAX_TERMS)) < 0.0;
    if unsuitable(a) || unsuitable(b) {
        return None;
    }

    // c = 1 + excess. Accept only when the residual of its three additions fits in two words.
    let (gap, gap_lo) = two_sum(w, -u);
    let (c_hi, c_mid) = two_sum(gap, b.0);
    let (lo, e1) = two_sum(gap_lo, b.1);
    let (lo, e2) = two_sum(c_mid, lo);
    if e1 != 0.0 || e2 != 0.0 {
        return None;
    }
    let c = two_sum(c_hi, lo);
    let excess = shift(c, -1.0);
    if excess <= 0.0 {
        return None;
    }

    let mut term = 1.0;
    let mut sum = term;
    let mut mass = 0.0;
    for k in 0..MAX_TERMS {
        let j = f64::from(k);
        let (x, y) = (shift(a, j), shift(b, j));
        // The nonzero shifted numerators are bounded by MAX_TERMS on a 2^-51 grid.
        // Both denominator factors are in [1, 1e30 + MAX_TERMS], so the ratio is normal.
        let next = term * ((x * y) / ((n + j) * shift(c, j)));
        if !next.is_normal() {
            return None;
        }

        // Once both factors are positive, the ratio is at most k/(k+2).
        // Telescoping (k+1)|t_k| bounds the tail from the first omitted term t_(k+1).
        let tail = (j + 2.0) * next.abs();
        if x > 0.0 && y > 0.0 && tail <= PFQ_TOL * sum.abs() {
            // Bound the inherited ratio errors and signed additions, and round the bound up.
            let rounding = f64::EPSILON * (j + 1.0) * (1.0 + 17.0 * mass);
            let outward = 1.0 + 32.0 * f64::EPSILON * (j + 2.0);
            let err = outward * ((tail + rounding) / sum.abs()) + 8.0 * f64::EPSILON;
            if err <= PFQ_TOL {
                let value = ((w - 1.0) / excess) * sum;
                return (value > 0.0 && value.mul_add(err, value).is_finite())
                    .then(|| Estimate::new(value, err));
            }
        }
        term = next;
        sum += term;
        mass += term.abs();
    }
    None
}

// Clausen's identity (DLMF 16.12.2), the quadratic transform (15.8.18), and the
// Legendre representation (18.5.7) give P_n(sqrt(1-z))^2.
pub(super) fn legendre_square(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    let (&[upper, 0.5, lower], &[1.0, 1.0]) = (a, b) else {
        return None;
    };
    let n = (-lower) as u32;
    if !(1..=PFQ_DD_MAXITER).contains(&n)
        || -lower != f64::from(n)
        || upper != f64::from(n + 1)
        || !(0.0 < z && z <= 1.0)
    {
        return None;
    }
    let odd = n % 2 == 1;
    if z == 1.0 && odd {
        return Some(Estimate::exact(0.0));
    }

    // Unlike in xsf, track even P_k and odd P_k / x in double-double, using exact x^2 = 1-z.
    let x2 = 1.0 - Dd::new(z);
    // The unscaled step contracts the norm with square (e_{k-1} - x e_k)^2 + z e_k^2.
    // Its local error enters that norm with its absolute magnitude. Allow 4 eps^2
    // for each of five DD operations and the reciprocal refinement, at a scale
    // of at most (2 + 1) * 2 while the accumulated error in P_k stays below 1.
    let step = 6.0 * (4.0 * EPSILON_SQ) * (3.0 * 2.0);
    let error = f64::from(n - 1) * step / z.sqrt();
    if error >= 1.0 {
        return None;
    }
    let (_, curr) = (1..n).fold((Dd::ONE, Dd::ONE), |(prev, curr), k| {
        let j = f64::from(k);
        let x_curr = if k % 2 == 0 { curr } else { x2 * curr };
        let d = j + 1.0;
        let inv = 1.0 / d;
        let inv = Dd::new(inv) + (-d).mul_add(inv, 1.0) / d;
        (curr, ((2.0 * j + 1.0) * x_curr - j * prev) * inv)
    });
    let square = curr * curr;
    let y = f64::from(if odd { square * x2 } else { square });
    let err = 3.0 * error / y.sqrt() + f64::EPSILON;
    let y = Estimate::rounded(y, err);
    y.is_accurate().then_some(y)
}

// Use the Pfaff–Saalschütz product for a terminating balanced 3F2 at one
// (DLMF 16.4.3) if the defining series loses accuracy.
pub(super) fn saalschutz(a: &[f64], b: &[f64], n: f64) -> Option<Estimate> {
    let [c, d] = *b else { return None };
    let [_, _, _] = *a else { return None };
    let i = a.iter().position(|&v| v == -n)?;
    let (u, v) = (a[(i + 1) % 3], a[(i + 2) % 3]);
    let exact = |x, y| {
        let (s, e) = two_sum(x, y);
        (e == 0.0 && s.is_finite()).then_some(s)
    };
    let cu = exact(c, -u)?;
    if exact(cu, exact(d, -v)?)? != 1.0 - n {
        return None;
    }
    let cv = exact(c, -v)?;
    // (c-u-v)_n = (-1)^n (d)_n. The unit upper parameter cancels the
    // implicit factorial in the coefficient ratio.
    if first_zero(&[cu, cv]) < n {
        return Some(Estimate::exact(0.0));
    }
    let (m, e) = (0..n as u32).fold((Dd::ONE, 0_i64), |(m, e), j| {
        let (r, re) = series::pfq_frexp_ratio::<Dd>(&[1.0, cu, cv], b, f64::from(j));
        let (m, shift) = (-m * r).frexp();
        (m, e + re + i64::from(shift))
    });
    let err = f64::EPSILON + 16.0 * <Dd as SeriesReal>::EPSILON * n;
    Estimate::scaled_certain(Scaled::new(m.hi, e.try_into().ok()?, err)?)
}
