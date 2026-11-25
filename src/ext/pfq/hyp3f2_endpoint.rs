// Positive 3F2 at unity by a telescoping remainder.
//
// For the positive series terms t_n, the residual of F(n) is
// 1 - F(n) + (t_(n+1)/t_n) F(n+1). Multiplying by t_n and summing from N
// gives the difference between the tail and t_N F(N), since the positive excess
// implies t_n F(n) tends to zero. If the residual is uniformly bounded by eta,
// the difference is at most eta times the tail.
use super::{Estimate, PFQ_TOL};
use crate::utils::{next_down_nonneg, next_up_nonneg};

const ORDER: usize = 6;
const N: u32 = 256;
const X: f64 = 1.0 / N as f64;
const U: f64 = f64::EPSILON / 2.0;

fn up(x: f64) -> f64 {
    if x < 0.0 {
        -next_down_nonneg(-x)
    } else {
        next_up_nonneg(x)
    }
}

#[derive(Clone, Copy)]
struct Bounds {
    lo: f64,
    hi: f64,
}

impl Bounds {
    const ZERO: Self = Self::point(0.0);
    const ONE: Self = Self::point(1.0);

    const fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }

    fn outward(lo: f64, hi: f64) -> Self {
        Self {
            lo: -up(-lo),
            hi: up(hi),
        }
    }

    fn is_zero(self) -> bool {
        self.lo == 0.0 && self.hi == 0.0
    }

    fn add(self, other: Self) -> Self {
        if self.is_zero() {
            return other;
        }
        if other.is_zero() {
            return self;
        }
        Self::outward(self.lo + other.lo, self.hi + other.hi)
    }

    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn scale_add(self, vector: Self, factor: f64) -> Self {
        if vector.is_zero() || factor == 0.0 {
            return self;
        }
        let (lo, hi) = if factor > 0.0 {
            (vector.lo, vector.hi)
        } else {
            (vector.hi, vector.lo)
        };
        Self::outward(lo.mul_add(factor, self.lo), hi.mul_add(factor, self.hi))
    }

    fn abs(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }

    fn mid(self) -> f64 {
        0.5 * (self.lo + self.hi)
    }
}

fn multiply(poly: &mut [Bounds], degree: usize, factor: f64) {
    for j in (1..=degree + 1).rev() {
        poly[j] = poly[j].scale_add(poly[j - 1], factor);
    }
}

fn product(factors: impl IntoIterator<Item = f64>) -> [Bounds; ORDER + 4] {
    let mut poly = [Bounds::ZERO; ORDER + 4];
    poly[0] = Bounds::ONE;
    for (degree, factor) in factors.into_iter().enumerate() {
        multiply(&mut poly, degree, factor);
    }
    poly
}

fn correction_polynomial(a: [f64; 3], b: [f64; 2], excess: f64) -> (Bounds, f64) {
    // With x=1/n and term ratio A(x)/B(x), the residual of F(n)=P(x)/x is
    // 1-F(n)+A(x)/B(x)*F(n+1) = R(x)/(x*B(x)*(1+x)^(ORDER-1)).
    let mut den = product(b);
    let upper = product(a);
    let mut residual = [Bounds::ZERO; 2 * ORDER + 3];
    residual[..3].copy_from_slice(&den[..3]);
    let mut p = [0.0; ORDER + 1];
    // At step j, store R(x)/x divided by the remaining factor (1+x)^(ORDER-j).
    // Its divided basis is (den - upper)/x and has constant excess+j.
    for (j, coefficient) in p.iter_mut().enumerate() {
        *coefficient = residual[j].mid() / (excess + j as f64);
        let degree = (j + 2).max(3);
        for k in 1..=degree {
            let basis = den[k].add(upper[k].neg());
            residual[j + k - 1] = residual[j + k - 1].scale_add(basis, -*coefficient);
        }
        if j < ORDER {
            multiply(&mut den, j + 2, 1.0);
            multiply(&mut residual, j + degree - 1, 1.0);
        }
    }
    let poly = p
        .into_iter()
        .rev()
        .fold(Bounds::ZERO, |v, c| Bounds::point(c).scale_add(v, X));
    let eta = residual
        .into_iter()
        .rev()
        .fold(0.0_f64, |v, c| next_up_nonneg(v.mul_add(X, c.abs())));
    (poly, eta)
}

// For positive sums with at most N additions and eleven roundings per term update, invert
// the relative rounding factors. A reported term t_k differs from its exact value by at most
// |t_k|*11k*u/(1-11N*u). The weighted sum accounts for the different errors of the terms.
fn prefix(a: [f64; 3], b: [f64; 2]) -> (f64, f64, f64) {
    let (mut term, mut sum, mut weighted) = (1.0_f64, 0.0, 0.0);
    for k in 0..N {
        let k = f64::from(k);
        sum += term;
        weighted = term.mul_add(k, weighted);
        let num = (a[0] + k) * (a[1] + k) * (a[2] + k);
        let den = (k + 1.0) * (b[0] + k) * (b[1] + k);
        term *= num / den;
    }
    let numerator = next_up_nonneg(weighted.mul_add(11.0 * U, sum * (f64::from(N) * U)));
    let sum_error = next_up_nonneg(numerator / (1.0 - f64::from(12 * N + 3) * U));
    (sum, term, sum_error)
}

#[cold]
pub(super) fn estimate(a: [f64; 3], b: [f64; 2]) -> Option<Estimate> {
    // This range keeps the positive prefix terms normal. The exact ratio is between
    // ((k+.5)/(k+32))^3 and its inverse. Bound these by powers of two on each dyadic
    // block of k to put the exact t_1..t_256 between 2^-613 and 2^613.
    if !a.into_iter().chain(b).all(|x| (0.5..=32.0).contains(&x)) {
        return None;
    }
    let excess = b
        .into_iter()
        .map(Bounds::point)
        .chain(a.into_iter().map(|x| Bounds::point(-x)))
        .fold(Bounds::ZERO, Bounds::add);
    if excess.lo <= 0.0 {
        return None;
    }
    let (poly, eta) = correction_polynomial(a, b, excess.mid());
    if !(poly.lo > 0.0 && eta < 1.0) {
        return None;
    }
    let (sum, term, sum_error) = prefix(a, b);
    let r = 11.0 * f64::from(N) * U;
    let scaled_term = term * f64::from(N);
    let factor_lo = next_down_nonneg(poly.lo * (1.0 - r));
    let factor_hi = next_up_nonneg(poly.hi / (1.0 - r));
    // The analytic remainder bounds the relative error of the tail correction against the
    // true tail. Enclose the true tail by dividing by 1+eta and 1-eta, respectively.
    let tail_lo =
        next_down_nonneg(next_down_nonneg(scaled_term * factor_lo) / next_up_nonneg(1.0 + eta));
    let tail_hi =
        next_up_nonneg(next_up_nonneg(scaled_term * factor_hi) / next_down_nonneg(1.0 - eta));
    let lo = next_down_nonneg(next_down_nonneg((sum - sum_error).max(0.0)) + tail_lo);
    let hi = next_up_nonneg(next_up_nonneg(sum + sum_error) + tail_hi);
    let y = sum + scaled_term * poly.mid();
    if !(lo > 0.0 && hi.is_finite() && lo <= y && y <= hi) {
        return None;
    }
    let err = next_up_nonneg(next_up_nonneg((y - lo).max(hi - y)) / lo);
    (0.0..=PFQ_TOL)
        .contains(&err)
        .then(|| Estimate::new(y, err))
}
