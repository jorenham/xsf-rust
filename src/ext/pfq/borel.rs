//! Divergent `p > q + 1` Borel summation and `2F0` evaluation via `hypu`.

use core::cell::OnceCell;
use core::iter::{once, successors};

use super::asymptotic::{Alg, asymptotic_sum, pfq_alg, pfq_alg_dd};
use super::series::{pfq_terms, pfq_weniger_dd};
use super::{
    Dd, Estimate, OrBetter, PFQ_DD_MAXITER, PFQ_INACCURATE, PFQ_TOL, SeriesReal, add_exp,
    ldexp_i64, pfq_settled,
};
use crate::utils::{SQRT_EPSILON, is_nonpos_int, pow2};
use crate::xsf::cephes::dd_real::{frexp_f64, two_prod, two_sum};
use crate::xsf::xsf_hypu_digits;

/// Digits of U beyond those of `PFQ_TOL` that `chgu` must estimate for the Borel sum of 2F0 to be
/// taken before the other methods, as these are often more accurate (measured)
const CHGU_DIGITS_MARGIN: f64 = 1.0;
/// pFq for `p > q + 1` and `z < 0` by `pfq_alg_dd`, or `None` if it may or may not overflow
fn pfq_borel_dd(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    Estimate::scaled_certain(pfq_alg_dd(a, b, z)?)
}

/// pFq for `p > q + 1` and `z < 0`, where the series diverges: its Borel sum
///
/// That is the Meijer G function of DLMF 16.5.1, whose residues at the poles of the
/// `Gamma(a_i + s)` sum to `pfq_alg`. For small `|z|`, where those cancel, it is the series up
/// to its smallest term, as for 2F0, or else its Weniger transformation in double-double. The
/// residues come first, as the steps of Weniger can stall where some terms nearly vanish, e.g.
/// for an `a_i` near a non-positive integer. Where none of these is accurate, the residues are
/// in double-double, as for 2F0.
pub(super) fn pfq_borel(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    // the term ratios are about n^d |z|, so that the terms are smallest near n = |z|^(-1/d)
    let d = (a.len() - b.len() - 1) as f64;
    // (capped as for 2F0)
    let scale = z.abs().powf(-1.0 / d).min(f64::EPSILON.recip());
    let terms = once(1.0).chain(pfq_terms(a, b, z));
    let settled = pfq_settled(a, b).ceil();
    let series = || Some(Estimate::absolute(asymptotic_sum(terms, settled, scale, |_| 0.0)));
    let residues = || Some(Estimate::scaled(add_exp(pfq_alg(a, b, z, Alg::Value)?)?, None));
    series()
        .or_better(residues)
        .or_better(|| pfq_weniger_dd(a, b, z, settled))
        .or_better(|| pfq_borel_dd(a, b, z))
}

// The Laplace kernel (1 + |z| X Y)^(-a) lies in [0,1] and is symmetric in X,Y.
#[cold]
pub(super) fn hyp3f0_near_coincident(a: &[f64], b: &[f64], z: f64) -> Option<Estimate> {
    let mut upper: [f64; 3] = a.try_into().ok()?;
    upper.sort_unstable_by(f64::total_cmp);
    let [lo, mid, hi] = upper;
    let delta = mid - lo;
    // The symmetrized Gamma density ratio is cosh(delta*ln(X/Y))/M,
    // with M = Gamma(lo)*Gamma(hi)/Gamma(mid)^2. The absolute error is at most
    // 1 - 1/M <= ln(M) <= delta² (1/lo² + 1/lo), using the trigamma series (DLMF 5.15.1).
    let ratio = delta / lo;
    let bound = 2.0 * ratio * (ratio + delta); // allow for rounding in the bound
    // The bound forces mid - lo < lo. Equal spacings then give exact subtractions by Sterbenz.
    if b.is_empty() && lo > 0.0 && delta > 0.0 && delta == hi - mid && bound <= PFQ_TOL {
        let mut centered = pfq_borel(&[mid; 3], &[], z).filter(|e| e.y.is_normal() && e.y > 0.0)?;
        centered.err += bound / centered.y;
        centered.is_accurate().then_some(centered)
    } else {
        None
    }
}

// The positive-shape 3F0 at -t is E[(1+tXY)^(-a)] for independent
// unit-scale Gamma variables X and Y with shapes b and c, respectively.
#[cold]
pub(super) fn hyp3f0_gamma_product(upper: &[f64], lower: &[f64], z: f64) -> Option<Estimate> {
    let [c, b, a] = super::descending(<[f64; 3]>::try_from(upper).ok()?);
    let s = b + c;
    // The density of 2 sqrt(XY) vanishes at zero for b >= 1. Keep its power and
    // gamma factors well within binary64 range on the finite integration interval.
    if !(lower.is_empty() && a > 0.0 && b >= 1.0 && s <= 64.0) {
        return None;
    }
    let root = (-z).sqrt() / 2.0;
    // U <= X+Y ~ Gamma(s), and E[exp((X+Y)/2)] = 2^s. At this endpoint
    // Markov's inequality bounds the omitted probability by less than epsilon^2.
    let end = 2.0 * core::f64::consts::LN_2 * (s + 106.0);
    // U = 2 sqrt(XY) has density 2 (u/2)^(b+c-1) K_(b-c)(u)/(Gamma(b) Gamma(c)).
    // root*u = sinh(w) resolves the scale of (1+(root*u)^2)^(-a).
    let hi = (root * end).asinh();
    // Allow for the rounded abscissa in the powers, the Bessel value and normalization.
    let rounding = f64::EPSILON * (64.0 + 3.0 * s + (0.5 * s + a) * hi);
    if rounding >= PFQ_TOL {
        return None;
    }
    let f = |w: f64| {
        let u = w.sinh() / root;
        let k = if b == c {
            crate::bessel_k0(u)
        } else {
            crate::bessel_k(c - b, u)
        };
        (u / 2.0).powf(s - 1.0) * k * w.cosh().powf(1.0 - 2.0 * a)
    };
    let (sum, error) = crate::ext::quadpack::quad(f, 0.0, hi, 0.0, PFQ_TOL / 4.0);
    let norm = 2.0 / (root * crate::gamma(b) * crate::gamma(c));
    let value = norm * sum;
    let tail = f64::EPSILON.powi(2) * hi.cosh().powf(-2.0 * a);
    let result = Estimate::absolute((value, norm * error + tail + rounding * value));
    ((f64::MIN_POSITIVE..=1.0).contains(&value) && result.is_accurate()).then_some(result)
}

/// Cached evaluation of the Borel sum `x^a1 U(a1, a1 - a2 + 1, x)` for `2F0(a1, a2;; -1/x)`,
/// with its finite-difference rounding sensitivity and pole proximity flag
struct RawBorel {
    p: f64,
    u: f64,
    digits: f64,
    rounding: f64,
    near_pole: bool,
}

impl RawBorel {
    fn eval(a1: f64, a2: f64, x: f64, dx: f64) -> Option<Self> {
        const T: f64 = pow2(20);
        const POLE_SPACINGS: f64 = 8.0;
        let (d, e1) = two_sum(a1, -a2);
        let (b, e2) = two_sum(d, 1.0);
        // chgu takes U as a polynomial where its `a1 - b + 1` (that should be a2) rounds to
        // a non-positive integer
        let (s, f1) = two_sum(a1, -b);
        let (a2_rounded, f2) = two_sum(s, 1.0);
        if is_nonpos_int(a2_rounded) && a2_rounded != a2 {
            return None;
        }
        let (p, (u, _, digits)) = (x.powf(a1), xsf_hypu_digits(a1, b, x));
        if !(p.is_normal() && u.is_normal()) {
            return None;
        }
        // The rounding errors of b and x change y by about (y_t - y) / T, with b and x moved
        // by T times them (far enough to stand out from the error of U, and near enough to be
        // linear). That of chgu's `a1 - b + 1` adds to that of b, as chgu then returns about U
        // at `b = a1 - a2_rounded + 1`. A move that would round away is a spacing instead.
        let moved = |v: f64, dv: f64| {
            let v_t = v + T * dv;
            let spacing = (f64::EPSILON * v.abs()).max(f64::MIN_POSITIVE);
            if v_t == v && dv != 0.0 {
                v + spacing.copysign(dv)
            } else {
                v_t
            }
        };
        let db = ((e1 + e2).abs() + (f1 + f2).abs()).copysign(e1 + e2);
        let (b_t, x_t) = (moved(b, db), moved(x, dx));
        let (digits, rounding) = if b_t == b && x_t == x {
            (f64::from(digits), 0.0)
        } else {
            let (p_t, (u_t, _, digits_t)) = (x_t.powf(a1), xsf_hypu_digits(a1, b_t, x_t));
            if !(p_t.is_normal() && u_t.is_normal()) {
                return None;
            }
            // relative to y, also where it overflows
            let rounding = ((p_t / p) * u_t - u).abs() / u.abs() / T;
            (f64::from(digits.min(digits_t)), rounding)
        };
        // Within a few spacings of a non-positive integer n (also 0), the rounding of chgu's
        // a2 can move a cancellation in U, so only an accurate y is kept there (unless U is
        // the polynomial of a2 = n).
        let n = a2_rounded.round().min(0.0);
        let spacing = f64::EPSILON * a1.abs().max(b.abs()).max(1.0);
        let near = (a2 - n).abs().min((a2_rounded - n).abs()) <= POLE_SPACINGS * spacing;
        let near_pole = near && (a2 != n || a2_rounded != n);
        Some(Self {
            p,
            u,
            digits,
            rounding,
            near_pole,
        })
    }

    fn estimate(&self, margin: f64) -> Option<Estimate> {
        let y = self.p * self.u;
        let digits = self.digits - margin;
        // with the rounding doubled, as y can curve within the move, e.g. for an a1 near a
        // non-positive integer
        let err = 10.0_f64.powf(-digits) + 2.0 * self.rounding;
        // if accurate, U has the digits of `PFQ_TOL` plus the margin, and the rounding moves
        // y by at most `PFQ_TOL / 2`; an infinite y overflows by more than `err`
        let accurate = digits >= -PFQ_TOL.log10()
            && if y.is_infinite() {
                (self.p / (1.0 + err) * self.u).is_infinite()
            } else {
                self.rounding <= PFQ_TOL / 2.0
            };
        let err = if accurate {
            PFQ_TOL
        } else if self.near_pole {
            return None;
        } else if y.is_infinite() {
            f64::NAN
        } else {
            err.max(PFQ_INACCURATE)
        };
        // with the spacing of a subnormal y
        Some(Estimate::rounded(y, err))
    }
}

/// 2F0 for `z < 0` by Miller's algorithm (DLMF 13.29.5 to 13.29.7), and an estimate of its
/// absolute error, or `None` if it would take more than `PFQ_DD_MAXITER` steps, or if the result
/// of a shifted run is subnormal
///
/// With `x = -1/z`, the `v_k = (a1)_k (a2)_k / k! U(a1 + k, a1 - a2 + 1, x)` are the minimal
/// solution of
///
/// `(a1 + k - 1) (a2 + k - 1) v_(k-1) = k ((a1 + a2 + 2 k - 1 + x) v_k - (k + 1) v_(k+1))`
///
/// and sum to `x^(-a1)`, so that 2F0 is `v_0 / sum_k v_k` (DLMF 13.6.21): without a gamma
/// function, a power, or the rounded `x` and `a1 - a2 + 1`. It is in double-double, with a
/// first-order estimate of its rounding errors. The estimate is not a bound: it propagates the
/// errors of the steps with one choice of signs.
///
/// As the sum cancels unless `a1, a2 > 0`, it is that of `a1 + shift` and `a2 + shift`, and the
/// `f_k = 2F0(a1 + k, a2 + k;; z)` then go down from `k = shift` by (DLMF 13.3.7)
///
/// `f_(k-1) = (1 - (a1 + a2 + 2 k - 1) z) f_k - (a1 + k) (a2 + k) z^2 f_(k+1)`
///
/// with the rounding errors of its steps and of `f_shift` and `f_(shift+1)`, each times the
/// sensitivity of `f_0` to it.
fn hyp2f0_miller(a1: f64, a2: f64, z: f64) -> Option<(f64, f64)> {
    // the least for a1 + shift, a2 + shift > 0
    let shift = ((-a1).max(-a2).floor() + 1.0).max(0.0);
    let x = Dd::new(-1.0) / z;
    let sum = x + a1 + a2 + (2.0 * shift - 1.0);
    // v_0 / sum_k v_k from v_n = 1 and v_(n+1) = 0, and its relative rounding and truncation
    // errors (with the rounding of the quotient and of `y`)
    let ratio = |n: f64| {
        let (mut v, mut e) = ([Dd::ONE, Dd::ZERO], [0.0; 2]);
        let (mut total, mut e_total) = (v[0], 0.0);
        for k in (1..=n as u32).rev().map(f64::from) {
            let (p, q) = ((sum + 2.0 * k) * v[0], v[1] * (k + 1.0));
            let den = (Dd::new(a1) + (shift + k - 1.0)) * (Dd::new(a2) + (shift + k - 1.0)) / k;
            // the rounding error of this step, with the sign of the errors so far
            let rounding = Dd::EPSILON * (p.hi.abs() + q.hi.abs());
            let propagated = (sum.hi + 2.0 * k) * e[0] - (k + 1.0) * e[1];
            e = [(propagated + rounding.copysign(propagated)) / den.hi, e[0]];
            v = [(p - q) / den, v[0]];
            total = total + v[0];
            e_total += e[0].abs();
        }
        // f_shift and f_(shift+1) times sum_k v_k (DLMF 13.6.21)
        let den = (Dd::new(a1) + shift) * (Dd::new(a2) + shift) * -z;
        let mut f = [v[0], v[1] / den];
        // and their rounding errors
        let e_f1 = (e[1] / den.hi).abs() + Dd::EPSILON * f[1].hi.abs();
        let mut steps = Vec::with_capacity(shift as usize);
        for k in (1..=shift as u32).rev().map(f64::from) {
            let c = (Dd::new(a1) + a2 + (2.0 * k - 1.0)) * -z + 1.0;
            let d = (Dd::new(a1) + k) * (Dd::new(a2) + k) * z * z;
            let (p, q) = (c * f[0], d * f[1]);
            steps.push((c.hi, d.hi, Dd::EPSILON * (p.hi.abs() + q.hi.abs())));
            f = [p - q, f[0]];
        }
        // the sensitivities of f_0 to f_(k-1) and f_k, up from k = 1
        let (mut w, mut e_f) = ([1.0, 0.0], 0.0);
        for (c, d, rounding) in steps.into_iter().rev() {
            e_f += (w[0] * rounding).abs();
            w = [w[0] * c + w[1], -w[0] * d];
        }
        e_f += (w[0] * e[0]).abs() + (w[1] * e_f1).abs();
        let rounding = (e_f / f[0].hi).abs() + (e_total / total.hi).abs() + f64::EPSILON;
        // The sum lacks the v_k past n, and v_(n+1) = 0 adds some of the dominant solution to
        // those up to n. With the roots r < 1 < R of the recurrence at k = n, by which the two
        // solutions change per step, these are r / (1 - r) and r / (R - 1) times v_n, which is
        // the 1 it starts from over 1 - r / R. In all, that is r R / ((1 - r) (R - 1)), or
        //
        // `b1 b2 / (n x - (b1 - n) (b2 - n))` for `b_i = a_i + n - 1`
        //
        // and infinite where the v_k don't decay yet.
        let [b1, b2] = [a1, a2].map(|a| a + (shift + n - 1.0));
        let decay = (n * x.hi - (b1 - n) * (b2 - n)).max(0.0);
        ((f[0] / total).hi, rounding, b1 * b2 / (decay * total.hi.abs()))
    };
    // A heuristic: v_k is about k^((a1 + a2) / 2 - 1) K_(a1 - a2)(2 sqrt(k x)) for k >> a1, a2
    // (DLMF 13.8.11), which decays by exp(-m) once 2 sqrt(k x) is |a1 - a2| + m, but for the power
    // of k. So the steps double while the truncation error exceeds the rounding errors. There is
    // no result otherwise, as its estimate is too small where it is large (the roots are those
    // of step n only).
    let (m, max) = (-f64::EPSILON.ln(), f64::from(PFQ_DD_MAXITER));
    let n = (m + ((a1 - a2).abs() + m).powi(2) / (4.0 * x.hi)).ceil();
    successors((shift <= max).then_some(n.min(max)), |&n| (n < max).then(|| (2.0 * n).min(max)))
        .map(ratio)
        .find(|&(_, rounding, truncation)| truncation <= rounding)
        // (a subnormal `y` of the shifted run has lost digits that its estimate misses)
        .filter(|&(y, ..)| shift == 0.0 || y.is_normal())
        .map(|(y, rounding, truncation)| (y, (rounding + truncation) * y.abs()))
}

/// Borel sum of `2F0(a1, a2;; z)` for `z < 0` (DLMF 13.6.21)
///
/// Tries `x^a1 U(a1, a1 - a2 + 1, x)` with `CHGU_DIGITS_MARGIN`, the truncated asymptotic
/// series, the swapped `U` orientation (DLMF 13.2.40), double-double Weniger and double-double
/// residues (`pfq_borel_dd`), Miller's algorithm (`hyp2f0_miller`) where the residues cancel,
/// e.g. for `a1, a2 < 0`, and finally the cached `U` evaluations at `margin = 0.0`.
/// For large arguments and nearly terminating series, prefer the residues.
pub(super) fn hyp2f0_borel(a1: f64, a2: f64, z: f64) -> Option<Estimate> {
    // the smaller |a_i| first, as x^a1 then overflows less often
    let (a1, a2) = if a1.abs() <= a2.abs() {
        (a1, a2)
    } else {
        (a2, a1)
    };
    let x = -1.0 / z;
    // the rounding error of x (as x z is near -1, x z + 1 is exact)
    let (xz, e) = two_prod(x, z);
    let dx = -((xz + 1.0) + e) / z;
    let raw1 = RawBorel::eval(a1, a2, x, dx);
    let raw2 = OnceCell::new();
    let eval2 = || {
        (a1 != a2)
            .then(|| raw2.get_or_init(|| RawBorel::eval(a2, a1, x, dx)).as_ref())
            .flatten()
    };
    let a = [a1, a2];
    // keep the 2-factor asymptotic terms split in `(mantissa, exponent)` so an underflowing or
    // subnormal term (e.g. tiny `a1 * a2 * z`) stays a conservative bound rather than ending the
    // series before `t_1`
    let (zm, ze) = frexp_f64(z);
    let series_terms = || {
        (0_u32..)
            .map(f64::from)
            .scan((1.0_f64, 0_i64), move |(m, e), n| {
                let t = ldexp_i64(*m, *e);
                let bound = if *m == 0.0 {
                    0.0
                } else if t.is_finite() {
                    t.abs().max(f64::MIN_POSITIVE)
                } else {
                    f64::INFINITY
                };
                let ((m1, e1), (m2, e2)) = (frexp_f64(a1 + n), frexp_f64(a2 + n));
                let (nm, ne) = frexp_f64(*m * m1 * m2 * zm / (n + 1.0));
                *m = nm;
                *e += i64::from(e1) + i64::from(e2) + i64::from(ze) + i64::from(ne);
                Some((t, bound))
            })
    };
    let settled = (-a1).max(-a2).max(0.0).ceil();
    // with x capped, as sqrt(2 pi x) in the error estimate could overflow; it stays above
    // 1e8, while past `settled` the terms alternate, so the first omitted one bounds the
    // remainder
    let scale = x.min(f64::EPSILON.recip());
    let series =
        || Some(Estimate::absolute(asymptotic_sum(series_terms(), settled, scale, |_| 0.0)));
    let residues = || pfq_borel_dd(&a, &[], z);
    let weniger_dd = || pfq_weniger_dd(&a, &[], z, settled);
    raw1.as_ref()
        .and_then(|r| r.estimate(CHGU_DIGITS_MARGIN))
        .or_better(series)
        .or_better(|| eval2()?.estimate(CHGU_DIGITS_MARGIN))
        .or_better(|| {
            // When |z| >= 1, a nearly vanishing Pochhammer factor can stall the delta
            // approximants even after the terms start alternating. Small steps alone do
            // not establish accuracy there, so keep Weniger as a nonstrict estimate.
            let near_integer = |v: f64| (v - v.round().min(0.0)).abs() <= SQRT_EPSILON;
            if x <= 1.0 && !a.into_iter().any(is_nonpos_int) && a.into_iter().any(near_integer) {
                residues().or_better(|| {
                    let mut e = weniger_dd()?;
                    e.err = e.err.max(PFQ_INACCURATE);
                    Some(e)
                })
            } else {
                weniger_dd().or_better(residues)
            }
        })
        .or_better(|| hyp2f0_miller(a1, a2, z).map(Estimate::absolute))
        .or_better(|| raw1.as_ref()?.estimate(0.0))
        .or_better(|| eval2()?.estimate(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hyp2f0_nearly_terminating_weniger_plateau() {
        // Arb at the exact binary64 inputs. The small-|z| control still accepts Weniger
        // when the residues cannot establish accuracy.
        for (a, z, expected) in [
            (
                [-2.552_159_512_381_629e-13, -58.999_999_998_911_83],
                -126.750_930_337_152_8,
                -7.172_705_902_498_113_5e189,
            ),
            (
                [-7.458_722_797_026_388e-12, -12.999_999_999_997_922],
                -1_363.858_526_276_311_5,
                -2.019_982_460_778_705_2e38,
            ),
            (
                [-90.000_000_000_000_53, -18.599_116_159_816_333],
                -0.002_618_368_943_112_718_4,
                0.005_619_573_644_285_804_5,
            ),
        ] {
            let y = crate::hyppfq(a, [], z).unwrap();
            assert!((y / expected - 1.0).abs() <= PFQ_TOL, "{y:e} != {expected:e}");
        }
    }

    #[test]
    fn test_hyp2f0_miller_err() {
        // the estimate covers the rounding of the result, that of the recurrence in a1 and a2
        // (1.7e-13), and that of twice the steps, up to `PFQ_DD_MAXITER` (Arb)
        for (a1, a2, z, expected) in [
            (-36.5, -51.8, -0.026, -6.120_387_741_382_293e-9),
            (
                -50.999_999_999_699_924,
                -104.536_535_990_336_39,
                -0.056_613_881_920_597_09,
                2.794_502_416_895_517_4e22,
            ),
            (7.0, 9.0, -1.8, 6.911_864_984_250_697e-8),
        ] {
            let (y, err) = hyp2f0_miller(a1, a2, z).unwrap();
            assert!(err >= f64::EPSILON * y.abs() && (y - expected).abs() <= err, "{y:e} {err:e}");
        }
        // and there is no result where those steps leave a truncation error (1.7e-11), or where
        // the v_k don't decay within them
        assert_eq!(hyp2f0_miller(7.6, 6.8, -3.3), None);
        assert_eq!(hyp2f0_miller(760.0, 280.0, -0.014), None);
    }

    #[test]
    fn test_pfq_borel_dd() {
        // the Borel sum for p > q + 1, whose series are in -1 / z (2F0; Arb) or in 1 / z (3F0;
        // mpmath meijerg)
        for (a, z, expected) in [
            (vec![0.3, 1.7], -2.5, 0.651_478_922_494_888_7),
            (vec![0.3, 1.7, 0.6], -2.5, 0.790_728_400_803_319_9),
        ] {
            let y = pfq_borel_dd(&a, &[], z).unwrap();
            assert!(y.err < 8.0 * f64::EPSILON, "{y:?}");
            assert!((y.y / expected - 1.0).abs() <= y.err, "{y:?}");
        }
    }

    #[test]
    fn test_raw_borel_exact_and_margin() {
        // exact dyadic parameters (db == 0, dx == 0) skip perturbed evaluation with zero rounding
        let raw = RawBorel::eval(0.5, 1.5, 4.0, 0.0).unwrap();
        assert_eq!(raw.rounding, 0.0);
        assert!(raw.estimate(CHGU_DIGITS_MARGIN).unwrap().is_accurate());
        // 13-digit cached evaluation is inaccurate (still trusted) at margin 1, accurate at margin 0
        let raw13 = RawBorel {
            p: 1.0,
            u: 2.0,
            digits: 13.0,
            rounding: 0.0,
            near_pole: false,
        };
        let strict = raw13.estimate(CHGU_DIGITS_MARGIN).unwrap();
        assert!(!strict.is_accurate() && strict.trusted().is_some());
        assert!(raw13.estimate(0.0).unwrap().is_accurate());
        // near-pole gate rejects inaccurate estimates while keeping accurate ones
        let raw_pole = RawBorel {
            near_pole: true,
            ..raw13
        };
        assert!(raw_pole.estimate(CHGU_DIGITS_MARGIN).is_none());
        assert!(raw_pole.estimate(0.0).unwrap().is_accurate());
        // overflow gate: overflows comfortably vs barely past f64::MAX
        let raw_over = RawBorel {
            p: 1e200,
            u: 1e200,
            digits: 13.0,
            rounding: 0.0,
            near_pole: false,
        };
        assert!(
            raw_over
                .estimate(CHGU_DIGITS_MARGIN)
                .unwrap()
                .trusted()
                .is_none()
        );
        assert_eq!(raw_over.estimate(0.0).unwrap().accurate(), Some(f64::INFINITY));
        let raw_edge = RawBorel {
            p: f64::MAX,
            u: 1.0 + 1e-14,
            digits: 13.0,
            rounding: 0.0,
            near_pole: false,
        };
        assert!(raw_edge.estimate(0.0).unwrap().trusted().is_none());
    }
}
