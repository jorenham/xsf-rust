//! Divergent `p > q + 1` Borel summation and `2F0` evaluation via `hypu`.

use core::cell::OnceCell;
use core::iter::once;

use super::asymptotic::{Alg, asymptotic_sum, pfq_alg, pfq_alg_dd};
use super::series::{pfq_terms, pfq_weniger, weniger_drift};
use super::{
    Dd, Estimate, OrBetter, PFQ_INACCURATE, PFQ_TOL, SeriesReal, add_exp, ldexp_i64, pfq_settled,
};
use crate::utils::{is_nonpos_int, pow2};
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
    let terms = || once(1.0).chain(pfq_terms(a, b, z));
    let settled = pfq_settled(a, b).ceil();
    let series = || Some(Estimate::absolute(asymptotic_sum(terms(), settled, scale, |_| 0.0)));
    // also with the rounding errors of the terms, as for 2F0
    let weniger_dd =
        || Some(weniger_drift(pfq_weniger::<Dd>(a, b, z)?, terms(), settled, Dd::EPSILON));
    let residues = || Some(Estimate::scaled(add_exp(pfq_alg(a, b, z, Alg::Value)?)?, None));
    series()
        .or_better(residues)
        .or_better(weniger_dd)
        .or_better(|| pfq_borel_dd(a, b, z))
}

// The Laplace kernel (1 + |z| X Y)^(-a) lies in [0,1] and is symmetric in X,Y.
#[cold]
pub(super) fn hyp3f0_near_coincident(mut a: [f64; 3], z: f64) -> Option<Estimate> {
    a.sort_unstable_by(f64::total_cmp);
    let [lo, mid, hi] = a;
    let delta = mid - lo;
    // Sterbenz's lemma gives exact differences when both ratios are at most two.
    if !(lo > 0.0 && delta > 0.0 && hi <= 2.0 * mid && mid <= 2.0 * lo && hi - mid == delta) {
        return None;
    }
    // The symmetrized Gamma densities have likelihood ratio cosh(delta*ln(X/Y))/M,
    // with M = Gamma(lo)*Gamma(hi)/Gamma(mid)^2. Their TV distance is at most
    // 1 - 1/M <= ln(M) <= delta² (1/lo² + 1/lo), using the trigamma series.
    let ratio = delta / lo;
    let bound = 2.0 * ratio * (ratio + delta); // allow for rounding in the bound
    if bound > PFQ_TOL {
        return None;
    }
    let centered = pfq_borel(&[mid; 3], &[], z)?;
    let err = centered.err + bound / centered.y;
    (centered.y.is_normal() && centered.y > 0.0 && err <= PFQ_TOL)
        .then_some(Estimate::new(centered.y, err))
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

/// Borel sum of `2F0(a1, a2;; z)` for `z < 0` (DLMF 13.6.21)
///
/// Tries `x^a1 U(a1, a1 - a2 + 1, x)` with `CHGU_DIGITS_MARGIN`, the truncated asymptotic
/// series, the swapped `U` orientation (DLMF 13.2.40), double-double Weniger, double-double
/// residues (`pfq_borel_dd`), and finally the cached `U` evaluations at `margin = 0.0`.
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
    let terms = || once(1.0).chain(pfq_terms(&a, &[], z));
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
    let weniger_dd =
        || Some(weniger_drift(pfq_weniger::<Dd>(&a, &[], z)?, terms(), settled, Dd::EPSILON));
    raw1.as_ref()
        .and_then(|r| r.estimate(CHGU_DIGITS_MARGIN))
        .or_better(series)
        .or_better(|| eval2()?.estimate(CHGU_DIGITS_MARGIN))
        .or_better(weniger_dd)
        .or_better(|| pfq_borel_dd(&a, &[], z))
        .or_better(|| raw1.as_ref()?.estimate(0.0))
        .or_better(|| eval2()?.estimate(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

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
