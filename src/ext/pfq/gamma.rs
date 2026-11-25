//! Stirling double-double log-gamma and scaled Gamma-product helpers for `pFq`.

use super::asymptotic::Alg;
use super::{Dd, Scaled, SeriesReal, ln_dd, ln_gamma};
use crate::xsf::cephes::dd_real::{LN_SQRT_2PI_DD, frexp_f64};

/// The `|x|` below which `ln_gamma` is by `lgam`, so that `mul_gammas` shifts `x` to Stirling's
/// formula
const PFQ_GAMMA_MAX: f64 = 13.0;

/// The coefficients `B_2k / (2k (2k - 1))` of Stirling's series (DLMF 5.11.1) in double-double,
/// for `k = 16, ..., 1`
#[rustfmt::skip]
const STIRLING_DD: [Dd; 16] = [
    Dd { hi: -15_238_221.539_407_415, lo: -8.767_745_224_906_25e-10 },
    Dd { hi: 691_472.268_851_313, lo: 2.558_529_630_515_8e-11 },
    Dd { hi: -36_108.771_253_724_99, lo: 5.897_583_353_514_365e-13 },
    Dd { hi: 2_193.103_333_333_333_5, lo: -1.333_925_562_600_294_8e-13 },
    Dd { hi: -156.848_284_626_002_03, lo: 9.391_823_141_715_389e-15 },
    Dd { hi: 13.402_864_044_168_393, lo: -6.154_114_101_993_966e-16 },
    Dd { hi: -1.392_432_216_905_901_1, lo: 1.583_705_698_923_030_3e-17 },
    Dd { hi: 0.179_644_372_368_830_57, lo: -6.401_600_482_710_946e-19 },
    Dd { hi: -0.029_550_653_594_771_242, lo: 4.861_760_957_508_855e-19 },
    Dd { hi: 0.006_410_256_410_256_41, lo: 2.224_004_456_380_521_7e-19 },
    Dd { hi: -0.001_917_526_917_526_917_6, lo: 1.067_570_277_687_247_5e-19 },
    Dd { hi: 0.000_841_750_841_750_841_7, lo: 3.687_017_488_923_769_4e-20 },
    Dd { hi: -0.000_595_238_095_238_095_3, lo: 5.369_382_187_547_26e-20 },
    Dd { hi: 0.000_793_650_793_650_793_7, lo: 6.883_823_317_368_282e-22 },
    Dd { hi: -0.002_777_777_777_777_778, lo: 1.060_108_790_874_715_4e-19 },
    Dd { hi: 0.083_333_333_333_333_33, lo: 4.625_929_269_271_485e-18 },
];

/// The term `k = 17` of Stirling's series at `PFQ_GAMMA_MAX`, which bounds the sum of the terms
/// after those of `STIRLING_DD` from there on
const STIRLING_DD_REST: f64 = 6.7e-29;

/// The `g < 0` down to which `mul_gammas_in` shifts `g` to `ln_gamma_stirling`, in about `|g|`
/// steps
const PFQ_GAMMA_DD_MIN: f64 = -1024.0;

/// `ln Gamma(g)` for `g >= PFQ_GAMMA_MAX` by Stirling's series (DLMF 5.11.1) in double-double,
/// its sign, and a bound on its absolute error
fn ln_gamma_stirling(g: Dd) -> (Dd, f64, f64) {
    let w = Dd::ONE / (g * g);
    let corr = STIRLING_DD.iter().fold(Dd::ZERO, |s, &c| s * w + c) / g;
    let v = (g - 0.5) * ln_dd(g) - g + LN_SQRT_2PI_DD + corr;
    (v, 1.0, STIRLING_DD_REST + Dd::EPSILON * v.hi.abs())
}

/// `y` times `v`, with the low part of the mantissa in the logarithm
fn mul_dd(mut y: Scaled, v: Dd) -> Scaled {
    let p = v * y.m;
    y.m = p.hi;
    y.mul_exp(Dd::new(p.lo / p.hi))
}

type LogGamma = (Dd, f64, f64);
pub(super) type LogGammaCache<'a> = &'a mut [Option<(Dd, LogGamma)>];

// The Gamma product still applies every sign, error, and scale update in its original order.
fn cached_ln_gamma(g: Dd, cache: LogGammaCache<'_>) -> Option<LogGamma> {
    for slot in cache {
        if let Some((key, value)) = *slot {
            if key.hi.to_bits() == g.hi.to_bits() && key.lo.to_bits() == g.lo.to_bits() {
                return Some(value);
            }
        } else {
            let value = ln_gamma(g)?;
            *slot = Some((g, value));
            return Some(value);
        }
    }
    ln_gamma(g)
}

/// `y` times `Gamma(num_0) ... / (Gamma(den_0) ...)`, with the error bounds of the factors
/// added, or `None` at a pole (or if the result is far out of range)
///
/// Each factor is by `ln_gamma`, which for `|g| < PFQ_GAMMA_MAX` is at `g + n >= PFQ_GAMMA_MAX`
/// instead, divided by the rising factorial `(g)_n` (DLMF 5.2.5). There it is Stirling's formula
/// in double-double, whose error is smaller than that of `lgam` at `g`. The `(g)_n` are multiplied
/// in double-double, and rounded once.
#[cfg(test)]
fn mul_gammas(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    y: Scaled,
) -> Option<Scaled> {
    mul_gammas_in(num, den, y, false, &mut [])
}

/// `mul_gammas` if not `dd`. Otherwise each factor is by `ln_gamma_stirling` at `g + n`, also for
/// `PFQ_GAMMA_DD_MIN < g <= -PFQ_GAMMA_MAX`, where `ln_gamma` reflects `g` in `f64`, and the
/// `(g)_n` are not rounded, so that the factors are as exact as double-double arithmetic.
fn mul_gammas_in(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    mut y: Scaled,
    dd: bool,
    cache: LogGammaCache<'_>,
) -> Option<Scaled> {
    // the product of the `(g)_n^-power`, as a mantissa and an exponent
    // (`rising_e` is an `f64`, which sums exactly, as many tiny factors can take it out of the
    // range of `i32`)
    let (mut rising, mut rising_e) = (Dd::ONE, 0.0);
    let gammas = num.into_iter().map(|g| (g, 1.0));
    for (mut g, power) in gammas.chain(den.into_iter().map(|g| (g, -1.0))) {
        let min = if dd { PFQ_GAMMA_DD_MIN } else { -PFQ_GAMMA_MAX };
        if min < g.hi && g.hi < PFQ_GAMMA_MAX {
            // (not the ceiling, which can give `g + n < PFQ_GAMMA_MAX` as the difference rounds)
            let n = (PFQ_GAMMA_MAX - g.hi).floor() + 1.0;
            let ks = (0_u32..).map(f64::from).take_while(|&k| k < n);
            // `(g)_n` as a product of mantissas and an exponent, so that it is normal also for
            // a tiny factor (the product is at least `2^-n`)
            let (r, r_e) = ks.fold((Dd::ONE, 0), |(p, e), k| {
                let (f, j) = (g + k).frexp();
                (p * f, e + j)
            });
            // (else a factor is 0: `g` is a pole, where `ln_gamma` is `None`)
            if r != 0.0 {
                let (q, s) = if power > 0.0 {
                    (rising / r, -r_e)
                } else {
                    (rising * r, r_e)
                };
                let (m, e) = q.frexp();
                (rising, rising_e, g) = (m, rising_e + f64::from(e + s), g + n);
                if dd {
                    // (the rounding of each factor of `(g)_n`)
                    y.err += n * Dd::EPSILON;
                }
            }
        }
        let (v, sign, v_err) = if dd && g.hi >= PFQ_GAMMA_MAX {
            ln_gamma_stirling(g)
        } else {
            cached_ln_gamma(g, cache)?
        };
        (y.m, y.err) = (sign * y.m, y.err + v_err);
        y = y.mul_exp(v * power);
    }
    if rising != 1.0 {
        // (one eps bounds the rounding of `rising` to `rising.hi` and of the product of the
        // mantissas, which can't overflow or underflow; `None` for an exponent beyond the range
        // that `fold_exp` takes, which needs about a million tiny Gamma arguments)
        let ((m, e), (r, k)) = (frexp_f64(y.m), frexp_f64(rising.hi));
        let e = f64::from(y.e + e + k) + rising_e;
        y.e = (e.abs() <= 1e9).then_some(e as i32)?;
        if dd {
            (y.m, y.err) = (m, y.err + Dd::EPSILON);
            y = mul_dd(y, rising.ldexp(-k));
        } else {
            (y.m, y.err) = (m * r, y.err + f64::EPSILON);
        }
    }
    Some(y)
}

/// `exp(l) sum` times `Gamma(num_0) ... / (Gamma(den_0) ...)`, where `err` is the absolute
/// error of `sum`
pub(super) fn gamma_series(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    (sum, err): (f64, f64),
    l: Dd,
    cache: LogGammaCache<'_>,
) -> Option<Scaled> {
    mul_gammas_in(num, den, Scaled::new(sum, 0, err / sum.abs())?.with_l(l), false, cache)
}

/// Appends `exp(l) sum` times `Gamma(num_0) ... / (Gamma(den_0) ...)` to `terms` for a term of
/// `pfq_alg`, skipping it if `(sum, err) == (0, 0)`, or returns `None` if `err` is not finite or
/// a Gamma factor fails.
///
/// If `sum` vanishes (or `err / |sum|` overflows) with a finite `err > 0`, pushes a zero term
/// carrying the scaled absolute error (or the scaled bound itself when `sizes` is set).
pub(super) fn gamma_alg_term(
    num: impl IntoIterator<Item = Dd>,
    den: impl IntoIterator<Item = Dd>,
    (sum, err): (Dd, f64),
    l: Dd,
    mode: Alg,
    terms: &mut Vec<Scaled>,
    cache: LogGammaCache<'_>,
) -> Option<()> {
    let dd = mode == Alg::Dd;
    let sizes = mode == Alg::Sizes;
    err.is_finite().then_some(())?;
    let rel = err / sum.hi.abs();
    let error_only = sum.hi == 0.0 || !rel.is_finite();
    let base = if error_only {
        let abs_err = err + sum.hi.abs();
        if abs_err == 0.0 {
            return Some(());
        }
        Scaled::new(abs_err, 0, 0.0)?.with_l(l)
    } else if dd {
        mul_dd(Scaled::new(1.0, 0, rel)?.with_l(l), sum)
    } else {
        Scaled::new(sum.hi, 0, rel)?.with_l(l)
    };
    let mut term = mul_gammas_in(num, den, base, dd, cache)?;
    if sizes {
        term.m = term.m.abs();
    } else if error_only {
        term = term.into_zero_err();
    }
    terms.push(term);
    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::pow2;

    #[test]
    fn test_mul_gammas() {
        // Gamma(g) Gamma(180.5) / (Gamma(178.25) Gamma(1.5)) for g = -3 + 2^-40 + lo, where the
        // large arguments are by `ln_gamma`, and where (lo 2^40)^2 is negligible or not (mpmath)
        let one = Scaled::new(1.0, 0, 0.0).unwrap();
        let (num, den) = ([Dd::new(180.5)], [178.25, 1.5].map(Dd::new));
        let check = |y: Option<Scaled>, expected: f64, tol: f64| {
            let (y, err) = y.unwrap().to_f64();
            assert!(err <= tol && (y / expected - 1.0).abs() <= err, "{y} != {expected}");
        };
        for (e, expected, tol) in
            [(-75, -2.419_553_175_200_6e16, 1e-14), (-60, -2.419_550_867_807_498_4e16, 1e-12)]
        {
            let g = Dd::new(-3.0 + pow2(-40)) + pow2(e);
            check(mul_gammas(num.into_iter().chain([g]), den, one), expected, tol);
        }
        // equal factors, whose errors are equal too
        let g = Dd::new(-8.000_005_264_759_276);
        check(mul_gammas([g; 3], [], one), -104.541_370_932_890_13, 1e-14);
        // the mantissa of `y` times `Gamma(g)` overflows
        let (g, big) = (Dd::new(-3.0 + pow2(-40)), Scaled::new(1e300, 0, 0.0));
        check(mul_gammas([g], [g], big.unwrap()), 1e300, 1e-14);
        // `PFQ_GAMMA_MAX - g` rounds to an integer, and `(g)_n` is below the normal range
        check(
            mul_gammas([Dd::new(-8.000_000_000_000_002)], [], one),
            -13_962_052_416.203_12,
            1e-15,
        );
        check(mul_gammas([Dd::new(1e-290)], [], one), 9.999_999_999_999_999e289, 1e-15);
        check(mul_gammas([], [Dd::new(-1e-290)], one), -1e-290, 1e-15);
        // a pole
        assert!(mul_gammas([Dd::new(-2.0)], den, one).is_none());
        assert!(mul_gammas(num, [Dd::new(-2.0)], one).is_none());
        // in double-double, below `-PFQ_GAMMA_MAX`: Gamma(-30.5) / Gamma(-28.5)
        let y = mul_gammas_in([Dd::new(-30.5)], [Dd::new(-28.5)], one, true, &mut []);
        check(y, 1.0 / (30.5 * 29.5), 1e-15);
    }

    #[test]
    fn test_ln_gamma_stirling() {
        // where its error is the largest: Gamma(13) = 12!
        let (v, sign, err) = ln_gamma_stirling(Dd::new(PFQ_GAMMA_MAX));
        assert!(sign == 1.0 && err <= 1e-28);
        assert!((v - ln_dd(479_001_600.0)).hi.abs() <= err);
    }
}
