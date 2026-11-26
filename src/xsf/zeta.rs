//! The complex `riemann_zeta` is translated into pure Rust from `xsf/zeta.h` (scipy/xsf@3cf1b75).
//!
//! Complex riemann-zeta function implementation based on Python implementation
//! written by Matt Haberland (@mdhaber) in:
//! <https://colab.research.google.com/drive/1zMDSAJlXCLRqMMtJ0e9nDGjQ8iZCnmn5?usp=sharing>
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_8, LN_2, PI, TAU};

use crate::utils::{LN_PI, LN_TAU, neg_if_odd, neg_one_pow};
use crate::xsf::c64::{C64, ViaC64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, FRAC_2_PI_DD, LN_SQRT_2PI_DD};
use crate::xsf::gamma::cgamma;
use crate::xsf::log::clog1p;
use crate::xsf::loggamma::cloggamma;
use crate::xsf::tools::kahan_step;

/// `std::pow(x, z)` for real `x > 0` and complex `z`, as in libstdc++
fn rpow(x: f64, z: C64) -> C64 {
    debug_assert!(x > 0.0);
    // std::polar(pow(x, z.real()), z.imag() * log(x))
    let rho = x.powf(z.re);
    let theta = z.im * x.ln();
    rho * C64::cis(theta)
}

/// Log of absolute value of expansion coefficients for Euler-Maclaurin
/// summation formula. log(|B2k / (2k)!|)
///
/// See <https://en.wikipedia.org/wiki/Riemann_zeta_function#Numerical_algorithms>
///
/// Generated with the script
///
/// ```python
/// import numpy as np
/// from mpmath import mp
///
/// mp.dps = 10000
///
/// results = []
/// for k in range(51):
///     results.append(
///         float(
///             mp.log(abs(mp.bernoulli(2*k)/(mp.factorial(2*k))))
///         )
///     )
/// ```
const ZETA_EM_LOG_ABS_COEFF_LOOKUP: [f64; 51] = [
    0.0,
    -2.484_906_649_788_000_4,
    -6.579_251_212_010_101,
    -1.031_692_083_029_347e1,
    -1.400_580_028_440_740_5e1,
    -1.768_462_940_266_784e1,
    -2.136_113_156_007_322_2e1,
    -2.503_707_050_291_142_3e1,
    -2.871_287_059_984_694_8e1,
    -3.238_863_619_752_229_5e1,
    -3.606_439_319_366_539e1,
    -3.974_014_804_199_518_4e1,
    -4.341_590_235_365_616e1,
    -4.709_165_653_118_148_5e1,
    -5.076_741_067_517_639e1,
    -5.444_316_481_078_909e1,
    -5.811_891_894_430_628e1,
    -6.179_467_307_729_959e1,
    -6.547_042_721_016_194e1,
    -6.914_618_134_299_154e1,
    -7.282_193_547_581_296e1,
    -7.649_768_960_863_234e1,
    -8.017_344_374_145_12e1,
    -8.384_919_787_426_993e1,
    -8.752_495_200_708_863e1,
    -9.120_070_613_990_733e1,
    -9.487_646_027_272_602e1,
    -9.855_221_440_554_472e1,
    -1.022_279_685_383_634e2,
    -1.059_037_226_711_821e2,
    -1.095_794_768_040_007_8e2,
    -1.132_552_309_368_194_7e2,
    -1.169_309_850_696_381_7e2,
    -1.206_067_392_024_568_5e2,
    -1.242_824_933_352_755_5e2,
    -1.279_582_474_680_942_4e2,
    -1.316_340_016_009_129_4e2,
    -1.353_097_557_337_316_3e2,
    -1.389_855_098_665_503e2,
    -1.426_612_639_993_69e2,
    -1.463_370_181_321_877e2,
    -1.500_127_722_650_064e2,
    -1.536_885_263_978_251e2,
    -1.573_642_805_306_437_5e2,
    -1.610_400_346_634_624_5e2,
    -1.647_157_887_962_811_5e2,
    -1.683_915_429_290_998_4e2,
    -1.720_672_970_619_185_4e2,
    -1.757_430_511_947_372_3e2,
    -1.794_188_053_275_559e2,
    -1.830_945_594_603_746e2,
];

/// Log of the absolute value of the expansion coefficients for Euler-Maclaurin summation formula.
///
/// Unlike in xsf, the sign of B2n/(2n)! (negative for even n) is not a term pi*i of the log,
/// whose rounding would give the terms a spurious imaginary part (wrong near the real axis).
fn zeta_em_log_abs_coeff(n: usize) -> f64 {
    if n < 50 {
        ZETA_EM_LOG_ABS_COEFF_LOOKUP[n]
    } else {
        /* Asymptotic formula
         * Uses https://dlmf.nist.gov/24.11#E1 to approximate B_{2n} and
         * Stirling's approximation for (2n)!.
         */
        LN_2 - 2.0 * n as f64 * LN_TAU
    }
}

/// Compute riemann_zeta for complex input z using the Euler-Maclaurin formula.
/// Computation of individual terms in expansion are logarithmized to avoid
/// overflow. TODO: only logarithmize when necessary.
fn zeta_euler_maclaurin(z: C64) -> C64 {
    if z == C64::ONE {
        /* Return NaN at pole since value depends on how z approaches 1.0. */
        return C64::NAN;
    }
    let n = (z.im.abs() / 4.0).max(50.0);
    if n >= usize::MAX as f64 {
        // Unlike in xsf, NaN if the number of terms doesn't fit in a `size_t` (UB in C++)
        return C64::NAN;
    }
    let n = n as usize;
    let m = n;
    let mut result = C64::ZERO;
    for i in 1..n {
        let term = rpow(i as f64, -z);
        result += term;
        // When z.real() > 1, series converges and we can consider early termination
        if z.re > 1.0 && term.norm() / result.norm() <= f64::EPSILON {
            return result;
        }
    }
    let big_n = n as f64;
    let b = rpow(big_n, -z);
    // Unlike in xsf, b n is divided by z - 1, as n / (z - 1) overflows near z = 1, and this pole
    // term is added after the tail: near z = 1 it dwarfs the real part, which would then end the
    // tail too early.
    let pole = b * big_n / (z - 1.0);
    result += 0.5 * b;
    /* The terms of the Euler-Maclaurin
     * expansion below are T(k, n) = B2k/(2k)! * n^(1 - z - 2k) * z(z+1)...(z+2k-2).
     * We work with logarithms to avoid overflow in all cases at the expense of
     * some accuracy. At the start of iteration k:
     *     log_poch will equal log(z(z+1)...(z+2k-2))
     *     log_factor will equal log(n^(1 - z - 2k))
     * These are updated one extra time after the loop completes for use in the
     * Euler-Maclaurin error estimate.
     */
    let mut log_poch = z.ln();
    let mut log_factor = -(z + 1.0) * big_n.ln();
    for k in 1..=m {
        let term = neg_if_odd(k + 1, (zeta_em_log_abs_coeff(k) + log_factor + log_poch).exp());
        result += term;
        if term.norm() / result.norm() <= f64::EPSILON {
            return result + pole;
        }
        log_poch += (z + (2 * k - 1) as f64).ln() + (z + (2 * k) as f64).ln();
        log_factor -= 2.0 * big_n.ln();
    }
    result += pole;
    /* Euler-maclaurin absolute error estimate.
     * The error is bounded above by |(z + 2m + 1)/(z.real + 2m + 1) * T(m+1, n)|
     * See https://en.wikipedia.org/wiki/Riemann_zeta_function#Numerical_algorithms
     */
    let mut error = (zeta_em_log_abs_coeff(m + 1) + log_factor + log_poch)
        .exp()
        .norm();
    error *= ((z + 2.0 * m as f64 + 1.0) / (z.re + 2.0 * m as f64 + 1.0)).norm();
    // convert to relative error estimate
    error /= result.norm();
    if error > 1e-1 {
        /* If error estimate predicts we don't even get 1 digit of precision, return NaN
         * and signal no result */
        return C64::NAN;
    }
    result
}

/// Lookup table of coefficients for Algorithm 2 from Borwein 1995
/// Borwein, Peter B.. "An efficient algorithm for the Riemann zeta function." (1995).
///
/// Stores coefficients as dk / dn, where dn is the final coefficient.
///
/// Generated with the Python script:
///
/// ```python
/// import numpy as np
/// import math
/// from mpmath import mp
///
/// mp.dps = 1000
/// n = 50
///
/// coeffs = []
/// S = mp.zero
/// for i in range(n + 1):
///     num = math.factorial(n + i - 1) * 4**i
///     den = math.factorial(n - i) * math.factorial(2*i)
///     S += mp.mpf(num) / mp.mpf(den)
///     coeffs.append(S*n)
///
/// dn = coeffs[-1]
/// coeffs = [float(dk/dn) for dk in coeffs[:-1]]
/// coeffs = np.asarray(coeffs)
/// ```
const ZETA_BORWEIN_COEFF: [f64; 50] = [
    1.055_507_836_138_287_8e-38,
    5.278_594_688_527_578e-35,
    4.401_468_732_204_496_3e-32,
    1.467_453_546_497_519e-29,
    2.617_862_196_688_831e-27,
    2.900_097_799_958_025e-25,
    2.184_440_361_492_933e-23,
    1.189_097_731_291_329_6e-21,
    4.887_139_616_687_227_6e-20,
    1.567_225_369_880_273_4e-18,
    4.022_931_234_264_572e-17,
    8.435_973_533_351_745e-16,
    1.469_296_379_914_116e-14,
    2.154_874_705_457_190_2e-13,
    2.691_953_053_753_512_4e-12,
    2.892_540_916_276_848_4e-11,
    2.695_750_569_369_985_6e-10,
    2.194_772_239_130_839e-9,
    1.570_782_293_700_57e-8,
    9.936_187_220_749_133e-8,
    5.581_721_578_217_702e-7,
    2.796_271_112_037_765e-6,
    1.253_886_254_275_813e-5,
    5.049_261_002_939_051e-5,
    1.831_288_445_970_366_6e-4,
    5.997_690_328_552_426e-4,
    1.778_050_108_246_096_8e-3,
    4.781_802_283_665_968e-3,
    1.169_043_228_713_167_1e-2,
    2.603_430_292_953_596_4e-2,
    5.292_298_247_275_485_6e-2,
    9.842_471_411_648_858e-2,
    1.678_961_079_654_034_4e-1,
    2.635_042_919_476_862_6e-1,
    3.819_442_810_600_665e-1,
    5.137_731_385_068_898e-1,
    6.452_925_385_418_66e-1,
    7.625_449_322_050_362e-1,
    8.556_063_057_019_102e-1,
    9.210_560_628_865_25e-1,
    9.616_100_580_028_147e-1,
    9.835_905_431_606_953e-1,
    9.939_187_229_336_753e-1,
    9.980_782_525_166_648e-1,
    9.994_930_141_459_889e-1,
    9.998_914_788_445_85e-1,
    9.999_819_091_990_203e-1,
    9.999_977_981_288_319e-1,
    9.999_998_260_580_315e-1,
    9.999_999_933_099_243e-1,
];

/// Compute riemann_zeta for complex input z using Algorithm 2 from Borwein 1995.
///
/// `eta_factor` is `1 - 2^(1 - z)`.
fn zeta_borwein(z: C64, eta_factor: C64) -> C64 {
    let mut result = C64::ZERO;
    // Sum in reverse order because smaller terms come later.
    for k in (0..50_i32).rev() {
        let sign = neg_one_pow(k);
        let den = rpow(f64::from(k + 1), z);
        let term = sign * (ZETA_BORWEIN_COEFF[k as usize] - 1.0) / den;
        result += term;
    }
    -result / eta_factor
}

/// Taylor coefficients `c_2k`, up to `k = 29`, of the even entire function
/// `F(p) = (exp(pi i (p^2/2 + 3/8)) - i sqrt(2) cos(pi p/2)) / (2 cos(pi p))` of the
/// Riemann-Siegel correction.
///
/// Generated with the python-flint script
///
/// ```python
/// from flint import acb, acb_series, arb, ctx
///
/// ctx.prec, ctx.cap = 256, 60
/// z, i = acb_series([0, 1]), acb(0, 1)
/// F = ((i * arb.pi() * (z * z / 2 + acb(3) / 8)).exp() - i * arb(2).sqrt() * (z / 2).cos_pi()) / (
///     2 * z.cos_pi()
/// )
/// coeffs = [complex(c) for c in F.coeffs()[::2]]
/// ```
const RS_F_COEFFS: [C64; 30] = [
    c64(0.191_341_716_182_544_9, -0.245_167_014_930_904_15),
    c64(0.218_620_234_038_760_21, -0.036_933_834_884_962_956),
    c64(0.066_188_287_740_171_76, 0.063_534_393_856_146_02),
    c64(-0.006_802_513_023_837_094, 0.027_223_912_663_570_066),
    c64(-0.006_783_810_985_051_790_5, 0.001_385_760_877_106_652),
    c64(-0.000_811_862_661_572_232_7, -0.001_189_449_446_101_378),
    c64(0.000_148_526_768_666_898_45, -0.000_212_698_201_928_933_23),
    c64(3.971_650_439_760_735e-5, 1.117_132_740_199_015_1e-5),
    c64(2.327_806_230_725_225_2e-7, 5.872_858_398_652_07e-6),
    c64(-7.163_625_815_477_553e-7, 2.498_212_552_923_518e-7),
    c64(-5.177_423_556_156_473e-8, -7.308_700_305_101_552e-8),
    c64(6.178_963_541_930_869e-9, -7.536_791_448_164_02e-9),
    c64(8.940_541_928_977_453e-10, 4.104_425_797_331_228_4e-10),
    c64(-1.695_707_194_963_518e-11, 9.106_559_550_294_084e-11),
    c64(-8.163_316_951_282_953e-12, 4.348_099_095_249_618_7e-13),
    c64(-1.892_554_659_270_610_3e-13, -6.520_913_261_540_13e-13),
    c64(4.663_711_629_600_862_5e-14, -2.574_698_839_194_822e-14),
    c64(2.610_921_507_989_068_5e-15, 2.978_335_196_086_287_4e-15),
    c64(-1.675_336_536_372_132e-16, 2.241_756_964_196_517e-16),
    c64(-1.706_213_261_405_863_2e-17, -7.993_661_378_773_457e-18),
    c64(2.875_601_670_716_199_6e-19, -1.176_894_351_646_754_5e-18),
    c64(7.447_650_681_605_753e-20, 3.424_141_569_957_982e-21),
    c64(6.282_686_358_510_708e-22, 4.354_432_465_678_006e-21),
    c64(-2.360_647_625_071_713e-22, 8.042_677_010_875_067e-23),
    c64(-6.634_534_681_519_81e-24, -1.187_404_814_328_495e-23),
    c64(5.526_719_997_560_709e-25, -4.533_873_522_499_421e-25),
    c64(2.749_823_188_763_732_5e-26, 2.362_213_241_963_689_5e-26),
    c64(-9.115_688_251_159_012e-28, 1.524_413_235_942_887e-27),
    c64(-7.844_701_868_860_44e-29, -3.053_767_696_325_618_7e-29),
    c64(7.919_817_544_119_006e-31, -3.781_590_707_554_057_6e-30),
];

/// Number of terms `C_k / a^k` of the Riemann-Siegel correction
const RS_TERMS: usize = 10;

/// Smallest `a = sqrt(|Im z| / (2 pi))` for the Riemann-Siegel formula. The cutoff is empirical:
/// on a grid of `|p| <= 1` and the dispatched `Re z`, the first omitted correction terms, weighted
/// as `a^(-Re z - 10) (|C_10(Re z, p)| + |C_10(1 - Re z, p)|)`, stay below eps.
const RS_MIN_A: f64 = 10.0;

/// Largest `|Im z|` with `0 <= Re z < 2.5` for which zeta is computed (NaN above)
const MAX_IM: f64 = 1e9;

/// Cost of a Dirichlet term (a power, a double-double logarithm and a sine and cosine) in steps
/// of `RS_CORRECTION_TERMS`. Calibrated by timing at large `|Im z|`; at small `|Im z|` the
/// Dirichlet sum stays cheaper somewhat beyond the resulting cutoff.
const DIRICHLET_TERM_STEPS: usize = 16;

/// Cost of the Riemann-Siegel correction in Dirichlet terms, counting `2 RS_F_COEFFS.len()` steps
/// for each of its `3 RS_TERMS - 2` derivatives of `F`.
const RS_CORRECTION_TERMS: u32 =
    ((3 * RS_TERMS - 2) * 2 * RS_F_COEFFS.len() / DIRICHLET_TERM_STEPS) as u32;

/// `sum_k C_k(sigma, p) / a^k` of the Riemann-Siegel formula, with
/// `C_k(sigma, p) = sum_l d_kl(sigma) F^(3k - 2l)(p) / (pi^(2k - l) (2i)^l)`, given the derivatives
/// `fp[m] = F^(m)(p)`.
fn rs_correction(sigma: f64, a: f64, fp: &[C64; 3 * RS_TERMS - 2]) -> C64 {
    // Row k of d_kl(sigma) by the recurrence of mpmath's `Rzeta_simul`, with d_kl at d[l + 2] so
    // that d_k,-2 = d_k,-1 = 0
    let mut d = [0.0; 3 * RS_TERMS / 2 + 3];
    d[2] = 1.0;
    let mut sum = fp[0];
    let mut scale = 1.0;
    for k in 1..RS_TERMS {
        let prev = d;
        // 1 / (a^k pi^(2k - l) (2i)^l) = (pi^2 a)^-k (-i pi / 2)^l
        scale /= PI * PI * a;
        let mut g = C64::from(scale);
        for l in 0..=3 * k / 2 {
            let m = 3 * k - 2 * l;
            d[l + 2] = if m == 0 {
                -(0..l)
                    .map(|r| {
                        let j = l - r;
                        neg_if_odd(j, d[r + 2])
                            * ((j + 1)..=2 * j).map(|i| i as f64).product::<f64>()
                    })
                    .sum::<f64>()
            } else {
                let m = m as f64;
                -(m + 1.0) * prev[l]
                    + prev[l + 2] / (4.0 * m)
                    + (1.0 - 2.0 * sigma) / (2.0 * m) * prev[l + 1]
            };
            sum += d[l + 2] * g * fp[m];
            g *= c64(0.0, -FRAC_PI_2);
        }
    }
    sum
}

/// `sum_{n <= n_max} n^-z` and `sum_{n <= n_max} n^(z - 1)`, compensated and from the smallest
/// terms.
fn dirichlet_sums(z: C64, n_max: u32) -> (C64, C64) {
    let (sigma, t) = (z.re, z.im);
    let (mut x, mut y) = (C64::ZERO, C64::ZERO);
    let (mut x_comp, mut y_comp) = (C64::ZERO, C64::ZERO);
    for n in (1..=n_max).rev().map(f64::from) {
        let w = n.powf(-sigma);
        // The f64 phase t ln(n) is off by about eps t ln(n), which this and the smaller terms
        // weigh by at most n^-sigma (1 + n / (sigma - 1)), so it is computed in double-double
        // where that exceeds eps. For sigma >= 1/2, the terms of the second sum weigh no more in
        // the Riemann-Siegel formula.
        let (sin, cos) = if sigma <= 1.0 || t * n.ln() * w * (1.0 + n / (sigma - 1.0)) > 1.0 {
            (Dd::new(n).log() * t).sin_cos()
        } else {
            (t * n.ln()).sin_cos()
        };
        kahan_step(&mut x, &mut x_comp, w * c64(cos, -sin));
        kahan_step(&mut y, &mut y_comp, c64(cos, sin) / (n * w));
    }
    (x, y)
}

/// Compute riemann_zeta for complex z with `|Im z| >= 2 pi RS_MIN_A^2` by the Riemann-Siegel
/// formula of Arias de Reyna (2011, Math. Comp. 80, 995-1009), as in mpmath's `rszeta`, or by
/// the Dirichlet series where that needs fewer terms.
///
/// For `s = sigma + i t`, `a = sqrt(t / (2 pi))`, `N = floor(a)` and `p = 1 - 2 (a - N)`,
/// `zeta(s) = R(s) + chi(s) conj(R(1 - conj(s)))`, with
/// `chi(s) = exp(-2i theta(t + i(1/2 - sigma)))` and
/// `R(s) = sum_{n <= N} n^-s + (-1)^(N-1) U a^-sigma sum_k C_k(sigma, p) / a^k`. Here
/// `U = exp(-i phi)`, and `phi = t ln(a) - t/2 - pi/8` is theta(t) without its O(1/t) terms.
fn zeta_riemann_siegel(z: C64) -> C64 {
    if z.im < 0.0 {
        return zeta_riemann_siegel(z.conj()).conj();
    }
    let (sigma, t) = (z.re, z.im);
    // a in double-double, as the error eps a of p = 1 - 2 (a - N) is not relative
    let a_dd = (FRAC_2_PI_DD.mul_pwr2(0.25) * t).sqrt();
    let a = a_dd.hi;
    let n_max = a as u32;
    let p = (1.0 - (a_dd - f64::from(n_max)).mul_pwr2(2.0)).hi;

    // The Dirichlet series is cheaper than the two sums below and the correction where its tail
    // after m terms, below m^(1 - sigma) / (sigma - 1), is below eps for some m <= N + the cost
    // of the correction.
    if sigma > 1.0 {
        let m = ((sigma - 1.0) * f64::EPSILON)
            .powf(-1.0 / (sigma - 1.0))
            .ceil();
        if m <= f64::from(n_max + RS_CORRECTION_TERMS) {
            return dirichlet_sums(z, m as u32).0;
        }
    }

    let (x, y) = dirichlet_sums(z, n_max);
    // phi in double-double, as its f64 rounding error grows like eps t ln(t)
    let phi = (Dd::new(t).log().mul_pwr2(0.5) - LN_SQRT_2PI_DD - 0.5) * t - FRAC_PI_8;
    let (sin, cos) = (-phi).sin_cos();
    let u = c64(cos, sin);

    // The derivatives F^(m)(p) = sum_k coeffs[k] p^(2k - m) for m <= 3 (RS_TERMS - 1)
    let q = p * p;
    let mut coeffs = RS_F_COEFFS;
    let mut fp = [C64::ZERO; 3 * RS_TERMS - 2];
    for (m, f) in fp.iter_mut().enumerate() {
        let even = coeffs[m.div_ceil(2)..]
            .iter()
            .rev()
            .fold(C64::ZERO, |acc, &c| acc * q + c);
        *f = if m % 2 == 0 { even } else { p * even };
        for (k, c) in coeffs.iter_mut().enumerate() {
            *c *= (2 * k) as f64 - m as f64;
        }
    }

    let cx = rs_correction(sigma, a, &fp);
    let cy = rs_correction(1.0 - sigma, a, &fp).conj();

    // chi(s) = u^2 exp(-2i delta), with delta = theta(t + i(1/2 - sigma)) - phi from the Stirling
    // series of theta (DLMF 5.11.1), whose next term 31 / (80640 t^5) is below eps / 50.
    let beta = sigma - 0.5;
    let w = c64(t, -beta);
    let delta = 0.5 * (w * clog1p(c64(0.0, -beta / t)) - C64::I * beta * (2.0 * a.ln() - 1.0))
        + (1.0 / 48.0 + 7.0 / (5760.0 * w * w)) / w;
    let e = (-2.0 * C64::I * delta).exp();
    let sign = neg_one_pow(n_max - 1);
    x + sign * u * a.powf(-sigma) * cx + u * e * (u * y + sign * a.powf(sigma - 1.0) * cy)
}

/// Smallest `|1 - 2^(1 - z)|` for the Borwein sum: dividing by it at most doubles the
/// absolute error of the eta sum.
const ETA_FACTOR_MIN: f64 = 0.5;

/// Compute riemann zeta for complex z and real part >= 0
fn zeta_right_halfplane(z: C64) -> C64 {
    if z == C64::ONE {
        return C64::NAN;
    }
    // Unlike in xsf, Euler-Maclaurin is used for |Im z| > 30 instead of 50: the error of the
    // 50-term Borwein sum grows like exp(pi |Im z| / 2) / (3 + sqrt(8))^50 (2e-7 at 1/2 + 48i).
    if z.re < 50.0 && z.im.abs() > 30.0 {
        if z.re >= 0.0 && z.re < 2.5 && z.im.abs() > MAX_IM {
            /* Euler-MacLaurin summation starts to take an unreasonable amount of time in this
             * region, so just give up and return NaN instead. */
            return C64::NAN;
        }
        // Unlike in xsf, the Riemann-Siegel formula (or a shorter Dirichlet sum) is used up to
        // |Im z| = MAX_IM: Euler-Maclaurin sums up to |Im z| / 4 terms, with phase errors
        // eps |Im z| ln(n) n^-Re z. Above MAX_IM, Euler-Maclaurin is kept for Re z >= 2.5, where
        // its sum stops early.
        if (TAU * RS_MIN_A * RS_MIN_A..=MAX_IM).contains(&z.im.abs()) {
            return zeta_riemann_siegel(z);
        }
        return zeta_euler_maclaurin(z);
    }
    // Unlike in xsf, Euler-Maclaurin is also used near the zeros 1 + 2 pi i k / ln(2) of
    // eta_factor: the Borwein sum gives eta(z) = eta_factor zeta(z), and the division
    // amplifies its rounding error by 1 / |eta_factor|.
    let eta_factor = 1.0 - rpow(2.0, 1.0 - z);
    if eta_factor.norm() < ETA_FACTOR_MIN {
        return zeta_euler_maclaurin(z);
    }
    zeta_borwein(z, eta_factor)
}

/// exp(pi*z) for complex z.
fn exppi(z: C64) -> C64 {
    let x = z.re;
    let y = z.im;
    let factor1 = c64(crate::cospi(y), crate::sinpi(y));
    let factor2 = (PI * x).exp();
    factor1 * factor2
}

/// log(sinpi(z)) using sin(z) = (exp(i*pi*z) - exp(-i*pi*z)) / 2i
///
/// No attempt is made to choose any particular branch of the logarithm.
/// This is an internal function and the intent is that this that the
/// result of log(sinpi(z)) will be added to other terms, and the sum
/// will then be exponentiated, making the choice of a specific branch
/// unnecessary.
fn logsinpi(z: C64) -> C64 {
    let result = z.sinpi().ln();
    // If it doesn't overflow, just do the regular calculation.
    //
    // NOTE: as in xsf, the condition is `!isfinite(result.imag())`, which only holds for a NaN
    // angle, so in effect the factored form below is always used.
    if result.re.is_finite() && !result.im.is_finite() {
        return result;
    }
    /* Otherwise factor before taking log. This is where we may end up
     * taking a branch other than the principal branch. */
    /* Calculating log((exp(i*pi*z) - exp(-i*pi*z)) / 2i). Factor out term
     * with larger magnitude before taking log. */
    if z.im > 0.0 {
        /* if z.imag() > 0 then, exp(-i*pi*z) has greatest magnitude. Factor it
         * out to get:
         * log(exp(-i*pi*z)*((exp(2*i*pi*z) - 1.0)/(2i)) =
         * log(exp(-i*pi*z)) + log((exp(2*i*pi*z) - 1.0)/(2i)) =
         * -i*pi*z + log((exp(2*i*pi*z) - 1.0)/(2i)) */
        return -C64::I * PI * z + ((exppi(2.0 * z * C64::I) - 1.0) / (2.0 * C64::I)).ln();
    }
    /* if z.imag() < 0 then, exp(i*pi*z) has greatest magnitude. Factor similarly
     * to above */
    C64::I * PI * z + ((1.0 - exppi(-2.0 * z * C64::I)) / (2.0 * C64::I)).ln()
}

/// Logarithm of the leading factor in the zeta reflection formula (DLMF 25.4.2):
/// zeta(z) = 2^z * pi^(z-1) * sin(pi*z/2) * gamma(1 - z) * zeta(1 - z)
/// TODO: Complexify the cephes zeta_reflection implementation, which uses
/// the lanczos approximation for the gamma function.
fn zeta_reflection_log_factor(z: C64) -> C64 {
    let t1 = z * LN_2;
    let t2 = (z - 1.0) * LN_PI;
    let t3 = logsinpi(z / 2.0);
    let t4 = cloggamma(1.0 - z);
    // Unlike in xsf, return the logarithm so the reflected zeta value can be
    // included before exponentiation.
    t1 + t2 + t3 + t4
}

fn zeta_reflection(z: C64) -> C64 {
    let mut factor = 2.0 * rpow(TAU, z - 1.0) * (z / 2.0).sinpi() * cgamma(1.0 - z);
    let result = zeta_right_halfplane(1.0 - z);
    if !factor.is_finite() {
        let log_factor = zeta_reflection_log_factor(z);
        factor = log_factor.exp();
        if !factor.is_finite() && result != C64::ONE {
            let product = factor * result;
            if product.is_nan() {
                // Unlike in xsf, combine both factors before exponentiating if the
                // overflowing product produces a NaN.
                return (log_factor + result.ln()).exp();
            }
            return product;
        }
    }
    // An exact unit result leaves the factor as is, preserving infinities in its components.
    if result == C64::ONE {
        factor
    } else {
        factor * result
    }
}

/// Taylor coefficients of zeta(z) - 1/(z - 1) at z = 0 (DLMF 25.2.4), highest degree first,
/// without the constant term 1/2.
///
/// Generated with mpmath as `[zeta(0, 1, k) / factorial(k) + 1 for k in range(13, 0, -1)]`.
const ZETA_REGULAR_TAYLOR_COEFFS: [f64; 13] = [
    6.539_904_805_891_015e-14,
    -4.092_582_630_415_831e-14,
    -9.406_895_665_666_177e-12,
    1.497_007_594_190_113_8e-10,
    -9.110_164_892_314_166e-10,
    -5.764_675_979_949_394e-9,
    1.686_158_263_892_200_7e-7,
    -1.301_146_013_959_624_2e-6,
    -1.940_896_320_456_038e-6,
    1.207_004_994_288_350_5e-4,
    -7.851_944_770_424_08e-4,
    -3.178_227_954_292_425_5e-3,
    0.081_061_466_795_327_26,
];

fn criemann_zeta(z: C64) -> C64 {
    if z.is_real() {
        return c64(crate::xsf::cephes::riemann_zeta(z.re), 0.0);
    }
    if z.re >= 0.5 {
        return zeta_right_halfplane(z);
    }
    if z.norm() < 0.5 {
        // The reflection evaluates zeta(1 - z) next to its pole. Sum the Taylor series instead,
        // with 1/(z - 1) = -1 + z/(z - 1) so that the constant -1/2 is exact.
        return -0.5 + z * (1.0 / (z - 1.0) + z.evalpoly(&ZETA_REGULAR_TAYLOR_COEFFS));
    }
    zeta_reflection(z)
}

pub trait ZetaArg: crate::sealed::Sealed {
    fn riemann_zeta(self) -> Self;
    fn zeta(self, q: f64) -> Self;
}

impl ZetaArg for f64 {
    #[inline]
    fn riemann_zeta(self) -> Self {
        crate::xsf::cephes::riemann_zeta(self)
    }

    #[inline]
    fn zeta(self, q: f64) -> Self {
        if q == 1.0 {
            return self.riemann_zeta();
        }
        crate::xsf::cephes::zeta(self, q)
    }
}

impl ZetaArg for num_complex::Complex<f64> {
    #[inline]
    fn riemann_zeta(self) -> Self {
        self.via(criemann_zeta)
    }

    #[inline]
    fn zeta(self, q: f64) -> Self {
        if q == 1.0 {
            return self.riemann_zeta();
        }
        if self.im == 0.0 {
            num_complex::Complex::new(self.re.zeta(q), 0.0)
        } else {
            // Complex input for Hurwitz Zeta is not currently supported.
            C64::NAN.into()
        }
    }
}

/// Riemann zeta function $\zeta(z)$ for real or complex $z$
///
/// Corresponds to [`scipy.special.zeta(z)`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html
///
/// # Definition
///
/// $$ \zeta(z) = \sum_{k=1}^\infty k^{-z} $$
///
/// # See also
/// - [`zeta`]: Hurwitz zeta function $\zeta(z, q)$
/// - [`zetac`]: Riemann zeta function minus one, $\zeta(x) - 1$
#[must_use]
#[inline]
pub fn riemann_zeta<T: ZetaArg>(z: T) -> T {
    z.riemann_zeta()
}

/// Hurwitz zeta function $\zeta(z, q)$ for real or complex $z$
///
/// Corresponds to [`scipy.special.zeta(z, q)`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html
///
/// # Definition
///
/// $$ \zeta(z, q) = \sum_{k=0}^\infty (k + q)^{-z} $$
///
/// # See also
/// - [`riemann_zeta`]: Riemann zeta function $\zeta(z)$
/// - [`zetac`]: Riemann zeta function minus one, $\zeta(x) - 1$
#[doc(alias = "hurwitz_zeta")]
#[must_use]
#[inline]
pub fn zeta<T: ZetaArg>(z: T, q: f64) -> T {
    z.zeta(q)
}

/// Riemann zeta function minus one, $\zeta(x) - 1$
///
/// Corresponds to [`scipy.special.zetac`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zetac.html
///
/// # Definition
///
/// $$ \zeta(x) - 1 = \sum_{k=2}^\infty k^{-x} $$
///
/// # See also
/// - [`riemann_zeta`]: Riemann zeta function $\zeta(z)$
/// - [`zeta`]: Hurwitz zeta function $\zeta(z, q)$
#[doc(alias = "riemann_zetac")]
#[must_use]
#[inline]
pub fn zetac(x: f64) -> f64 {
    crate::xsf::cephes::zetac(x)
}
