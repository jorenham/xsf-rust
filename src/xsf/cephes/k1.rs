//! Translated into pure Rust from `xsf/cephes/k1.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 2000 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the error conditions are corrected in the comments.

use crate::xsf::cephes::chbevl::chbevl;
use crate::xsf::cephes::i1::i1;

/* Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
 * in the interval [0,2].
 *
 * lim(x->0){ x(K1(x) - log(x/2) I1(x)) } = 1.
 */
const K1_A: [f64; 11] = [
    -7.023_863_479_386_288e-18,
    -2.427_449_850_519_366e-15,
    -6.666_901_694_199_329e-13,
    -1.411_488_392_633_527_8e-10,
    -2.213_387_630_734_726e-8,
    -2.433_406_141_565_968_4e-6,
    -1.730_288_957_513_052e-4,
    -6.975_723_859_639_864e-3,
    -1.226_111_808_226_571_5e-1,
    -3.531_559_607_765_449e-1,
    1.525_300_227_338_947_8,
];

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
 * in the interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
 */
const K1_B: [f64; 25] = [
    -5.756_744_483_665_017e-18,
    1.794_050_873_147_559_2e-17,
    -5.689_462_558_442_859e-17,
    1.838_093_544_366_639e-16,
    -6.057_047_248_373_319e-16,
    2.038_703_165_624_334e-15,
    -7.019_837_090_418_314e-15,
    2.477_154_424_481_304_3e-14,
    -8.976_705_182_324_994e-14,
    3.348_419_666_078_429_3e-13,
    -1.289_173_960_951_028_9e-12,
    5.139_639_673_481_73e-12,
    -2.129_967_838_427_568_3e-11,
    9.218_315_187_605_006e-11,
    -4.190_354_759_341_896_5e-10,
    2.015_049_755_197_033e-9,
    -1.034_576_246_567_809_7e-8,
    5.741_084_125_450_049_5e-8,
    -3.501_960_603_087_812_6e-7,
    2.406_484_947_837_217e-6,
    -1.936_197_974_166_083e-5,
    1.952_155_184_713_516_2e-4,
    -2.857_816_859_622_779e-3,
    1.039_237_365_768_172_4e-1,
    2.720_626_190_484_442_7,
];

/// Modified Bessel function, third kind, order one
///
/// Computes the modified Bessel function of the third kind
/// of order one of the argument.
///
/// The range is partitioned into the two intervals \[0,2\] and
/// (2, infinity).  Chebyshev polynomial expansions are employed
/// in each interval.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       8.9e-16     1.3e-16
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// k1 singular         x == 0          INFINITY
/// k1 domain           x < 0           NaN
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn k1(x: f64) -> f64 {
    if x == 0.0 {
        return f64::INFINITY;
    } else if x < 0.0 {
        return f64::NAN;
    }
    let z = 0.5 * x;

    if x <= 2.0 {
        let y = x * x - 2.0;
        return z.ln() * i1(x) + chbevl(y, &K1_A) / x;
    }

    (-x).exp() * chbevl(8.0 / x - 2.0, &K1_B) / x.sqrt()
}

/// Modified Bessel function, third kind, order one,
/// exponentially scaled
///
/// Returns exponentially scaled modified Bessel function
/// of the third kind of order one of the argument:
///
/// ```text
///      k1e(x) = exp(x) * k1(x).
/// ```
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0, 30       20000       9.4e-16     1.1e-16
/// ```
///
/// See k1().
#[allow(clippy::float_cmp)]
pub(crate) fn k1e(x: f64) -> f64 {
    if x == 0.0 {
        return f64::INFINITY;
    } else if x < 0.0 {
        return f64::NAN;
    }

    if x <= 2.0 {
        let y = x * x - 2.0;
        let y = (0.5 * x).ln() * i1(x) + chbevl(y, &K1_A) / x;
        return y * x.exp();
    }

    chbevl(8.0 / x - 2.0, &K1_B) / x.sqrt()
}
