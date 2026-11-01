//! Translated into pure Rust from `xsf/cephes/zetac.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.1:  January, 1989
//! Copyright 1984, 1987, 1989 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Unlike in xsf (scipy/xsf#346), `zetac_positive` only returns 0 for `x >= 1075` (where `2**-x`
//! rounds to 0), instead of for `x >= 127` (a leftover from single precision).

use std::f64::consts::{E, PI};

use crate::xsf::cephes::lanczos::{LANCZOS_G, lanczos_sum_expg_scaled};
use crate::xsf::cephes::polevl::{p1evl, polevl};
use crate::xsf::cephes::zeta;

// cephes/const.h
const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53
const SQRT2OPI: f64 = 0.797_884_560_802_865_4; // sqrt(2/pi)

/// Riemann zeta(x) - 1
/// for integer arguments between 0 and 30.
const AZETAC: [f64; 31] = [
    -1.5,
    0.0, /* Not used; zetac(1.0) is infinity. */
    6.449_340_668_482_264e-1,
    2.020_569_031_595_943e-1,
    8.232_323_371_113_819e-2,
    3.692_775_514_336_993e-2,
    1.734_306_198_444_914e-2,
    8.349_277_381_922_827e-3,
    4.077_356_197_944_34e-3,
    2.008_392_826_082_214_3e-3,
    9.945_751_278_180_853e-4,
    4.941_886_041_194_645e-4,
    2.460_865_533_080_483e-4,
    1.227_133_475_784_891_5e-4,
    6.124_813_505_870_483e-5,
    3.058_823_630_702_049e-5,
    1.528_225_940_865_187e-5,
    7.637_197_637_899_763e-6,
    3.817_293_264_999_84e-6,
    1.908_212_716_553_939e-6,
    9.539_620_338_727_962e-7,
    4.769_329_867_878_064e-7,
    2.384_505_027_277_33e-7,
    1.192_199_259_653_110_6e-7,
    5.960_818_905_125_948e-8,
    2.980_350_351_465_228e-8,
    1.490_155_482_836_504_3e-8,
    7.450_711_789_835_43e-9,
    3.725_334_024_788_457e-9,
    1.862_659_723_513_049e-9,
    9.313_274_324_196_682e-10,
];

/// 2**x (1 - 1/x) (zeta(x) - 1) = P(1/x)/Q(1/x), 1 <= x <= 10
const ZETAC_P: [f64; 9] = [
    5.857_465_145_697_253e11,
    2.575_341_277_561_025_7e11,
    4.877_811_595_679_482e10,
    5.153_995_380_238_857e9,
    3.416_460_735_147_541e8,
    1.608_370_068_806_565e7,
    5.927_854_673_421_095e5,
    1.511_291_699_649_388_3e4,
    2.018_224_444_859_979_5e2,
];

const ZETAC_Q: [f64; 8] = [
    /*  1.00000000000000000000E0, */
    3.904_976_763_733_711_5e11,
    5.228_582_353_682_722e10,
    5.644_515_172_712_806e9,
    3.390_067_460_153_504e8,
    1.794_103_715_001_264_6e7,
    5.666_668_251_313_848e5,
    1.603_829_768_109_441_3e4,
    1.964_362_372_233_873e2,
];

/// log(zeta(x) - 1 - 2**-x), 10 <= x <= 50
const ZETAC_A: [f64; 11] = [
    8.707_285_674_845_902e6,
    1.765_068_656_703_464_7e8,
    2.608_895_067_074_832_5e10,
    5.298_063_740_098_948e11,
    2.268_881_561_192_382_4e13,
    3.318_844_029_327_050_6e14,
    5.137_789_979_758_682e15,
    -1.981_236_881_339_071_8e15,
    -9.927_638_100_399_835e16,
    7.829_053_761_808_706e16,
    9.267_862_757_689_277e16,
];

const ZETAC_B: [f64; 10] = [
    /* 1.00000000000000000000E0, */
    -7.926_254_105_637_411e6,
    -1.605_299_699_329_202_2e8,
    -2.376_692_609_755_432e10,
    -4.803_195_843_504_552e11,
    -2.078_209_617_541_733_2e13,
    -2.960_754_045_072_722_5e14,
    -4.862_991_036_946_091e15,
    5.345_895_096_757_899e15,
    5.714_641_110_922_976e16,
    -1.799_155_976_586_765_6e16,
];

/// (1-x) (zeta(x) - 1), 0 <= x <= 1
const ZETAC_R: [f64; 6] = [
    -3.287_174_745_065_627_5e-1,
    1.551_625_287_426_239_6e1,
    -2.487_628_316_808_219_6e2,
    1.010_503_680_532_376_7e3,
    1.267_260_614_102_351_4e4,
    -1.115_780_947_705_151_8e5,
];

const ZETAC_S: [f64; 5] = [
    /* 1.00000000000000000000E0, */
    1.951_076_749_140_605_3e1,
    3.177_103_117_506_47e2,
    3.038_355_008_744_457_3e3,
    2.036_658_764_357_706e4,
    7.438_539_651_367_678e4,
];

const ZETAC_TAYLOR0: [f64; 10] = [
    -1.000_000_000_911_016_6,
    -1.000_000_005_764_676,
    -9.999_998_313_841_736e-1,
    -1.000_001_301_146_014,
    -1.000_001_940_896_320_6,
    -9.998_792_995_005_712e-1,
    -1.000_785_194_477_042_5,
    -1.003_178_227_954_292_5,
    -9.189_385_332_046_728e-1,
    -1.5,
];

const ZETAC_MAXL2: f64 = 1075.0;

/// Compute zetac for positive arguments
#[allow(
    clippy::float_cmp,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::neg_cmp_op_on_partial_ord
)]
fn zetac_positive(x: f64) -> f64 {
    let mut a;
    let mut b;
    let mut s;
    let mut w;

    if x == 1.0 {
        return f64::INFINITY;
    }

    if x >= ZETAC_MAXL2 {
        /* because first term is 2**-x */
        return 0.0;
    }

    /* Tabulated values for integer argument */
    w = x.floor();
    if w == x {
        let i = x as usize;
        if i < 31 {
            return AZETAC[i];
        }
    }

    if x < 1.0 {
        w = 1.0 - x;
        a = polevl(x, &ZETAC_R) / (w * p1evl(x, &ZETAC_S));
        return a;
    }

    if x <= 10.0 {
        b = x.exp2() * (x - 1.0);
        w = 1.0 / x;
        s = (x * polevl(w, &ZETAC_P)) / (b * p1evl(w, &ZETAC_Q));
        return s;
    }

    if x <= 50.0 {
        b = (-x).exp2();
        w = polevl(x, &ZETAC_A) / p1evl(x, &ZETAC_B);
        w = w.exp() + b;
        return w;
    }

    /* Basic sum of inverse powers */
    s = 0.0;
    a = 1.0;
    loop {
        a += 2.0;
        b = a.powf(-x);
        s += b;
        if !(b / s > MACHEP) {
            break;
        }
    }

    b = (-x).exp2();
    s = (s + b) / (1.0 - b);
    s
}

/// Compute zetac for small negative x. We can't use the reflection
/// formula because to double precision 1 - x = 1 and zetac(1) = inf.
fn zetac_smallneg(x: f64) -> f64 {
    polevl(x, &ZETAC_TAYLOR0)
}

/// Compute zetac using the reflection formula (see DLMF 25.4.2) plus
/// the Lanczos approximation for Gamma to avoid overflow.
#[allow(clippy::float_cmp)]
fn zeta_reflection(x: f64) -> f64 {
    let hx = x / 2.0;
    if hx == hx.floor() {
        /* Hit a zero of the sine factor */
        return 0.0;
    }

    /* Reduce the argument to sine */
    let x_shift = x % 4.0;
    let mut small_term = -SQRT2OPI * (0.5 * PI * x_shift).sin();
    small_term *= lanczos_sum_expg_scaled(x + 1.0) * zeta(x + 1.0, 1.0);

    /* Group large terms together to prevent overflow */
    let base = (x + LANCZOS_G + 0.5) / (2.0 * PI * E);
    let mut large_term = base.powf(x + 0.5);
    if large_term.is_finite() {
        return large_term * small_term;
    }
    /*
     * We overflowed, but we might be able to stave off overflow by
     * factoring in the small term earlier. To do this we compute
     *
     * (sqrt(large_term) * small_term) * sqrt(large_term)
     *
     * Since we only call this method for negative x bounded away from
     * zero, the small term can only be as small sine on that region;
     * i.e. about machine epsilon. This means that if the above still
     * overflows, then there was truly no avoiding it.
     */
    large_term = base.powf(0.5 * x + 0.25);
    (large_term * small_term) * large_term
}

/// Riemann zeta function, minus one
///
/// ```text
///                inf.
///                 -    -x
///   zetac(x)  =   >   k   ,   x > 1,
///                 -
///                k=2
/// ```
///
/// is related to the Riemann zeta function by
///
/// ```text
///     Riemann zeta(x) = zetac(x) + 1.
/// ```
///
/// Extension of the function definition for x < 1 is implemented.
/// Zero is returned for x > log2(INFINITY).
/// (Unlike in xsf, where this is for x >= 127 instead of x >= 1075, see scipy/xsf#346.)
///
/// ACCURACY:
///
/// Tabulated values have full machine accuracy.
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      1,50        10000       9.6e-16            1.3e-16
/// ```
pub(crate) fn zetac(x: f64) -> f64 {
    if x.is_nan() {
        x
    } else if x == f64::NEG_INFINITY {
        f64::NAN
    } else if x < 0.0 && x > -0.01 {
        zetac_smallneg(x)
    } else if x < 0.0 {
        zeta_reflection(-x) - 1.0
    } else {
        zetac_positive(x)
    }
}

/// Riemann zeta function
pub(crate) fn riemann_zeta(x: f64) -> f64 {
    if x.is_nan() {
        x
    } else if x == f64::NEG_INFINITY {
        f64::NAN
    } else if x < 0.0 && x > -0.01 {
        1.0 + zetac_smallneg(x)
    } else if x < 0.0 {
        zeta_reflection(-x)
    } else {
        1.0 + zetac_positive(x)
    }
}
