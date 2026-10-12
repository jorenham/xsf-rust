//! Translated into pure Rust from `xsf/cephes/shichi.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Unlike in xsf, `(Shi, Chi)` is returned instead of written to output references, the (unused)
//! status code is dropped, and the range for which INFINITY is returned is corrected in the
//! comments.

use crate::xsf::cephes::chbevl::chbevl;

const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53
const SCIPY_EULER: f64 = 0.577_215_664_901_532_9;

/* x exp(-x) shi(x), inverted interval 8 to 18 */
const SHICHI_S1: [f64; 22] = [
    1.838_892_301_733_994_7e-17,
    -9.554_855_322_796_556e-17,
    2.043_261_059_808_798_8e-16,
    1.098_969_490_749_053_5e-15,
    -1.313_135_343_440_926e-14,
    5.939_762_262_643_143e-14,
    -3.471_970_104_977_491_4e-14,
    -1.400_597_646_131_171_2e-12,
    9.490_446_262_242_236e-12,
    -1.615_961_811_454_354_6e-11,
    -1.778_997_844_364_303_2e-10,
    1.354_554_697_672_469_6e-9,
    -1.032_571_217_928_195e-9,
    -3.566_996_111_149_825_4e-8,
    1.448_188_773_842_673_5e-7,
    7.820_182_151_840_513e-7,
    -5.399_191_184_038_050_5e-6,
    -3.124_582_021_689_598_6e-5,
    8.901_367_419_507_276e-5,
    2.025_584_747_438_468_7e-3,
    2.960_644_408_556_332_6e-2,
    1.118_477_510_472_570_4,
];

/* x exp(-x) shi(x), inverted interval 18 to 88 */
const SHICHI_S2: [f64; 23] = [
    -1.053_115_741_548_509_4e-17,
    2.624_460_955_963_552e-17,
    8.820_901_356_253_682e-17,
    -3.384_598_118_781_030_5e-16,
    -8.306_080_263_669_358e-16,
    3.933_978_754_370_5e-15,
    1.017_655_659_697_290_5e-14,
    -4.211_281_703_076_408e-14,
    -1.608_182_045_198_024_7e-13,
    3.347_149_541_759_945e-13,
    2.726_003_521_291_530_7e-12,
    1.668_949_547_528_390_8e-12,
    -3.492_781_410_247_309e-11,
    -1.585_806_616_664_827e-10,
    -1.792_894_371_833_556_3e-10,
    1.762_816_291_442_645_3e-9,
    1.690_502_288_794_213e-8,
    1.253_917_712_284_870_4e-7,
    1.162_299_470_686_773_3e-6,
    1.610_382_601_173_763e-5,
    3.498_103_756_010_54e-4,
    1.284_780_652_596_476_1e-2,
    1.036_657_225_887_983_4,
];

/* x exp(-x) chin(x), inverted interval 8 to 18 */
const SHICHI_C1: [f64; 23] = [
    -8.124_353_852_258_64e-18,
    2.175_864_132_903_392e-17,
    5.226_243_949_240_722e-17,
    -9.488_121_105_916_906e-16,
    5.355_463_116_474_652e-15,
    -1.210_099_701_137_329_2e-14,
    -6.008_651_785_534_474e-14,
    7.163_396_491_560_285e-13,
    -2.934_960_726_075_998_7e-12,
    -1.403_594_381_364_912_7e-12,
    8.763_022_886_090_549e-11,
    -4.400_924_762_132_823_5e-10,
    -1.879_920_756_405_693e-10,
    1.314_581_509_894_745_9e-8,
    -4.755_139_309_247_654_6e-8,
    -2.217_750_188_018_488_7e-7,
    1.946_355_313_732_725e-6,
    4.335_058_892_573_164e-6,
    -6.133_870_010_764_944e-5,
    -3.130_854_774_929_974_5e-4,
    4.971_647_898_231_161e-4,
    2.643_474_960_313_745_4e-2,
    1.114_461_508_766_992_2,
];

/* x exp(-x) chin(x), inverted interval 18 to 88 */
const SHICHI_C2: [f64; 24] = [
    8.069_134_082_551_556e-18,
    -2.080_741_681_801_481_7e-17,
    -5.981_113_296_582_723e-17,
    2.685_339_510_859_457_5e-16,
    4.523_139_416_989_047e-16,
    -3.107_349_173_352_994_6e-15,
    -4.428_232_073_325_319_5e-15,
    3.496_396_954_108_069_4e-14,
    6.634_067_317_189_116e-14,
    -3.719_024_480_931_192_4e-13,
    -1.271_354_181_323_383e-12,
    2.748_511_419_353_153_8e-12,
    2.337_818_439_854_534_4e-11,
    2.714_360_063_776_124_4e-11,
    -2.566_001_800_003_56e-10,
    -1.610_213_751_638_034_5e-9,
    -4.725_430_648_762_717e-9,
    -3.000_951_780_286_816_7e-9,
    7.793_874_743_909_149e-8,
    1.069_427_655_664_015e-6,
    1.595_031_648_023_132e-5,
    3.495_925_751_537_78e-4,
    1.284_753_875_300_652_6e-2,
    1.036_656_939_179_342_8,
];

/*
 * Evaluate 3F0(a1, a2, a3; z)
 *
 * The series is only asymptotic, so this requires z large enough.
 */
#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn hyp3f0(a1: f64, a2: f64, a3: f64, z: f64) -> f64 {
    let m = z.powf(-1.0 / 3.0);
    let maxiter = if m < 50.0 { m as i32 } else { 50 };

    let mut term = 1.0;
    let mut sum = term;
    for n in 0..maxiter {
        let n = f64::from(n);
        term *= (a1 + n) * (a2 + n) * (a3 + n) * z / (n + 1.0);
        sum += term;
        if term.abs() < 1e-13 * sum.abs() || term == 0.0 {
            break;
        }
    }

    let err = term.abs();

    if err > 1e-13 * sum.abs() {
        return f64::NAN;
    }

    sum
}

/// Hyperbolic sine and cosine integrals
///
/// Approximates the integrals
///
/// ```text
///                            x
///                            -
///                           | |   cosh t - 1
///   Chi(x) = eul + ln x +   |    -----------  dt,
///                         | |          t
///                          -
///                          0
///
///               x
///               -
///              | |  sinh t
///   Shi(x) =   |    ------  dt
///            | |       t
///             -
///             0
/// ```
///
/// where eul = 0.57721566490153286061 is Euler's constant.
/// The integrals are evaluated by power series for x < 8
/// and by Chebyshev expansions for x between 8 and 88.
/// For large x, both functions approach exp(x)/2x.
/// Arguments greater than 1000 in magnitude return INFINITY.
///
/// ACCURACY:
///
/// Test interval 0 to 88.
///
/// ```text
///                      Relative error:
/// arithmetic   function  # trials      peak         rms
///    IEEE         Shi      20000       7.7e-16     1.6e-16
///        Absolute error, except relative when |Chi| > 1:
///    IEEE         Chi      20000       9.5e-16     1.3e-16
/// ```
#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::neg_cmp_op_on_partial_ord
)]
pub(crate) fn shichi(mut x: f64) -> (f64, f64) {
    let sign: i16 = if x < 0.0 {
        x = -x;
        -1
    } else {
        0
    };

    if x == 0.0 {
        return (0.0, f64::NEG_INFINITY);
    }

    let (mut s, c) = 'done: {
        'asymp: {
            'chb: {
                if x >= 8.0 {
                    break 'chb;
                }

                if x >= 88.0 {
                    break 'asymp;
                }

                let z = x * x;

                /*     Direct power series expansion   */
                let mut a = 1.0;
                let mut s = 1.0;
                let mut c = 0.0;
                let mut k = 2.0;

                loop {
                    a *= z / k;
                    c += a / k;
                    k += 1.0;
                    a /= k;
                    s += a / k;
                    k += 1.0;
                    if !((a / s).abs() > MACHEP) {
                        break;
                    }
                }

                s *= x;
                break 'done (s, c);
            }

            // chb:
            /* Chebyshev series expansions */
            if x < 18.0 {
                let a = (576.0 / x - 52.0) / 10.0;
                let k = x.exp() / x;
                let s = k * chbevl(a, &SHICHI_S1);
                let c = k * chbevl(a, &SHICHI_C1);
                break 'done (s, c);
            }

            if x <= 88.0 {
                let a = (6336.0 / x - 212.0) / 70.0;
                let k = x.exp() / x;
                let s = k * chbevl(a, &SHICHI_S2);
                let c = k * chbevl(a, &SHICHI_C2);
                break 'done (s, c);
            }
        }

        // asymp:
        let (mut si, ci);
        if x > 1000.0 {
            si = f64::INFINITY;
            ci = f64::INFINITY;
        } else {
            /* Asymptotic expansions
             * http://functions.wolfram.com/GammaBetaErf/CoshIntegral/06/02/
             * http://functions.wolfram.com/GammaBetaErf/SinhIntegral/06/02/0001/
             */
            let a = hyp3f0(0.5, 1.0, 1.0, 4.0 / (x * x));
            let b = hyp3f0(1.0, 1.0, 1.5, 4.0 / (x * x));
            si = x.cosh() / x * a + x.sinh() / (x * x) * b;
            ci = x.sinh() / x * a + x.cosh() / (x * x) * b;
        }
        if sign != 0 {
            si = -si;
        }
        return (si, ci);
    };

    // done:
    if sign != 0 {
        s = -s;
    }

    let si = s;

    let ci = SCIPY_EULER + x.ln() + c;
    (si, ci)
}
