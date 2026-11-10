use core::f64::consts::PI;

use xsf::{cosine_cdf, cosine_invcdf};

// reference values from scipy/special/tests/test_cosine_distr.py

#[test]
#[allow(clippy::approx_constant)]
fn test_cosine_cdf() {
    assert_eq!(cosine_cdf(-4.0), 0.0);
    assert_eq!(cosine_cdf(0.0), 0.5);
    assert_eq!(cosine_cdf(PI), 1.0);
    assert_eq!(cosine_cdf(4.0), 1.0);
    assert!(cosine_cdf(f64::NAN).is_nan());

    let xs = [
        3.1409, 2.25, -1.599, -1.601, -2.0, -3.0, -3.1415, -3.141_59, -PI,
    ];
    let expected = [
        0.999_999_999_991_185,
        0.981_932_817_328_790_7,
        0.086_419_598_383_825_53,
        0.086_110_582_992_713,
        0.036_970_933_596_161_1,
        7.522_387_241_801_384e-5,
        2.109_869_685_443_648e-14,
        4.956_444_476_505_336e-19,
        4.871_934_450_264_861e-50,
    ];
    let actual = xs.map(cosine_cdf);
    xsref::np_assert_allclose!(&actual, &expected, rtol = 5e-15);
}

#[test]
fn test_cosine_invcdf() {
    assert_eq!(cosine_invcdf(0.0), -PI);
    assert_eq!(cosine_invcdf(0.5), 0.0);
    assert_eq!(cosine_invcdf(1.0), PI);
    assert!(cosine_invcdf(-0.1).is_nan());
    assert!(cosine_invcdf(1.1).is_nan());
    assert!(cosine_invcdf(f64::NAN).is_nan());

    let ps = [
        1e-50,
        1e-14,
        1e-8,
        0.001_800_1,
        0.01,
        0.06,
        0.125,
        0.25,
        0.4,
        0.419,
        0.421,
        0.75,
        0.94,
        0.999_999_999_6,
    ];
    let expected = [
        -PI,
        -3.141_520_413_705_845_4,
        -3.134_368_658_912_452_4,
        -2.732_563_923_138_336,
        -2.412_765_890_086_78,
        -1.788_124_497_533_015_7,
        -1.375_252_366_986_927_4,
        -0.831_711_193_579_736,
        -0.316_795_451_239_528_9,
        -0.255_860_256_269_199_06,
        -0.249_475_707_504_456_63,
        0.831_711_193_579_736,
        1.788_124_497_533_015_3,
        3.139_122_083_991_716_7,
    ];
    let actual = ps.map(cosine_invcdf);
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-14);
}
