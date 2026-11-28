use core::f64::consts::{FRAC_PI_2, PI};

use num_complex::c64;
use xsref::np_assert_allclose;

xsref::tests! {
    test_sici_f64: "sici", "d-d_d", |x| xsf::sici(x[0]);
    test_sici_c64: "sici", "cd-cd_cd", |x| xsf::sici(c64(x[0], x[1]));
    test_shichi_f64: "shichi", "d-d_d", |x| xsf::shichi(x[0]);
    test_shichi_c64: "shichi", "cd-cd_cd", |x| xsf::shichi(c64(x[0], x[1]));
}

#[test]
fn test_sici_c64_negative_real_cut() {
    // the sign of `Im z` picks the side of the cut, also for `Im z = +-0` (DLMF 6.4.4)
    for x in [-1.0, f64::NEG_INFINITY] {
        for (y, im) in [(0.0, PI), (-0.0, -PI)] {
            assert_eq!(xsf::sici(c64(x, y)).1.im, im);
        }
    }
}

#[test]
fn test_shichi_c64_negative_real_cut() {
    for x in [-1.0, -10.0] {
        for (y, im) in [(0.0, PI), (-0.0, -PI)] {
            assert_eq!(xsf::shichi(c64(x, y)).1.im, im);
        }
    }
}

#[test]
fn test_shichi_c64_beyond_expi_overflow() {
    // Arb at the exact binary64 inputs.
    let (re, im) = (1.252_999_571_039_782_9e308, 1.164_840_763_997_151e308);
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            let (shi, chi) = xsf::shichi(c64(sx * 717.0, sy * 0.75));
            np_assert_allclose!(
                [shi.re, shi.im, chi.re, chi.im],
                [sx * re, sy * im, re, sx * sy * im],
                rtol = 1e-15
            );
        }
    }

    let (shi, chi) = xsf::shichi(c64(718.0, 1.25));
    np_assert_allclose!([shi.re, chi.re], [1.472_027_305_516_963_8e308; 2], rtol = 1e-15);
    assert_eq!((shi.im, chi.im), (f64::INFINITY, f64::INFINITY));

    // A part much smaller than the overflowing other part.
    let (shi, chi) = xsf::shichi(c64(721.519_149_195_927_9, -6.283_185_307_183_285_5));
    np_assert_allclose!([shi.im, chi.im], [1.360_233_275_528_933e308; 2], rtol = 1e-15);
    assert_eq!((shi.re, chi.re), (f64::INFINITY, f64::INFINITY));
    let (shi, chi) = xsf::shichi(c64(-722.975_452_114_241_2, 1.570_635_707_205_021_7));
    let re = 1.560_223_251_349_315e308;
    np_assert_allclose!([shi.re, chi.re], [-re, re], rtol = 1e-15);
    assert_eq!((shi.im, chi.im), (f64::INFINITY, f64::NEG_INFINITY));
    // A part whose two product magnitudes sum to more than `|1/w|`.
    let (shi, chi) = xsf::shichi(c64(719.0, 618.159_024_192_202_4));
    np_assert_allclose!([shi.re, chi.re], [-1.193_437_846_459_754e308; 2], rtol = 2e-15);
    assert_eq!((shi.im, chi.im), (f64::INFINITY, f64::INFINITY));

    // Near a phase zero, the bound cannot certify a formerly overflowing part.
    let (shi, chi) = xsf::shichi(c64(755.999_999_999_540_4, 9.437_277_058_566_693));
    assert_eq!((shi.im, chi.im), (f64::NEG_INFINITY, f64::NEG_INFINITY));
}

#[test]
fn test_sici_c64_huge_re() {
    // `Ci` keeps the imaginary part of the small `E1(+-iz)`
    let z = c64(2.981_344_618_416_097e61, 55.611_119_890_015_31);
    let (si, ci) = xsf::sici(z);
    let (_, ci_neg) = xsf::sici(-z);
    let actual = [si, ci, ci_neg];
    let expected = [
        c64(FRAC_PI_2, -7.448_138_578_702_911e-39),
        c64(-7.448_138_578_702_911e-39, -2.258_051_079_572_012_8e-38),
        c64(-7.448_138_578_702_911e-39, -PI),
    ];
    np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}

#[test]
fn test_sici_c64_negative_imag_axis() {
    // `Si(-iy) = -i Shi(y)`, `Ci(-iy) = Chi(y) - pi i / 2` (DLMF 6.5.3/6.5.4), also for `y >= 40`
    for y in [39.0, 40.0, 100.0] {
        let (si, ci) = xsf::sici(c64(0.0, -y));
        assert_eq!((si.re, ci.im), (0.0, -FRAC_PI_2));
    }
}

#[test]
fn test_sici_c64_beyond_exp_overflow() {
    // `exp(710)` overflows, but Si and Ci do not (Arb)
    let (si, ci) = xsf::sici(c64(0.0, 710.0));
    let desired = [0.0, 1.575_457_844_103_100_7e305, 1.575_457_844_103_100_7e305, FRAC_PI_2];
    np_assert_allclose!([si.re, si.im, ci.re, ci.im], desired, rtol = 1e-15);
    let (si, ci) = xsf::sici(c64(1.0, 710.0));
    let (a, b) = (1.324_498_826_579_209_2e305, 8.530_916_358_913_208e304);
    np_assert_allclose!([si.re, si.im, ci.re, ci.im], [a, b, b, -a], rtol = 1e-15);
}

#[test]
fn test_sici_c64_beyond_exp1_overflow() {
    // E1 overflows at these inputs, but its half can still be finite (Arb).
    let axis = 1.710_808_329_990_549_1e308;
    let (a, b) = (1.438_301_762_813_376_4e308, 9.263_624_948_034_014e307);
    for sign_y in [-1.0_f64, 1.0] {
        for x in [-0.0, 0.0] {
            let (si, ci) = xsf::sici(c64(x, sign_y * 717.0));
            assert_eq!(si.re.to_bits(), (-sign_y * 0.0).to_bits());
            assert_eq!(ci.im, sign_y * FRAC_PI_2);
            np_assert_allclose!([si.im, ci.re], [sign_y * axis, axis], rtol = 1e-15);
        }

        for sign_x in [-1.0, 1.0] {
            let (si, ci) = xsf::sici(c64(sign_x, sign_y * 717.0));
            np_assert_allclose!(
                [si.re, si.im, ci.re, ci.im],
                [sign_x * a, sign_y * b, b, -sign_x * sign_y * a],
                rtol = 1e-15
            );
        }
    }
}

#[test]
fn test_sici_c64_partial_overflow() {
    // Each output component overflows independently (Arb).
    let (si, ci) = xsf::sici(c64(0.5, 717.2));
    let a = 1.000_240_791_103_672_7e308;
    assert_eq!((si.im, ci.re), (f64::INFINITY, f64::INFINITY));
    np_assert_allclose!([si.re, ci.im], [a, -a], rtol = 1e-15);

    let (si, ci) = xsf::sici(c64(0.0, -718.0));
    assert_eq!(si, c64(0.0, f64::NEG_INFINITY));
    assert_eq!(ci, c64(f64::INFINITY, -FRAC_PI_2));
}
