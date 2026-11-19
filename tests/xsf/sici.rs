use core::f64::consts::{FRAC_PI_2, PI};

use num_complex::c64;

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
    xsref::np_assert_allclose!(&actual, &expected, rtol = 1e-15);
}
