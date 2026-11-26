use num_complex::c64;

use crate::DENORM_MIN;

xsref::tests! {
    test_expm1_f64: "expm1", "d-d", |x| xsf::expm1(x[0]);
    test_expm1_c64: "expm1", "cd-cd", |x| xsf::expm1(c64(x[0], x[1]));
    test_exp2_f64: "exp2", "d-d", |x| xsf::exp2(x[0]);
    test_exp10_f64: "exp10", "d-d", |x| xsf::exp10(x[0]);
}

#[test]
fn test_expm1_c64_edges() {
    // small |z|, where exp(z) - 1 would lose precision
    let w = xsf::expm1(c64(1e-10, 1e-10));
    xsref::np_assert_allclose!([w.re, w.im], [1e-10, 1e-10 + 1e-20], rtol = 1e-15);

    // each branch: re <= -40, -40 < re <= -1, and re > -1
    for re in [-41.0_f64, -40.0, -1.5, -1.0, -0.25, 2.0] {
        let w = xsf::expm1(c64(re, 2.0));
        let expected = [re.exp() * 2.0_f64.cos() - 1.0, re.exp() * 2.0_f64.sin()];
        xsref::np_assert_allclose!([w.re, w.im], expected, rtol = 1e-14);
    }
    assert_eq!(xsf::expm1(c64(-40.0, 1.0)).re, -1.0);
}

#[test]
fn test_expm1_c64_overflow() {
    let inf = f64::INFINITY;

    // e^x overflows, but not e^x cos(y) or e^x sin(y); references from Arb
    for (z, expected) in [
        (c64(710.0, 3.0), [-inf, 3.152_613_594_064_418e307]),
        (c64(709.9, 2.0), [-8.412_000_710_440_457e307, inf]),
        (c64(710.0, 1e-300), [inf, 223_399_476.616_171_1]),
    ] {
        let w = xsf::expm1(z);
        xsref::np_assert_allclose!([w.re, w.im], expected, rtol = 1e-14);
    }
    // near pi / 2, the real part is finite but the imaginary part is not
    let w = xsf::expm1(c64(710.0, std::f64::consts::FRAC_PI_2));
    xsref::np_assert_allclose!([w.re], [1.367_927_269_845_94e292], rtol = 1e-14);
    assert_eq!(w.im, inf);

    // zero imaginary part, with its sign
    let w = xsf::expm1(c64(710.0, 0.0));
    assert_eq!((w.re, w.im.to_bits()), (inf, 0.0_f64.to_bits()));
    let w = xsf::expm1(c64(710.0, -0.0));
    assert_eq!((w.re, w.im.to_bits()), (inf, (-0.0_f64).to_bits()));
    let w = xsf::expm1(c64(1e300, 0.0));
    assert_eq!((w.re, w.im.to_bits()), (inf, 0.0_f64.to_bits()));

    // e^(x / 2) overflows too (MPFR)
    let w = xsf::expm1(c64(1420.0, DENORM_MIN));
    assert_eq!(w.re, inf);
    xsref::np_assert_allclose!([w.im], [2.465_749_532_768_495e293], rtol = 1e-14);
    let w = xsf::expm1(c64(1e300, 1.0));
    assert_eq!((w.re, w.im), (inf, inf));
}
