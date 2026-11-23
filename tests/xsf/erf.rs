use num_complex::c64;

xsref::tests! {
    // erf
    test_erf_f64: "erf", "d-d", |x| xsf::erf(x[0]);
    test_erf_c64: "erf", "cd-cd", |x| xsf::erf(c64(x[0], x[1]));

    // erfc
    test_erfc_f64: "erfc", "d-d", |x| xsf::erfc(x[0]);
    test_erfc_c64: "erfc", "cd-cd", |x| xsf::erfc(c64(x[0], x[1]));

    // erfcx
    test_erfcx_f64: "erfcx", "d-d", |x| xsf::erfcx(x[0]);
    test_erfcx_c64: "erfcx", "cd-cd", |x| xsf::erfcx(c64(x[0], x[1]));

    // erfi
    test_erfi_f64: "erfi", "d-d", |x| xsf::erfi(x[0]);
    test_erfi_c64: "erfi", "cd-cd", |x| xsf::erfi(c64(x[0], x[1]));
}

#[test]
fn test_erfc_f64_tail() {
    // mpmath and Arb, at the f64 values of x, where rounding x^2 would cost ~x^2 ulp
    let x = [10.5, 15.3, 20.1, 25.7, 26.5];
    let desired = [
        7.035_928_090_177_523e-50,
        7.976_560_125_420_7e-104,
        9.735_977_140_428_878e-178,
        3.118_899_933_007_384e-289,
        2.210_907_664_263_734e-307,
    ];
    xsref::np_assert_allclose!(x.map(xsf::erfc), desired, rtol = 2e-15);
}

#[test]
fn test_erfc_f64_subnormal() {
    // mpmath
    let x = [26.65, 26.7, 27.0, 27.2, 27.3];
    let desired = [7.580_711_477_038_6e-311, 5.253_110_413_596e-312, 5.237_05e-319, 1e-323, 0.0];
    xsref::np_assert_allclose!(x.map(xsf::erfc), desired, rtol = 1e-12, atol = 5e-324);
    assert_eq!(x.map(|x| xsf::erfc(-x)), [2.0; 5]);
}

#[test]
fn test_erfi_f64_near_overflow() {
    // mpmath
    let x = [26.65, 26.7, 26.71];
    let desired = [5.912_168_581_618_709e306, 8.499_867_261_268_985e307, 1.449_459_118_932_731e308];
    xsref::np_assert_allclose!(x.map(xsf::erfi), desired, rtol = 1e-13);
    xsref::np_assert_allclose!(x.map(|x| xsf::erfi(-x)), desired.map(|y| -y), rtol = 1e-13);
    assert_eq!(
        [26.72, -26.72, 1e200].map(xsf::erfi),
        [f64::INFINITY, -f64::INFINITY, f64::INFINITY]
    );
}

#[test]
fn test_complex_erf_axis_near_overflow() {
    // The exact f64 inputs were evaluated with mpmath and Arb.
    for (y, expected) in [
        (26.672_633_112_078_042, 1.974_735_818_006_546e307),
        (26.675_331_171_194_614, 2.280_213_994_097_188_6e307),
        (26.704_663_756_633_916, 1.090_195_355_963_869_6e308),
        (26.714_033_109_640_933, 1.797_693_134_861_981e308),
    ] {
        for sign in [-1.0, 1.0] {
            let y = sign * y;
            let expected = sign * expected;
            let w = xsf::erf(c64(-0.0, y));
            let wc = xsf::erfc(c64(-0.0, y));
            let wi = xsf::erfi(c64(y, -0.0));
            assert_eq!(w.re.to_bits(), (-0.0f64).to_bits());
            assert_eq!(wc.re, 1.0);
            assert_eq!(wi.im.to_bits(), (-0.0f64).to_bits());
            for value in [w.im, -wc.im, wi.re, xsf::erfi(y)] {
                assert!((value / expected - 1.0).abs() < 2e-15, "y = {y}: {value:e}");
            }
        }
    }

    for y in [26.714_033_109_640_937, f64::INFINITY] {
        assert_eq!(xsf::erf(c64(0.0, y)), c64(0.0, f64::INFINITY));
        assert_eq!(xsf::erfc(c64(0.0, -y)), c64(1.0, f64::INFINITY));
        assert_eq!(xsf::erfi(c64(-y, 0.0)), c64(f64::NEG_INFINITY, 0.0));
        assert_eq!(xsf::erfi(y), f64::INFINITY);
    }
}

#[test]
fn test_complex_erf_overflow_phase() {
    // The signs at the exact f64 inputs agree in mpmath and Arb.
    for (z, expected) in [
        (
            c64(9.569_339_454_174_042e-8, 2.860_739_775_716_457_3e25),
            c64(f64::INFINITY, f64::NEG_INFINITY),
        ),
        (
            c64(-15_587_454_079_821.066, 1.024_736_831_829_340_5e46),
            c64(f64::INFINITY, f64::NEG_INFINITY),
        ),
        (
            c64(165_111.046_468_009_44, -129_479_747_721.267_38),
            c64(f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
    ] {
        assert_eq!(xsf::erf(z), expected);
        assert_eq!(xsf::erfi(c64(z.im, -z.re)), c64(expected.im, -expected.re));
    }
}

// erf_zeros (from `scipy.speceial.tests.test_basic.TestErf.test_erf_zeros`)

#[test]
fn test_erf_zeros() {
    let erz = xsf::erf_zeros::<5>();
    let erzr = [
        c64(1.450_616_16, 1.880_943_00),
        c64(2.244_659_28, 2.616_575_14),
        c64(2.839_741_05, 3.175_628_10),
        c64(3.335_460_74, 3.646_174_38),
        c64(3.769_005_57, 4.060_697_23),
    ];
    xsref::np_assert_allclose!(&erz, &erzr, atol = 1.5e-4);
}

#[test]
fn test_erf_zeros_large() {
    // Zeros from mpmath.findroot of mpmath.erf, starting at the asymptotic estimates.
    let erz = xsf::erf_zeros::<1000>();
    xsref::np_assert_allclose!(
        &[erz[76], erz[99], erz[399], erz[999]],
        &[
            c64(15.481_740_840_635_416, 15.599_573_700_766_374),
            c64(17.659_970_331_721_766, 17.767_043_069_336_868),
            c64(35.411_898_828_807_104, 35.475_194_975_359_294),
            c64(56.024_353_733_827_17, 56.068_470_016_470_09),
        ],
        rtol = 1e-15
    );
}

xsref::tests! {
    // dawsn
    test_dawsn_f64: "dawsn", "d-d", |x| xsf::dawsn(x[0]);
    test_dawsn_c64: "dawsn", "cd-cd", |x| xsf::dawsn(c64(x[0], x[1]));

    // wofz
    test_wofz: "wofz", "cd-cd", |x| xsf::wofz(c64(x[0], x[1]));

    // voigt_profile
    test_voigt_profile: "voigt_profile", "d_d_d-d", |x| xsf::voigt_profile(x[0], x[1], x[2]);
}
