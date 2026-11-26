use num_complex::c64;

use crate::DENORM_MIN;

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
fn test_complex_erfc_underflow_signs() {
    for (y, negative_re) in [(0.03125, true), (DENORM_MIN, false)] {
        let z = c64(40.0, y);
        let value = xsf::erfc(z);
        assert_eq!(value, c64(0.0, 0.0));
        assert_eq!(value.re.is_sign_negative(), negative_re);
        assert!(value.im.is_sign_negative());

        let conjugate = xsf::erfc(z.conj());
        assert_eq!(conjugate, c64(0.0, 0.0));
        assert_eq!(conjugate.re.to_bits(), value.re.to_bits());
        assert_eq!(conjugate.im.to_bits(), 0);

        let reflected = xsf::erfc(-z);
        assert_eq!(reflected.re.to_bits(), 2.0_f64.to_bits());
        assert_eq!(reflected.im.to_bits(), 0);
        assert_eq!(xsf::erfc(-z.conj()).im.to_bits(), value.im.to_bits());
    }
}

#[test]
fn test_complex_gaussian_underflow_wide_phase() {
    // Signs from the leading terms at the exact binary64 coordinates.
    for (z, tail, cdf) in [
        (c64(1e8, 9.1e7), c64(0.0, -0.0), c64(-0.0, 0.0)),
        (c64(1e200, 1e190), c64(-0.0, 0.0), c64(-0.0, -0.0)),
    ] {
        let bits = |w: num_complex::Complex64| (w.re.to_bits(), w.im.to_bits());
        assert_eq!(bits(xsf::erfc(z)), bits(tail));
        assert_eq!(bits(xsf::erfc(z.conj())), bits(tail.conj()));
        assert_eq!(bits(xsf::erfc(-z)), bits(c64(2.0, -tail.im)));
        assert_eq!(bits(xsf::ndtr(-z)), bits(cdf));
        assert_eq!(bits(xsf::ndtr(-z.conj())), bits(cdf.conj()));
        assert_eq!(bits(xsf::ndtr(z)), bits(c64(1.0, -cdf.im)));
    }
}

#[test]
fn test_complex_erfc_underflow_recentered_phase() {
    // Signs from the leading terms at the exact binary64 coordinates.
    for (z, expected) in [
        (c64(-11_347_509_302.992_588, 172_748.156_114_021_35), c64(2.0, -0.0)),
        (c64(174_188_969_231_383.62, -6.154_904_348_482_96), c64(-0.0, -0.0)),
    ] {
        let bits = |w: num_complex::Complex64| (w.re.to_bits(), w.im.to_bits());
        assert_eq!(bits(xsf::erfc(z)), bits(expected));
        assert_eq!(bits(xsf::erfc(z.conj())), bits(expected.conj()));
    }
}

#[test]
fn test_complex_erfc_overflow_signs() {
    for (z, expected) in [
        (
            c64(-0.686_256_021_888_918, 777.639_300_098_700_8),
            c64(f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        (
            c64(148.186_901_700_005_93, 21_603_009_296_368.42),
            c64(f64::NEG_INFINITY, f64::INFINITY),
        ),
        (c64(-78_937.031_298_532_31, 56_581_067_547.780_93), c64(f64::INFINITY, f64::INFINITY)),
        (c64(15.953_961_398_340_361, -2_903_588.443_042_992_6), c64(f64::INFINITY, f64::INFINITY)),
    ] {
        assert_eq!(xsf::erfc(z), expected);
        assert_eq!(xsf::erfc(z.conj()), expected.conj());
        assert_eq!(xsf::erfc(-z), -expected);
        assert_eq!(xsf::erfc(-z.conj()), -expected.conj());
    }
}

#[test]
fn test_complex_erfc_overflow_exponent_range() {
    for (z, expected) in [
        (c64(-11.443_961_205_562_948, 29.224_816_573_696_327), c64(f64::INFINITY, f64::INFINITY)),
        (
            c64(208.485_367_286_785_53, -412.523_131_637_695_6),
            c64(f64::NEG_INFINITY, f64::INFINITY),
        ),
    ] {
        assert_eq!(xsf::erfc(z), expected);
        assert_eq!(xsf::erfc(z.conj()), expected.conj());
        assert_eq!(xsf::erfc(-z), -expected);
        assert_eq!(xsf::erfc(-z.conj()), -expected.conj());
    }
}

#[test]
fn test_complex_erfc_overflow_phase_signs() {
    for (z, expected) in [
        (
            c64(22_254_854.105_771_724, 108_004_952.648_255_74),
            c64(f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        (
            c64(-98_072_459.290_448_9, -179_629_882.932_427_44),
            c64(f64::NEG_INFINITY, f64::INFINITY),
        ),
        (c64(1e-30, -8.988_465_674_311_05e307), c64(f64::NEG_INFINITY, f64::INFINITY)),
    ] {
        assert_eq!(xsf::erfc(z), expected);
        assert_eq!(xsf::erfc(z.conj()), expected.conj());
        assert_eq!(xsf::erfc(-z), -expected);
        assert_eq!(xsf::erfc(-z.conj()), -expected.conj());
    }
}

#[test]
fn test_complex_erfc_near_component_overflow() {
    let actual = xsf::erfc(c64(13.411_601_942_156_58, 29.896_363_007_581_32));
    assert!(!actual.re.is_infinite());
    assert_eq!(actual.im, f64::INFINITY);
}

#[test]
fn test_complex_erfc_rescales_finite_overflow_companion() {
    let x = 63.495_043_197_677_92;
    let y = 68.897_756_302_938_25;
    // Arb at the exact binary64 inputs, including the actual two-component overflow nearby.
    for sign_re in [-1.0, 1.0] {
        for sign_im in [-1.0, 1.0] {
            let actual = xsf::erfc(c64(sign_re * x, sign_im * y));
            xsref::np_assert_allclose!(
                [actual.re],
                [-sign_re * 1.797_693_134_861_679_2e308],
                rtol = 2e-15
            );
            assert_eq!(actual.im, sign_im * f64::INFINITY);
            assert_eq!(
                xsf::erfc(c64(sign_re * x, sign_im * 68.897_756_302_938_65)),
                c64(-sign_re * f64::INFINITY, sign_im * f64::INFINITY)
            );
        }
    }
}

#[test]
fn test_complex_erfc_rescales_nan_companion() {
    let x = 13.411_601_942_156_58;
    // Arb at the exact binary64 inputs.
    for (y, re, im) in [
        (29.896_363_007_581_32, 7.349_056_666_811_188e307, f64::INFINITY),
        (29.896_363_007_581_318, 7.349_056_666_807_914e307, 1.797_693_134_862_006_4e308),
        (29.896_363_007_581_11, 7.349_056_666_618_002e307, 1.797_693_134_843_931_6e308),
    ] {
        let actual = xsf::erfc(c64(x, y));
        xsref::np_assert_allclose!([actual.re], [re], rtol = 2e-15);
        if im.is_infinite() {
            assert_eq!(actual.im, im);
        } else {
            xsref::np_assert_allclose!([actual.im], [im], rtol = 2e-15);
        }
    }

    let swapped = xsf::erfc(c64(-36.088_703_318_231_68, 44.909_029_817_434_586));
    xsref::np_assert_allclose!(
        [swapped.re, swapped.im],
        [-1.797_693_134_862_158e308, -3.409_022_451_000_588e307],
        rtol = 2e-15
    );

    let z = c64(50.149_973_482_398_91, 56.836_113_599_097_594);
    let actual = xsf::erfc(z);
    assert_eq!(actual.re, f64::NEG_INFINITY);
    xsref::np_assert_allclose!([actual.im], [-1.587_012_934_191_036_6e308], rtol = 2e-15);
}

#[test]
fn test_complex_erfc_rescales_infinite_companion() {
    // Arb at the exact binary64 inputs.
    for sign_re in [-1.0, 1.0] {
        for sign_im in [-1.0, 1.0] {
            let z = c64(sign_re * 63.495_043_197_677_92, sign_im * 68.897_469_852_466_05);
            let actual = xsf::erfc(z);
            let expected_re = -sign_re * 1.794_691_780_639_764_7e308;
            xsref::np_assert_allclose!([actual.re], [expected_re], rtol = 2e-15);
            assert_eq!(actual.im, sign_im * f64::INFINITY);
        }
    }
}

#[test]
fn test_complex_erfc_overflow_shared_with_ndtr() {
    // Arb at the exact binary64 inputs: a small phase, and both components just above MAX.
    for (z, expected) in [
        (c64(1e-30, 300.0), c64(f64::NEG_INFINITY, f64::NEG_INFINITY)),
        (c64(-250.0, -251.435_876_517_254_4), c64(f64::NEG_INFINITY, f64::NEG_INFINITY)),
    ] {
        assert_eq!(xsf::erfc(z), expected);
        assert_eq!(xsf::erfc(z.conj()), expected.conj());
        assert_eq!(xsf::erfc(-z), -expected);
        assert_eq!(xsf::erfc(-z.conj()), -expected.conj());
    }
}

#[test]
fn test_complex_erfc_phase_product_overflow() {
    let z = c64(5.0, 8.988_465_674_311_05e307);
    let expected = c64(f64::INFINITY, f64::NEG_INFINITY);
    assert_eq!(xsf::erfc(z), expected);
    assert_eq!(xsf::erfc(z.conj()), expected.conj());
    assert_eq!(xsf::erfc(-z), -expected);
    assert_eq!(xsf::erfc(-z.conj()), -expected.conj());
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

#[test]
fn test_complex_erf_phase_product_overflow() {
    for (z, expected) in [
        (c64(-2.5, -8.988_465_674_311_05e307), c64(f64::INFINITY, f64::NEG_INFINITY)),
        (c64(10.0, -8.988_465_674_311_05e307), c64(f64::NEG_INFINITY, f64::NEG_INFINITY)),
    ] {
        assert_eq!(xsf::erf(z), expected);
        assert_eq!(xsf::erf(z.conj()), expected.conj());
        assert_eq!(xsf::erfi(c64(z.im, -z.re)), c64(expected.im, -expected.re));
    }
}

#[test]
fn test_complex_erf_diagonal_phase() {
    let boundary = f64::from_bits(0x5fe6_a09e_667f_3bcd);
    let boundary_tail = 4.162_221_635_295_825e-155_f64;
    let ulp = f64::from_bits(boundary_tail.to_bits() + 1) - boundary_tail;
    // Arb at the exact inputs, including a finite square just before overflow.
    let finite_square = f64::from_bits(0x5fe6_a09e_667f_3bcb);
    for (x, tail, atol) in [
        (2.0_f64.powi(55), -1.101_181_121_889_039_2e-17, 8.0e-33),
        (36_028_797_018_963_976.0, -1.104_582_757_127_772_3e-17, 8.0e-33),
        (2.0_f64.powi(82), 4.273_248_405_561_108e-27, 4.0e-42),
        (4.777_560_510_820_643_5e25, 5.534_338_826_722_963e-27, 4.0e-42),
        (2.0_f64.powi(86), 5.137_282_552_751_807e-27, 4.0e-42),
        (finite_square, 3.462_935_624_014_001_4e-155, 4.0 * ulp),
        (8.988_465_674_311_05e307, -1.427_172_851_852_344e-309, 0.0),
        (boundary, boundary_tail, 4.0 * ulp),
    ] {
        for sx in [-1.0, 1.0] {
            for sy in [-1.0, 1.0] {
                let z = c64(sx * x, sy * x);
                let erf = xsf::erf(z);
                let erfi = xsf::erfi(z);
                assert_eq!(erf.re, sx);
                assert_eq!(erfi.im, sy);
                assert!((erf.im - sy * tail).abs() <= atol, "erf({z}) = {erf}");
                assert!((erfi.re - sx * tail).abs() <= atol, "erfi({z}) = {erfi}");
            }
        }
    }
}

#[test]
fn test_complex_erfc_overrange_diagonal() {
    let x = 8.988_465_674_311_05e307;
    let r = f64::from_bits(0x0003_05a4_7c24_0ed3);
    let i = f64::from_bits(0x0001_06b8_2891_590f);
    let tol = 2.0 * DENORM_MIN;
    for sy in [-1.0, 1.0] {
        let right = xsf::erfc(c64(x, sy * x));
        assert!(right.re.is_subnormal() && right.im.is_subnormal());
        assert!((right.re - r).abs() <= tol);
        assert!((right.im - sy * i).abs() <= tol);
        let left = xsf::erfc(c64(-x, sy * x));
        assert_eq!(left, c64(2.0, right.im));
    }

    let x = f64::from_bits(0x5fe6_a09e_667f_3bcd);
    let z = xsf::erfc(c64(x, x));
    let r = f64::from_bits(0x9fb5_3a53_34fc_2b35);
    let i = f64::from_bits(0x9fe1_dba6_454b_b024);
    assert!((z.re - r).abs() <= 4.0 * f64::EPSILON * r.abs());
    assert!((z.im - i).abs() <= 4.0 * f64::EPSILON * i.abs());
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
