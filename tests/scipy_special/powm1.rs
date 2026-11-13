// Ported from `scipy.special.tests.test_powm1`

const INF: f64 = f64::INFINITY;
const NAN: f64 = f64::NAN;

#[test]
fn test_powm1() {
    for (x, y, expected, rtol) in [
        (1.25, 0.75, 0.182_177_011_253_969_76, 1e-15),
        (2.0, 1e-7, 6.931_472_045_825_965e-08, 1e-15),
        (25.0, 5e-11, 1.609_437_912_563_614_8e-10, 1e-15),
        (0.99996, 0.75, -3.000_015_000_253_005_8e-05, 1e-15),
        (
            0.999_999_999_999_090_5,
            20.0,
            -1.818_989_403_530_14e-11,
            1e-15,
        ),
        (-1.25, 751.0, -6.017_550_852_453_444e+72, 2e-15),
    ] {
        xsref::np_assert_allclose!(&[xsf::powm1(x, y)], &[expected], rtol = rtol);
    }
}

#[test]
fn test_powm1_exact_cases() {
    for (x, y, expected) in [
        (0.0, 0.0, 0.0),
        (0.0, -1.5, INF),
        (0.0, 1.75, -1.0),
        (-1.5, 2.0, 1.25),
        (-1.5, 3.0, -4.375),
        (NAN, 0.0, 0.0),
        (1.0, NAN, 0.0),
        (1.0, INF, 0.0),
        (1.0, -INF, 0.0),
        (INF, 7.5, INF),
        (INF, -7.5, -1.0),
        (3.25, INF, INF),
        (INF, INF, INF),
        (INF, -INF, -1.0),
        (INF, 0.0, 0.0),
        (-INF, 0.0, 0.0),
        (-INF, 2.0, INF),
        (-INF, 3.0, -INF),
        (-1.0, 2.0_f64.powi(53) - 1.0, -2.0),
    ] {
        assert_eq!(xsf::powm1(x, y), expected, "powm1({x}, {y})");
    }
}

#[test]
fn test_powm1_return_nan() {
    for (x, y) in [
        (-1.25, 751.03),
        (-1.25, INF),
        (NAN, NAN),
        (-INF, -INF),
        (-INF, 2.5),
    ] {
        assert!(xsf::powm1(x, y).is_nan(), "powm1({x}, {y})");
    }
}

#[test]
fn test_powm1_overflow() {
    // SciPy returns 0 for y < 0
    for (x, y, expected) in [
        (2.0, 2000.0, INF),
        (0.5, -2000.0, INF),
        (-0.5, -2001.0, -INF),
        (0.5, -INF, INF),
    ] {
        assert_eq!(xsf::powm1(x, y), expected, "powm1({x}, {y})");
    }
}
