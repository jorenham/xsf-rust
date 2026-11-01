use core::f64::consts::{FRAC_PI_4, PI};

use num_complex::c64;

#[test]
fn test_sinpi_f64() {
    xsref::test("sinpi", "d-d", |x| xsf::sinpi(x[0]));
}

#[test]
fn test_sinpi_c64() {
    xsref::test("sinpi", "cd-cd", |x| xsf::sinpi(c64(x[0], x[1])));
}

#[test]
fn test_cospi_f64() {
    xsref::test("cospi", "d-d", |x| xsf::cospi(x[0]));
}

#[test]
fn test_cospi_c64() {
    xsref::test("cospi", "cd-cd", |x| xsf::cospi(c64(x[0], x[1])));
}

#[test]
fn test_sindg() {
    xsref::test("sindg", "d-d", |x| xsf::sindg(x[0]));
}

#[test]
fn test_cosdg() {
    xsref::test("cosdg", "d-d", |x| xsf::cosdg(x[0]));
}

#[test]
fn test_tandg() {
    xsref::test("tandg", "d-d", |x| xsf::tandg(x[0]));
}

#[test]
fn test_cotdg() {
    xsref::test("cotdg", "d-d", |x| xsf::cotdg(x[0]));
}

#[test]
fn test_cosm1() {
    xsref::test("cosm1", "d-d", |x| xsf::cosm1(x[0]));
}

#[test]
fn test_radian() {
    xsref::test("radian", "d_d_d-d", |x| xsf::radian(x[0], x[1], x[2]));
}

// Edge cases that the xsref tables don't cover (they skip large inputs and outputs, and most
// NaNs, and don't distinguish signed zeros).

fn assert_close(actual: f64, expected: f64) {
    let ok = if expected.is_finite() {
        (actual - expected).abs() <= 1e-14 * expected.abs()
    } else {
        actual.to_bits() == expected.to_bits() || (actual.is_nan() && expected.is_nan())
    };
    assert!(ok, "actual: {actual:e}, expected: {expected:e}");
}

fn next_up(x: f64) -> f64 {
    f64::from_bits(x.to_bits() + 1)
}

#[test]
fn test_nan() {
    for f in [
        xsf::sinpi::<f64>,
        xsf::cospi::<f64>,
        xsf::sindg,
        xsf::cosdg,
        xsf::tandg,
        xsf::cotdg,
        xsf::cosm1,
    ] {
        assert!(f(f64::NAN).is_nan());
    }
    assert!(xsf::radian(f64::NAN, 0.0, 0.0).is_nan());
    assert!(xsf::sinpi(f64::INFINITY).is_nan());
    assert!(xsf::cospi(f64::NEG_INFINITY).is_nan());
    assert!(xsf::cosm1(f64::INFINITY).is_nan());
}

#[test]
fn test_sinpi_cospi_exact() {
    assert_eq!(xsf::sinpi(0.0_f64).to_bits(), 0.0_f64.to_bits());
    assert_eq!(xsf::sinpi(-0.0_f64).to_bits(), (-0.0_f64).to_bits());
    assert_eq!(xsf::sinpi(0.5_f64), 1.0);
    assert_eq!(xsf::sinpi(-0.5_f64), -1.0);
    assert_eq!(xsf::sinpi(1.5_f64), -1.0);
    assert_eq!(xsf::sinpi(1e300_f64), 0.0);
    assert_eq!(xsf::cospi(0.0_f64), 1.0);
    assert_eq!(xsf::cospi(-0.0_f64), 1.0);
    assert_eq!(xsf::cospi(1.0_f64), -1.0);
    assert_eq!(xsf::cospi(-3.0_f64), -1.0);
    // never -0.0
    for x in [0.5, -0.5, 2.5, -2.5, 1e15 + 0.5] {
        assert_eq!(xsf::cospi(x).to_bits(), 0.0_f64.to_bits(), "cospi({x})");
    }
}

#[test]
fn test_sindg_cosdg_exact() {
    for (x, s, c) in [
        (0.0, 0.0, 1.0),
        (90.0, 1.0, 0.0),
        (180.0, 0.0, -1.0),
        (270.0, -1.0, 0.0),
        (360.0, 0.0, 1.0),
        (-90.0, -1.0, 0.0),
        (-180.0, 0.0, -1.0),
        (450.0, 1.0, 0.0),
        (-720.0, 0.0, 1.0),
    ] {
        assert_eq!(xsf::sindg(x), s, "sindg({x})");
        assert_eq!(xsf::cosdg(x), c, "cosdg({x})");
    }
    assert_close(xsf::sindg(30.0), 0.5);
    assert_close(xsf::cosdg(60.0), 0.5);
    assert_close(xsf::sindg(-45.0), -FRAC_PI_4.sin());
    assert_close(xsf::cosdg(135.0), -FRAC_PI_4.cos());
}

#[test]
fn test_tandg_cotdg_exact() {
    for (x, t) in [
        (0.0, 0.0),
        (45.0, 1.0),
        (135.0, -1.0),
        (-45.0, -1.0),
        (180.0, 0.0),
        (225.0, 1.0),
    ] {
        assert_eq!(xsf::tandg(x), t, "tandg({x})");
    }
    for (x, c) in [
        (45.0, 1.0),
        (90.0, 0.0),
        (135.0, -1.0),
        (-45.0, -1.0),
        (270.0, 0.0),
    ] {
        assert_eq!(xsf::cotdg(x), c, "cotdg({x})");
    }
    for x in [90.0, 270.0, -90.0] {
        assert!(xsf::tandg(x).is_infinite(), "tandg({x})");
    }
    for x in [0.0, 180.0, -180.0] {
        assert!(xsf::cotdg(x).is_infinite(), "cotdg({x})");
    }
    assert_close(xsf::tandg(30.0), (PI / 6.0).tan());
    assert_close(xsf::cotdg(60.0), 1.0 / (PI / 3.0).tan());
}

#[test]
fn test_tandg_cotdg_near_poles() {
    // mpmath reference values
    let x = [89.999_999_999_9, 270.000_000_001];
    let t = [572_947_570_131.580_4, -57_296_385_445.315_09];
    xsref::np_assert_allclose!(&x.map(xsf::tandg), &t, rtol = 1e-15);

    let x = [1e-10, 1e-15, -1e-300, 179.999_999_999_9];
    let c = [
        572_957_795_130.823_2,
        5.729_577_951_308_232e16,
        -5.729_577_951_308_232e301,
        -573_029_000_997.147_8,
    ];
    xsref::np_assert_allclose!(&x.map(xsf::cotdg), &c, rtol = 1e-15);
}

#[test]
fn test_degrees_loss_threshold() {
    let x = 1e14;
    assert!(xsf::sindg(x).abs() <= 1.0 && xsf::sindg(x) != 0.0);
    assert!(xsf::cosdg(x).abs() <= 1.0 && xsf::cosdg(x) != 0.0);
    assert!(xsf::tandg(x).is_finite() && xsf::tandg(x) != 0.0);
    assert!(xsf::cotdg(x).is_finite() && xsf::cotdg(x) != 0.0);
    // no result: total loss of precision
    for x in [next_up(1e14), -next_up(1e14), 1e300] {
        assert_eq!(xsf::sindg(x), 0.0);
        assert_eq!(xsf::cosdg(x), 0.0);
        assert_eq!(xsf::tandg(x), 0.0);
        assert_eq!(xsf::cotdg(x), 0.0);
    }
}

#[test]
fn test_cosm1_edges() {
    assert_eq!(xsf::cosm1(0.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(xsf::cosm1(-0.0).to_bits(), 0.0_f64.to_bits());
    // the polynomial is used for |x| <= pi/4, and cos(x) - 1 beyond
    for x in [
        FRAC_PI_4,
        next_up(FRAC_PI_4),
        -FRAC_PI_4,
        -next_up(FRAC_PI_4),
    ] {
        assert_close(xsf::cosm1(x), x.cos() - 1.0);
    }
    assert_close(xsf::cosm1(1e-5), -5e-11 + 1e-20 / 24.0);
}

#[test]
fn test_radian_edges() {
    assert_eq!(xsf::radian(0.0, 0.0, 0.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        xsf::radian(-0.0, -0.0, -0.0).to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(xsf::radian(1.0, -60.0, 0.0), 0.0);
    assert_eq!(xsf::radian(-1.0, 59.0, 60.0), 0.0);
    assert_close(xsf::radian(180.0, 0.0, 0.0), PI);
    assert_close(xsf::radian(-90.0, 0.0, 0.0), -PI / 2.0);
    assert_close(xsf::radian(0.0, 0.0, 3600.0), PI / 180.0);
    assert_close(xsf::radian(0.0, -30.0, 0.0), -PI / 360.0);
}

#[test]
fn test_sinpi_cospi_complex_small_imag() {
    // |pi y| < 700, both signs of x and y, just below the threshold
    for (x, y) in [
        (0.25, 222.8),
        (-0.25, -222.8),
        (0.75, -100.0),
        (-1.25, 3.0),
        (0.1, 0.0),
    ] {
        let (s, c) = (xsf::sinpi(c64(x, y)), xsf::cospi(c64(x, y)));
        let (sx, cx) = ((PI * x).sin(), (PI * x).cos());
        let (chy, shy) = ((PI * y).cosh(), (PI * y).sinh());
        assert_close(s.re, sx * chy);
        assert_close(s.im, cx * shy);
        assert_close(c.re, cx * chy);
        assert_close(c.im, -sx * shy);
    }
}

/// `a * cosh(t)` (`sgn = 1`) or `a * sinh(t)` (`sgn = sgn(t)`) for `|t| >= 700`, in the log
/// domain, so that it's independent of the `exp(|t|/2)` scaling in `csinpi` and `ccospi`
fn mul_cosh_sinh_large(a: f64, t: f64, sgn: f64) -> f64 {
    if a == 0.0 {
        0.0_f64.copysign(a * sgn)
    } else {
        a.signum() * sgn * (a.abs().ln() + t.abs() - core::f64::consts::LN_2).exp()
    }
}

fn assert_same_or_close(actual: f64, expected: f64, what: &str) {
    let ok = if expected.is_finite() && expected != 0.0 {
        (actual - expected).abs() <= 1e-12 * expected.abs()
    } else {
        actual.to_bits() == expected.to_bits() || (actual.is_nan() && expected.is_nan())
    };
    assert!(ok, "{what}: actual: {actual:e}, expected: {expected:e}");
}

#[test]
fn test_sinpi_cospi_complex_large_imag() {
    // 700 <= |pi y|: below and beyond the overflow of exp(|pi y| / 2) at |y| ~ 451.9
    let xs = [
        0.0, -0.0, 0.25, -0.25, 0.5, -0.5, 0.75, 1.0, -1.0, 1.75, 2.5, 1e-10, -3e-200,
    ];
    let ys = [
        222.9,
        225.0,
        300.0,
        451.0,
        452.0,
        1000.0,
        1e300,
        f64::INFINITY,
    ];
    for x in xs {
        for y in ys.into_iter().flat_map(|y| [y, -y]) {
            let z = c64(x, y);
            let (sx, cx) = (xsf::sinpi(x), xsf::cospi(x));
            let (t, sgn) = (PI * y, y.signum());
            let (s, c) = (xsf::sinpi(z), xsf::cospi(z));
            assert_same_or_close(
                s.re,
                mul_cosh_sinh_large(sx, t, 1.0),
                &format!("sinpi({z}).re"),
            );
            assert_same_or_close(
                s.im,
                mul_cosh_sinh_large(cx, t, sgn),
                &format!("sinpi({z}).im"),
            );
            assert_same_or_close(
                c.re,
                mul_cosh_sinh_large(cx, t, 1.0),
                &format!("cospi({z}).re"),
            );
            assert_same_or_close(
                c.im,
                mul_cosh_sinh_large(-sx, t, sgn),
                &format!("cospi({z}).im"),
            );
        }
    }

    // the cases that xsf v0.2.2 got wrong
    let z = xsf::cospi(c64(0.0, 1000.0));
    assert_eq!(
        (z.re, z.im.to_bits()),
        (f64::INFINITY, (-0.0_f64).to_bits())
    );
    let z = xsf::sinpi(c64(0.0, -300.0));
    assert_eq!((z.re, z.im), (0.0, f64::NEG_INFINITY));
    let z = xsf::cospi(c64(0.5, 1000.0));
    assert_eq!(
        (z.re.to_bits(), z.im),
        (0.0_f64.to_bits(), f64::NEG_INFINITY)
    );
    let z = xsf::cospi(c64(0.25, 225.0));
    assert!(z.re > 0.0 && z.im < 0.0);

    // subnormal sin(pi x)
    let x = f64::from_bits(1);
    let (sx, z) = (xsf::sinpi(x), c64(x, 222.9));
    let expected = mul_cosh_sinh_large(sx, PI * z.im, 1.0);
    assert_same_or_close(xsf::sinpi(z).re, expected, "sinpi(5e-324+222.9i).re");
    assert_same_or_close(xsf::cospi(z.conj()).im, expected, "cospi(5e-324-222.9i).im");

    // NaN or infinite real part, or NaN imaginary part
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for y in [300.0, -1000.0, f64::INFINITY] {
            let (s, c) = (xsf::sinpi(c64(x, y)), xsf::cospi(c64(x, y)));
            assert!(s.re.is_nan() && s.im.is_nan(), "sinpi({x}+{y}i) = {s}");
            assert!(c.re.is_nan() && c.im.is_nan(), "cospi({x}+{y}i) = {c}");
        }
    }
    for x in [0.0, -0.5, 0.25, f64::NAN] {
        let (s, c) = (xsf::sinpi(c64(x, f64::NAN)), xsf::cospi(c64(x, f64::NAN)));
        assert!(s.re.is_nan() && s.im.is_nan(), "sinpi({x}+NaNi) = {s}");
        assert!(c.re.is_nan() && c.im.is_nan(), "cospi({x}+NaNi) = {c}");
    }
}
