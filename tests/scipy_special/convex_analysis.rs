use xsref::np_assert_allclose;

fn get_test_xy() -> [(f64, f64); 25] {
    let z1 = [-1.0, -0.5, 0.0, 0.5, 1.0];
    let z2 = z1.map(|x| z1.map(|y| (x, y)));
    z2.as_flattened().try_into().unwrap()
}

fn map2<const N: usize>(func: impl Fn(f64, f64) -> f64, xy: [(f64, f64); N]) -> [f64; N] {
    xy.map(|(x, y)| func(x, y))
}

/// translated from `scipy.special.tests.test_basic.test_entr`
#[test]
fn test_entr() {
    fn xfunc(x: f64) -> f64 {
        if x < 0.0 {
            f64::NEG_INFINITY
        } else {
            -xsf::xlogy(x, x)
        }
    }

    let z = [f64::NEG_INFINITY, -1.0, -0.5, 0.0, 0.5, 1.0, f64::INFINITY];
    let w = z.map(xfunc);
    np_assert_allclose!(z.map(xsf::entr), w, rtol = 1e-13, atol = 1e-13);
}

/// translated from `scipy.special.tests.test_basic.test_kl_div`
#[test]
fn test_kl_div() {
    fn xfunc(x: f64, y: f64) -> f64 {
        if x < 0.0 || y < 0.0 || (y == 0.0 && x != 0.0) {
            // extension of natural domain to preserve convexity
            f64::INFINITY
        } else if (x.is_infinite() && x > 0.0) || (y.is_infinite() && y > 0.0) {
            // limits within the natural domain
            f64::INFINITY
        } else if x == 0.0 {
            y
        } else {
            xsf::xlogy(x, x / y) - x + y
        }
    }

    let z = get_test_xy();
    let w = map2(xfunc, z);
    np_assert_allclose!(map2(xsf::kl_div, z), w, rtol = 1e-13, atol = 1e-13);
}

#[test]
fn test_kl_div_ratio_overflow() {
    // x / y underflows or overflows
    let inputs = [
        (1.126_513_958_653_338_2e-184, 1.521_698_452_194_790_4e258),
        (1e300, 1e-10),
    ];
    // mpmath, x * (log(x) - log(y)) - x + y
    let expected = [1.521_698_452_194_790_4e258, 7.128_013_788_281_542e302];
    np_assert_allclose!(map2(xsf::kl_div, inputs), expected, rtol = 1e-15);
}

/// translated from `scipy.special.tests.test_basic.test_rel_entr`
#[test]
fn test_rel_entr() {
    fn xfunc(x: f64, y: f64) -> f64 {
        if x > 0.0 && y > 0.0 {
            xsf::xlogy(x, x / y)
        } else if x == 0.0 && y >= 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    }

    let z = get_test_xy();
    let w = map2(xfunc, z);
    np_assert_allclose!(map2(xsf::rel_entr, z), w, rtol = 1e-13, atol = 1e-13);
}

/// translated from `scipy.special.tests.test_basic.test_rel_entr_gh_20710_near_zero`
#[test]
fn test_rel_entr_gh_20710_near_zero() {
    // Check accuracy of inputs which are very close
    let inputs = [
        // x, y
        (0.945_665_771_343_000_1, 0.945_665_771_343_009_4),
        (0.480_660_985_647_915_15, 0.480_660_985_647_947_74),
        (0.786_048_657_854_401, 0.786_048_657_854_236_7),
    ];
    // Known values produced using `x * mpmath.log(x / y)` with dps=30
    let expected = [
        -9.325_873_406_851_269e-15,
        -3.258_504_577_274_724e-14,
        1.643_130_076_445_403_3e-13,
    ];
    np_assert_allclose!(map2(xsf::rel_entr, inputs), expected, rtol = 1e-13);
}

/// translated from `scipy.special.tests.test_basic.test_rel_entr_gh_20710_overflow`
#[test]
fn test_rel_entr_gh_20710_overflow() {
    let inputs = [
        // x, y
        // Overflow
        (4.0, 2.22e-308),
        // Underflow
        (1e-200, 1e+200),
        // Subnormal
        (2.22e-308, 1e15),
    ];
    // Known values produced using `x * mpmath.log(x / y)` with dps=30
    let expected = [
        2_839.139_983_229_607,
        -9.210_340_371_976_183e-198,
        -1.649_321_200_807_447_5e-305,
    ];
    np_assert_allclose!(map2(xsf::rel_entr, inputs), expected, rtol = 1e-13);
}

/// translated from `scipy.special.tests.test_basic.test_huber`
#[test]
#[allow(clippy::float_cmp)]
fn test_huber() {
    assert_eq!(xsf::huber(-1.0, 1.5), f64::INFINITY);
    np_assert_allclose!([xsf::huber(2.0, 1.5)], [0.5 * 1.5_f64.powi(2)]);
    np_assert_allclose!([xsf::huber(2.0, 2.5)], [2.0 * (2.5 - 0.5 * 2.0)]);
}

#[test]
#[allow(clippy::float_cmp)]
fn test_pseudo_huber() {
    assert_eq!(xsf::pseudo_huber(-1.0, 1.5), f64::INFINITY);
    np_assert_allclose!([xsf::pseudo_huber(2.0, 1.5)], [1.0]);
    np_assert_allclose!([xsf::pseudo_huber(2.0, 2.5)], [41.0_f64.sqrt() - 4.0]);
    // delta^2 underflows to 0, but r / delta doesn't
    assert_eq!(xsf::pseudo_huber(1e-170, 1e-160), 0.0);
}

#[test]
fn test_pseudo_huber_extreme() {
    // delta^2 or (r / delta)^2 over- or underflows, or |r / delta| is large
    let inputs = [
        (5.632_860_139_466_413e-99, 1.158_620_610_502_460_7e211),
        (1e200, 1.0),
        (1e160, -1e150),
        (1e-170, 1e-100),
        (1.0, 1e100),
        (f64::from_bits(1), 1e300),
        (f64::INFINITY, 1.0),
    ];
    // mpmath, r^2 / (sqrt(1 + (r / delta)^2) + 1)
    let expected = [
        6.526_347_853_663_551e112,
        0.5,
        5e299,
        1e-270,
        1e100,
        4.940_656_458_412_466e-24,
        0.5,
    ];
    np_assert_allclose!(map2(xsf::pseudo_huber, inputs), expected, rtol = 1e-15);
}
