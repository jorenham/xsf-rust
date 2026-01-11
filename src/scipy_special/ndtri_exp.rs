//! Ported from `scipy/special/_ndtri_exp.pxd`

use core::f64::consts::SQRT_2;

const P1: [f64; 9] = [
    4.055_448_923_059_624_5,
    3.152_510_945_998_938_8e1,
    5.716_281_922_464_213e1,
    4.408_050_738_932_008e1,
    1.468_495_619_288_580_3e1,
    2.186_633_068_507_902_5,
    -1.402_560_791_713_545e-1,
    -3.504_246_268_278_482e-2,
    -8.574_567_851_546_854e-4,
];
const Q1: [f64; 8] = [
    1.577_998_832_564_667_5e1,
    4.539_076_351_288_792e1,
    4.131_720_382_546_72e1,
    1.504_253_856_929_075e1,
    2.504_649_462_083_094,
    -1.421_829_228_547_877_9e-1,
    -3.808_064_076_915_783e-2,
    -9.332_594_808_954_574e-4,
];
const P2: [f64; 9] = [
    3.237_748_917_769_460_3,
    6.915_228_890_689_842,
    3.938_810_252_924_744_4,
    1.333_034_608_158_075_5,
    2.014_853_895_491_790_8e-1,
    1.237_166_348_178_200_3e-2,
    3.015_815_535_082_354_3e-4,
    2.658_069_746_867_375_5e-6,
    6.239_745_391_849_833e-9,
];
const Q2: [f64; 8] = [
    6.024_270_393_647_42,
    3.679_835_638_561_608_7,
    1.377_020_994_890_813_2,
    2.162_369_935_944_966_3e-1,
    1.342_040_060_885_431_8e-2,
    3.280_144_646_821_277_4e-4,
    2.892_478_647_453_806_8e-6,
    6.790_194_080_099_813e-9,
];

/// Translated from cephes/polevl.h
#[inline]
fn polevl(x: f64, coef: &[f64], n: usize) -> f64 {
    // double ans;
    // int i;
    // const double *p;
    // p = coef;
    // ans = *p++;
    // i = N;
    // do {
    //     ans = ans * x + *p++;
    // } while (--i);
    // return (ans);

    coef.iter()
        .take(n + 1)
        .copied()
        .reduce(|acc, c| acc * x + c)
        .unwrap_or(0.0)
}

/// Translated from cephes/polevl.h
#[inline]
fn p1evl(x: f64, coef: &[f64], n: usize) -> f64 {
    // double ans;
    // const double *p;
    // int i;
    // p = coef;
    // ans = x + *p++;
    // i = N - 1;
    // do
    //     ans = ans * x + *p++;
    // while (--i);
    // return (ans);
    core::iter::once(x + coef[0])
        .chain(coef.iter().take(n + 1).skip(1).copied())
        .reduce(|acc, c| acc * x + c)
        .unwrap_or(x)
}

/// Return inverse of log CDF of normal distribution for very small y
///
/// For p sufficiently small, the inverse of the CDF of the normal
/// distribution can be approximated to high precision as a rational function
/// in sqrt(-2.0 * log(p)).
#[inline]
fn ndtri_exp_small_y(y: f64) -> f64 {
    // sqrt(-2 * y) is faster and has more precision but overflows when y < -DBL_MAX * 0.5
    let x = if y >= -f64::MAX * 0.5 {
        (-2.0 * y).sqrt()
    } else {
        SQRT_2 * (-y).sqrt()
    };
    let x0 = x - x.ln() / x;
    let z = x.recip();
    let x1 = z * if x < 8.0 {
        polevl(z, &P1, 8) / p1evl(z, &Q1, 8)
    } else {
        polevl(z, &P2, 8) / p1evl(z, &Q2, 8)
    };
    x1 - x0
}

/// Inverse of [`log_ndtr`](crate::log_ndtr) vs x.
///
/// Allows for greater precision than [`ndtri`](crate::ndtri) composed with [`f64::exp`] for very
/// small values of `y` and for `y` close to 0.
///
/// Corresponds to [`scipy.special.ndtri_exp`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri_exp.html
///
/// # See also
/// - [`log_ndtr`](crate::log_ndtr): log of the standard normal cumulative distribution function
/// - [`ndtr`](crate::ndtr): standard normal cumulative distribution function
/// - [`ndtri`](crate::ndtri): standard normal percentile function
#[must_use]
#[inline]
pub fn ndtri_exp(y: f64) -> f64 {
    if y < -f64::MAX {
        f64::NEG_INFINITY
    } else if y < -2.0 {
        ndtri_exp_small_y(y)
    } else if y > -0.145_413_457_868_859_06 {
        // y > log1p(-exp(-2))
        unsafe { -crate::ffi::xsf::ndtri(-y.exp_m1()) }
    } else {
        unsafe { crate::ffi::xsf::ndtri(y.exp()) }
    }
}

#[cfg(test)]
mod tests {
    //! Translated from `scipy.special.tests.TestNdtriExp` at
    //! <https://github.com/scipy/scipy/blob/5a7df53/scipy/special/tests/test_ndtri_exp.py>

    use core::f64;

    use crate::np_assert_allclose;

    const UNIFORM_RANDOM_POINTS: [f64; 20] = [
        0.191_519_45,
        0.622_108_77,
        0.437_727_74,
        0.785_358_58,
        0.779_975_81,
        0.272_592_61,
        0.276_464_26,
        0.801_872_18,
        0.958_139_35,
        0.875_932_63,
        0.357_817_27,
        0.500_995_13,
        0.683_462_94,
        0.712_702_03,
        0.370_250_75,
        0.561_196_19,
        0.503_083_17,
        0.013_768_45,
        0.772_826_62,
        0.882_641_19,
    ];

    fn log_ndtr_ndtri_exp(y: f64) -> f64 {
        crate::log_ndtr(crate::ndtri_exp(y))
    }

    #[test]
    fn test_very_small_arg() {
        for scale in [-1e1, -1e2, -1e10, -1e20, -f64::MAX] {
            let points = UNIFORM_RANDOM_POINTS.map(|p| scale * (0.5 * p + 0.5));
            np_assert_allclose!(points.map(log_ndtr_ndtri_exp), points, rtol = 1e-14);
        }
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_asymptotes() {
        assert_eq!(crate::ndtri_exp(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(crate::ndtri_exp(0.0), f64::INFINITY);
    }

    #[test]
    fn test_outside_domain() {
        assert!(crate::ndtri_exp(1.0).is_nan());
    }
}
