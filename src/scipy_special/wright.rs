//! Rust translation of <https://github.com/scipy/scipy/blob/d3ac276/scipy/special/wright.cc>

use crate::numpy::nextafter;
use num_complex::{Complex64, c64};
use num_traits::Zero;

const TWOITERTOL: f64 = f64::EPSILON;
const I: Complex64 = Complex64 { re: 0.0, im: 1.0 };

/// Computes fl(a+b) and err(a+b).
///
/// Translated from <https://github.com/scipy/scipy/blob/3dcd7b9/scipy/special/_round.h#L16-L25>
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let c = s - a;
    let d = b - c;
    let e = s - c;
    let err = (a - e) + d;
    (s, err)
}

/// Translated from <https://github.com/scipy/scipy/blob/3dcd7b9/scipy/special/_round.h#L28-L45>
fn add_round_up(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }

    let (s, err) = two_sum(a, b);
    if err > 0.0 {
        // fl(a + b) rounded down
        nextafter(s, f64::INFINITY)
    } else {
        // fl(a + b) rounded up or didn't round
        s
    }
}

/// Translated from <https://github.com/scipy/scipy/blob/3dcd7b9/scipy/special/_round.h#L48-L63>
fn add_round_down(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }

    let (s, err) = two_sum(a, b);
    if err < 0.0 {
        nextafter(s, f64::NEG_INFINITY)
    } else {
        s
    }
}

#[allow(clippy::too_many_lines)]
fn wrightomega_ext(z: Complex64) -> Complex64 {
    let pi = core::f64::consts::PI;
    let mut s = 1.0;

    // extract real and imaginary parts of z
    let (x, y) = (z.re, z.im);

    // compute if we are near the branch cuts
    let mut ympi = y - pi;
    let mut yppi = y + pi;
    let near = 0.1e-1;

    // Test for floating point exceptions
    if x.is_nan() || y.is_nan() {
        // NaN output for NaN inpu
        return c64(f64::NAN, f64::NAN);
    } else if x.is_infinite() && x.is_sign_negative() && (-pi < y) && (y <= pi) {
        // Signed zeros between branches
        #[allow(clippy::collapsible_else_if)]
        return if y.abs() <= pi / 2.0 {
            if y >= 0.0 {
                c64(0.0, 0.0)
            } else {
                c64(0.0, -0.0)
            }
        } else {
            if y >= 0.0 {
                c64(-0.0, 0.0)
            } else {
                c64(-0.0, -0.0)
            }
        };
    } else if x.is_infinite() || y.is_infinite() {
        // Asymptotic for large z
        return c64(x, y);
    }

    // Test if exactly on the singular points
    #[allow(clippy::float_cmp)]
    if x == -1.0 && y.abs() == pi {
        return c64(-1.0, 0.0);
    }

    // Choose approximation based on region
    let mut w: Complex64;
    #[allow(clippy::manual_range_contains)]
    if -2.0 < x && x <= 1.0 && 1.0 < y && y < 2.0 * pi {
        // Region 1: upper branch point
        // Series about z = -1 + Pi*I
        let pz = (2.0 * (z + 1.0 - I * pi)).conj().sqrt().conj();
        w = -1.0
            + (I + (1.0 / 3.0
                + (-1.0 / 36.0 * I + (1.0 / 270.0 + 1.0 / 4320.0 * I * pz) * pz) * pz)
                * pz)
                * pz;
    } else if -2.0 < x && x <= 1.0 && -2.0 * pi < y && y < -1.0 {
        // Region 2: lower branch point
        // Series about z = -1 - Pi*I
        let pz = (2.0 * (z + 1.0 + I * pi)).conj().sqrt().conj();
        w = -1.0
            + (-I
                + (1.0 / 3.0 + (1.0 / 36.0 * I + (1.0 / 270.0 - 1.0 / 4320.0 * I * pz) * pz) * pz)
                    * pz)
                * pz;
    } else if x <= -2.0 && -pi < y && y <= pi {
        // Region 3: between branch cuts
        // Series: About -infinity
        let pz = z.exp();
        w = (1.0 + (-1.0 + (3.0 / 2.0 + (-8.0 / 3.0 + 125.0 / 24.0 * pz) * pz) * pz) * pz) * pz;
        if w.is_zero() {
            // sf_error(underflow)
            return w;
        }
    } else if (-2.0 < x && x <= 1.0 && -1.0 <= y && y <= 1.0)
        || ((-2.0 < x) && (x - 0.1e1) * (x - 0.1e1) + y * y <= pi * pi)
    {
        // Region 4: Mushroom
        // Series about z=1
        let pz = z - 1.0;
        w = 1.0 / 2.0
            + 1.0 / 2.0 * z
            + (1.0 / 16.0 + (-1.0 / 192.0 + (-1.0 / 3072.0 + 13.0 / 61440.0 * pz) * pz) * pz)
                * pz
                * pz;
    } else if x <= -0.105e1 && pi < y && y - pi <= -0.75e0 * (x + 0.1e1) {
        // Region 5: Top wing
        // Negative log series
        let t = z - I * pi;
        let pz = (-t).ln();
        w = t - pz;
        let mut fac = pz / t;
        w += fac;
        fac /= t;
        w += fac * (0.5 * pz - 1.0);
        fac /= t;
        w += fac * (pz * pz / 3.0 - 3.0 * pz / 2.0 + 1.0);
        if z.norm() > 1e50 {
            // Series is accurate and the iterative scheme could overflow
            return w;
        }
    } else if x <= -0.105e1 && 0.75e0 * (x + 0.1e1) < y + pi && y + pi <= 0.0e0 {
        // Region 6: Bottom wing
        // Negative log series
        let t = z + I * pi;
        let pz = (-t).ln();
        w = t - pz;
        let mut fac = pz / t;
        w += fac;
        fac /= t;
        w += fac * (0.5 * pz - 1.0);
        fac /= t;
        w += fac * (pz * pz / 3.0 - 3.0 * pz / 2.0 + 1.0);
        if z.norm() > 1e50 {
            // Series is accurate and the iterative scheme could overflow
            return w;
        }
    } else {
        // Region 7: Everywhere else
        // Series solution about infinity
        let pz = z.ln();
        w = z - pz;
        let mut fac = pz / z;
        w += fac;
        fac /= z;
        w += fac * (0.5 * pz - 1.0);
        fac /= z;
        w += fac * (pz * pz / 3.0 - 3.0 * pz / 2.0 + 1.0);
        if z.norm() > 1e50 {
            // Series is accurate and the iterative scheme could overflow */
            return z;
        }
    }

    // Regularize if near branch cuts
    let z = if x <= -0.1e1 + near && (ympi.abs() <= near || yppi.abs() <= near) {
        s = -1.0;
        if ympi.abs() <= near {
            /* Recompute ympi with directed rounding */
            ympi = add_round_up(y, -pi);

            if ympi <= 0.0 {
                ympi = add_round_down(y, -pi);
            }

            x + I * ympi
        } else {
            // Recompute yppi with directed rounding
            yppi = add_round_up(y, pi);

            if yppi <= 0.0 {
                yppi = add_round_down(y, pi);
            }

            x + I * yppi
        }
    } else {
        z
    };

    // Iteration one
    w = s * w;
    let r = z - s * w - w.ln();
    let wp1 = s * w + 1.0;
    let e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
        / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
    w *= 1.0 + e;

    // Iteration two
    if ((2.0 * w * w - 8.0 * w - 1.0) * r.norm().powi(4)).norm()
        >= TWOITERTOL * 72.0 * wp1.norm().powi(6)
    {
        let r = z - s * w - w.ln();
        let wp1 = s * w + 1.0;
        let e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
            / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
        w *= 1.0 + e;
    }

    // Undo regularization
    s * w
}

fn wrightomega_real(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }

    // if (isinf(x)) {
    //     if (x > 0.0) {
    //         return x;
    //     } else {
    //         return 0.0;
    //     }
    // }
    if x.is_infinite() {
        // Positive infinity is asymptotically x
        // Negative infinity is zero
        return if x.is_sign_positive() { x } else { 0.0 };
    }

    // if (x < -50.0) {
    //     w = exp(x);
    //     if (w == 0.0) {
    //         sf_error("wrightomega", SF_ERROR_UNDERFLOW, "underflow in exponential series");
    //     }
    //     return w;
    // }
    if x < -50.0 {
        // Skip the iterative scheme because it exp(x) is already accurate to double precision.
        let w = x.exp();
        if w == 0.0 {
            // underflow in exponential series
        }
        return w;
    }
    // if (x > 1e20) {
    //     return x;
    // }
    if x > 1e20 {
        // Skip the iterative scheme because the result is just x to double precision
        return x;
    }

    // Split into three distinct intervals (-inf,-2), [-2,1), [1,inf)
    // if (x < -2.0) {
    //     w = exp(x);
    // } else if (x < 1) {
    //     w = exp(2.0*(x-1.0)/3.0);
    // } else {
    //     w = log(x);
    //     w = x - w + w/x;
    // }
    let w = if x < -2.0 {
        // exponential is approx < 1.3-1 accurate
        x.exp()
    } else if x < 1.0 {
        // on [-2,1) approx < 1.5e-1 accurate
        (2.0 * (x - 1.0) / 3.0).exp()
    } else {
        // infinite series with 2 terms approx <1.7e-1 accurate
        let w = x.ln();
        x - w + w / x
    };

    // Iteration one of Fritsch, Shafer, and Crowley (FSC) iteration
    // r = x - w - log(w);
    let r = x - w - w.ln();
    // wp1 = w + 1.0;
    let wp1 = w + 1.0;
    // e = r/wp1*(2.0*wp1*(wp1+2.0/3.0*r)-r)/(2.0*wp1*(wp1+2.0/3.0*r)-2.0*r);
    let e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
        / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
    // w = w*(1.0+e);
    let w = w * (1.0 + e);

    // Iteration two (if needed based on the condition number)
    // if (fabs((2.0*w*w-8.0*w-1.0)*pow(fabs(r),4.0)) >= TWOITERTOL*72.0*pow(fabs(wp1),6.0)) {
    //     r = x - w - log(w);
    //     wp1 = w+1.0;
    //     e = r/wp1*(2.0*wp1*(wp1+2.0/3.0*r)-r)/(2.0*wp1*(wp1+2.0/3.0*r)-2.0*r);
    //     w = w*(1.0+e);
    // }
    // return w;
    if ((2.0 * w * w - 8.0 * w - 1.0) * r.abs().powi(4)).abs()
        >= TWOITERTOL * 72.0 * wp1.abs().powi(6)
    {
        let r = x - w - w.ln();
        let wp1 = w + 1.0;
        let e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
            / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
        w * (1.0 + e)
    } else {
        w
    }
}

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait WrightOmegaArg: sealed::Sealed {
    fn wrightomega(self) -> Self;
}

impl WrightOmegaArg for f64 {
    #[inline]
    fn wrightomega(self) -> Self {
        wrightomega_real(self)
    }
}

impl WrightOmegaArg for Complex64 {
    #[inline]
    fn wrightomega(self) -> Self {
        wrightomega_ext(self)
    }
}

/// Wright Omega function, $\omega(z)$
///
/// Translated into pure Rust from [`scipy.special.wrightomega`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wrightomega.html
///
/// # Definition
///
/// Defined as the solution to
///
/// $$
/// \omega + \ln(\omega) = z \\ ,
/// $$
///
/// where $\ln$ is the principal branch of the complex logarithm.
///
/// # See also
/// - [`lambertw`](crate::lambertw): Lambert W function
#[must_use]
#[inline]
pub fn wrightomega<T: WrightOmegaArg>(z: T) -> T {
    z.wrightomega()
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    //! based on <https://github.com/scipy/scipy/blob/main/scipy/special/tests/test_wrightomega.py>

    use core::f64::{self, consts::PI};
    use num_complex::c64;
    use xsref::np_assert_allclose;

    #[test]
    fn test_wrightomega_nan() {
        let pts = [
            c64(f64::NAN, 0.0),
            c64(0.0, f64::NAN),
            c64(f64::NAN, f64::NAN),
            c64(f64::NAN, 1.0),
            c64(1.0, f64::NAN),
        ];
        for p in pts {
            let res = crate::wrightomega(p);
            assert!(res.re.is_nan());
            assert!(res.im.is_nan());
        }
    }

    #[test]
    fn test_wrightomega_inf_branch() {
        let pts = [
            c64(f64::NEG_INFINITY, PI / 4.0),
            c64(f64::NEG_INFINITY, -PI / 4.0),
            c64(f64::NEG_INFINITY, 3.0 * PI / 4.0),
            c64(f64::NEG_INFINITY, -3.0 * PI / 4.0),
        ];
        let expected_results = [
            c64(0.0, 0.0),
            c64(0.0, -0.0),
            c64(-0.0, 0.0),
            c64(-0.0, -0.0),
        ];
        for i in 0..4 {
            let (p, expected) = (pts[i], expected_results[i]);
            let res = crate::wrightomega(p);
            assert_eq!(res.re.signum(), expected.re.signum());
            assert_eq!(res.im.signum(), expected.im.signum());
        }
    }

    #[test]
    fn test_wrightomega_inf() {
        let pts = [
            c64(f64::INFINITY, 10.0),
            c64(f64::NEG_INFINITY, 10.0),
            c64(10.0, f64::INFINITY),
            c64(10.0, f64::NEG_INFINITY),
        ];
        for p in pts {
            assert_eq!(crate::wrightomega(p), p);
        }
    }

    #[test]
    fn test_wrightomega_singular() {
        let pts = [c64(-1.0, PI), c64(-1.0, -PI)];
        for p in pts {
            let res = crate::wrightomega(p);
            assert_eq!(res, c64(-1.0, 0.0));
            assert!(res.im.is_sign_positive());
        }
    }

    #[test]
    fn test_wrightomega_real_infinities() {
        assert_eq!(crate::wrightomega(f64::NEG_INFINITY), 0.0);
        assert_eq!(crate::wrightomega(f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn test_wrightomega_real_nan() {
        assert!(crate::wrightomega(f64::NAN).is_nan());
    }

    #[test]
    fn test_wrightomega_real_series_crossover() {
        let desired_error = 2.0 * f64::EPSILON;
        let crossover = 1e20;
        let x_before_crossover = super::nextafter(crossover, f64::NEG_INFINITY);
        let x_after_crossover = super::nextafter(crossover, f64::INFINITY);
        // Computed using Mpmath
        let desired_before_crossover = 99_999_999_999_999_983_569.948;
        let desired_after_crossover = 100_000_000_000_000_016_337.948;
        np_assert_allclose!(
            [crate::wrightomega(x_before_crossover)],
            [desired_before_crossover],
            atol = 0.0,
            rtol = desired_error
        );
        np_assert_allclose!(
            [crate::wrightomega(x_after_crossover)],
            [desired_after_crossover],
            atol = 0.0,
            rtol = desired_error
        );
    }

    #[test]
    fn test_wrightomega_exp_approximation_crossover() {
        let desired_error = 2.0 * f64::EPSILON;
        let crossover = -50.0;
        let x_before_crossover = super::nextafter(crossover, f64::INFINITY);
        let x_after_crossover = super::nextafter(crossover, f64::NEG_INFINITY);
        // Computed using Mpmath
        let desired_before_crossover = 1.928_749_847_963_931_5e-22;
        let desired_after_crossover = 1.928_749_847_963_904_2e-22;
        np_assert_allclose!(
            [crate::wrightomega(x_before_crossover)],
            [desired_before_crossover],
            atol = 0.0,
            rtol = desired_error
        );
        np_assert_allclose!(
            [crate::wrightomega(x_after_crossover)],
            [desired_after_crossover],
            atol = 0.0,
            rtol = desired_error
        );
    }

    #[test]
    fn test_wrightomega_real_versus_complex() {
        // x = np.linspace(-500, 500, 1001)
        let res_f: Vec<f64> = (-500..=500)
            .map(|x| crate::wrightomega(f64::from(x)))
            .collect();
        let res_c: Vec<f64> = (-500..=500)
            .map(|x| crate::wrightomega(c64(f64::from(x), 0.0)).re)
            .collect();
        np_assert_allclose!(res_c, res_f, atol = 0.0, rtol = 1e-14);
    }
}
