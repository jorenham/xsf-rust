use super::power::{ScaledC64, ldexp_sat, one_minus_z_pow};
use crate::ext::{Scaled, fold_exp, mul_gamma_ratio};
use crate::utils::pow2;
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::dd_real::DoubleDouble as Dd;

const U: f64 = f64::EPSILON / 2.0;
const TOL: f64 = 1e-13;

#[derive(Clone, Copy)]
struct Estimate {
    value: C64,
    abs_error: f64,
}

#[inline]
pub(super) fn evaluate(a: f64, b: f64, c: f64, z: C64, radius: f64) -> Option<C64> {
    // Skip short Maclaurin sums and angles too small to certify both components.
    // For Re z < 1, p|Im z|/(1-Re z) bounds the phases in the beta integral.
    let p = -a.min(b);
    let imag = z.im.abs();
    if p >= 12.0
        && radius <= pow2(200)
        && imag >= pow2(-100)
        && p * radius >= 6.0
        && p * imag > 0.005 * (1.0 - z.re)
    {
        contour(a.max(b), p, c, z)
    } else {
        None
    }
}

// Unlike in xsf, deform the Euler beta integral from [0, 1] through t = 1/z.
// For Im z > 0 the contour stays below the t cuts and the power 1-zt stays in
// the sector from 1 to 1-z. The arc at its zero vanishes for positive p.
// With q = 1/z and d = c-a, the first contribution is q^a B(a,p+1)/B(a,d)
// times the beta(a,p+1) mean of (1-q*u)^(d-1). Reflecting t about 1/2 gives the
// second contribution with a,d swapped and q replaced by 1-q, times (1-z)^p.
#[cold]
fn contour(a: f64, p: f64, c: f64, z: C64) -> Option<C64> {
    if !(a > 0.0 && a < c && c.is_finite()) {
        return None;
    }
    let a = Dd::new(a);
    let c = Dd::new(c);
    let d = c - a;
    // The dispatch bounds |z| between 2^-100 and 2^200, so its squared norm
    // needs no rescaling. Form (z-1)/z before rounding, especially near z=1.
    let x = Dd::new(z.re);
    let y = Dd::new(z.im);
    let yy = y * y;
    let inverse_norm = Dd::ONE / (x * x + yy);
    let imag = (y * inverse_norm).hi;
    let q = c64((x * inverse_norm).hi, -imag);
    let r = c64((((x - 1.0) * x + yy) * inverse_norm).hi, imag);
    let p1 = Dd::new(p) + 1.0;
    let left = (a, d, q);
    let right = (d, a, r);
    let left_mean = endpoint_mean(left, p1)?;
    let right_mean = endpoint_mean(right, p1)?;
    let numerator = mul_gamma_ratio([c, p1], [], || Scaled::new(1.0, 0, 0.0))?;
    let shared = (p1, numerator);
    let left = scale_endpoint(left, left_mean, shared, ScaledC64::ONE)?;
    let right = scale_endpoint(right, right_mean, shared, one_minus_z_pow(z, p))?;
    let value = left.value + right.value;
    let error =
        left.abs_error + right.abs_error + 4.0 * U * (left.value.norm() + right.value.norm());
    let magnitude = value.norm();
    // The modulus estimate alone does not protect a tiny real or imaginary part.
    let component = value.re.abs().min(value.im.abs());
    (magnitude.is_finite() && error <= TOL * magnitude && error <= 1e-12 * component)
        .then_some(value)
}

fn endpoint_mean((a, d, q): (Dd, Dd, C64), p1: Dd) -> Option<Estimate> {
    // Taylor expand (1-q*u)^(d-1) against the beta(a,p+1) weight. The integral
    // remainder is bounded by the next term's magnitude times min|1-q*u|^(d-1-n).
    // The distance from 1 to the full line through 0 and q is a lower bound.
    // Argument error can be amplified by the beta shape in q^a.
    let roundoff = (128.0 + 20.0 * a.hi) * U;
    if roundoff >= TOL {
        return None;
    }
    let radius = q.norm();
    let min_distance = q.im.abs() / radius * (1.0 - 64.0 * U);
    let mut term = C64::ONE;
    let (mut sum, mut abs_term, mut roundoff_sum) = (term, 1.0, roundoff);
    let beta_sum = p1 + a;
    for k in 0..48 {
        let k = f64::from(k);
        let n = k + 1.0;
        let exponent = (d - n).hi;
        let ratio = -exponent * (a + k).hi / (n * (beta_sum + k).hi);
        term = term * q * ratio;
        abs_term *= ratio.abs() * radius;
        if exponent <= 0.0 && abs_term <= 2.0 * U * sum.l1_norm() {
            let tail = if exponent == 0.0 {
                0.0
            } else {
                abs_term * min_distance.powf(exponent - 1.0) * (1.0 + 4096.0 * U)
            };
            if tail <= 2.0 * U * sum.norm() {
                // Bound the input and multiplication roundings of each term, all
                // additions, and the rounding of the complex prefactors.
                return Some(Estimate {
                    value: sum,
                    abs_error: roundoff_sum + tail,
                });
            }
        }
        sum += term;
        roundoff_sum += (roundoff + 40.0 * U * n) * abs_term;
    }
    None
}

fn scale_endpoint(
    (shape, other_shape, q): (Dd, Dd, C64),
    mean: Estimate,
    (p1, numerator): (Dd, Scaled),
    power: ScaledC64,
) -> Option<Estimate> {
    let beta = fold_exp(mul_gamma_ratio([], [other_shape, p1 + shape], || Some(numerator))?);
    let scale = q.powf(shape.hi) * q.powf(shape.lo) * power.mantissa * beta.m;
    let scale_abs = scale.norm();
    if !scale_abs.is_normal() {
        return None;
    }
    let mantissa = scale * mean.value;
    let error = scale_abs * (mean.abs_error + beta.err * mean.value.norm());
    if !(mantissa.is_finite() && error.is_finite()) {
        return None;
    }
    let exponent = power.exponent + f64::from(beta.e);
    Some(Estimate {
        value: ScaledC64 { mantissa, exponent }.to_c64(),
        abs_error: ldexp_sat(error, exponent).max(f64::MIN_POSITIVE),
    })
}
