use core::ops::Mul;

use crate::ext::HYP2F1_RECUR_TOL;
use crate::utils::pow2;
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, frexp_f64, ldexp_f64};

// (1-z)^p, on the lip of the cut that the sign of a zero Im z selects.
pub(super) fn pow_1mz(z: C64, p: f64) -> C64 {
    let w = 1.0 - z;
    if w == C64::ZERO {
        return c64(0.0_f64.powf(p), 0.0);
    }
    // The power at an infinite or NaN z is infinite, zero or NaN.
    if !z.is_finite() {
        return w.powf(p);
    }
    // Unlike in xsf, where the square of p times the ratio of the components of 1-z underflows,
    // take the power of the larger one, on its lip of the cut, times 1 + i p t for the tilt t.
    // (compared as mantissas and exponents, as the scaled bound itself underflows)
    let [(pm, pe), (re, re_e), (im, im_e)] =
        [p.abs().max(1.0), w.re.abs(), w.im.abs()].map(frexp_f64);
    let below = |small: f64, small_e: i32, large: f64, large_e: i32| {
        pm * small < ldexp_f64(large, large_e - small_e - pe + (f64::MIN_EXP - 1) / 2)
    };
    let (axis, tilt) = if below(im, im_e, re, re_e) {
        (c64(z.re, 0.0_f64.copysign(z.im)), Some((w.im, w.re)))
    } else if below(re, re_e, im, im_e) {
        (c64(1.0, z.im), Some((-w.re, w.im)))
    } else {
        (z, None)
    };
    let scaled = one_minus_z_pow_sqrt(axis, p);
    let ScaledC64 { mantissa, exponent } = scaled;
    // The squarings lose about |p| 2^-104 of the power, or |p z| 2^-53 where 1-z rounds to 1.
    // Unlike in xsf, beyond the accepted error only a real or imaginary power, or an
    // underflow, is certain: the rest is NaN.
    let error = p.abs() * pow2(-104).min(axis.l1_norm() * pow2(-53));
    let certain = mantissa.re == 0.0 || mantissa.im == 0.0 || exponent + 2.0 * error < -1075.0;
    if error > HYP2F1_RECUR_TOL && !certain {
        return C64::NAN;
    }
    let mut value = scaled.to_c64();
    // The linear term only reaches a component that vanishes on the axis.
    if let Some((num, den)) = tilt {
        let linear = |component: f64| {
            let [(p, ep), (num, en), (den, ed)] = [p, num, den].map(frexp_f64);
            ldexp_sat(p * component * num / den, exponent + f64::from(ep + en - ed))
        };
        // (a zero stays positive)
        if mantissa.im == 0.0 {
            value.im = linear(mantissa.re) + 0.0;
        } else if mantissa.re == 0.0 {
            value.re = 0.0 - linear(mantissa.im);
        }
    }
    value
}

// m 2^e, where any |e| > 4096 under- or overflows.
pub(super) fn ldexp_sat(m: f64, e: f64) -> f64 {
    ldexp_f64(m, e.clamp(-4096.0, 4096.0) as i32)
}

// mantissa 2^exponent, where the exponent can exceed i32.
#[derive(Clone, Copy)]
pub(super) struct ScaledC64 {
    pub(super) mantissa: C64,
    pub(super) exponent: f64,
}

impl ScaledC64 {
    pub(super) const ONE: Self = Self {
        mantissa: C64::ONE,
        exponent: 0.0,
    };

    pub(super) const ZERO: Self = Self {
        mantissa: C64::ZERO,
        exponent: 0.0,
    };

    pub(super) fn to_c64(self) -> C64 {
        let Self { mantissa, exponent } = self;
        c64(ldexp_sat(mantissa.re, exponent), ldexp_sat(mantissa.im, exponent))
    }

    // The same value, exactly, with the larger component of a finite mantissa in [0.5, 1).
    pub(super) fn normalized(self) -> Self {
        let Self { mantissa, exponent } = self;
        let (_, shift) = frexp_f64(mantissa.linf_norm());
        Self {
            mantissa: c64(ldexp_f64(mantissa.re, -shift), ldexp_f64(mantissa.im, -shift)),
            exponent: exponent + f64::from(shift),
        }
    }
}

impl Mul for ScaledC64 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            mantissa: self.mantissa * rhs.mantissa,
            exponent: self.exponent + rhs.exponent,
        }
    }
}

impl Mul<C64> for ScaledC64 {
    type Output = Self;

    // Scale the factor by its larger component, which then stays in range until the final
    // scaling.
    fn mul(self, rhs: C64) -> Self {
        self * Self {
            mantissa: rhs,
            exponent: 0.0,
        }
        .normalized()
    }
}

// (-z)^p, for |z| > 1, with a separate scale.
// (1 - (z + 1) is -z to about an ulp of |z|.)
pub(super) fn neg_z_pow(z: C64, p: f64) -> ScaledC64 {
    let power = ScaledC64 {
        mantissa: (-z).powf(p),
        exponent: 0.0,
    };
    if !z.is_finite() {
        power
    } else if power.mantissa.norm().is_normal() {
        power.normalized()
    } else {
        one_minus_z_pow(z + 1.0, p)
    }
}

// (1-z)^p, for finite z != 1.
pub(super) fn one_minus_z_pow_sqrt(z: C64, p: f64) -> ScaledC64 {
    // Unlike in xsf, a half power is a square root, which keeps a component that vanishes on
    // the cut.
    let half = p - p.round();
    if half.abs() == 0.5 {
        let power = one_minus_z_pow(z, p.round());
        let root = c64(1.0 - z.re, -z.im).sqrt();
        let mantissa = power.mantissa * if half > 0.0 { root } else { 1.0 / root };
        ScaledC64 { mantissa, ..power }
    } else {
        one_minus_z_pow(z, p)
    }
}

// (1-z)^p, for finite z != 1, or NaN for a nonfinite p.
pub(super) fn one_minus_z_pow(z: C64, p: f64) -> ScaledC64 {
    // (halving a nonfinite degree never reaches zero)
    if !p.is_finite() {
        return ScaledC64 {
            mantissa: C64::NAN,
            exponent: 0.0,
        };
    }
    // Exponentiation by squaring with double-double components and a separate scale.
    let n = p.round();
    let real = Dd::ONE - z.re;
    let w = c64(real.hi, -z.im);
    let normalize = |re: Dd, im: Dd, exponent: f64| {
        let (_, shift) = frexp_f64(re.hi.abs().max(im.hi.abs()));
        (re.ldexp(-shift), im.ldexp(-shift), exponent + f64::from(shift))
    };
    let multiply = |(a, b, i): (Dd, Dd, f64), (c, d, j): (Dd, Dd, f64)| {
        normalize(a * c - b * d, a * d + b * c, i + j)
    };
    let mut base = normalize(real, Dd::new(w.im), 0.0);
    let mut value = None;
    let mut degree = n.abs();
    while degree != 0.0 {
        let half = (degree / 2.0).floor();
        if degree != 2.0 * half {
            value = Some(value.map_or(base, |value| multiply(value, base)));
        }
        degree = half;
        if degree != 0.0 {
            base = multiply(base, base);
        }
    }
    let power = value.unwrap_or((Dd::ONE, Dd::ZERO, 0.0));
    let (real, imag, exponent) = if n < 0.0 {
        let (re, im, exponent) = power;
        let norm = re * re + im * im;
        (re / norm, -im / norm, -exponent)
    } else {
        power
    };
    // The remaining exponent is small and barely amplifies rounding in 1-z.
    let mantissa = c64(real.hi, imag.hi) * w.powf(p - n);
    ScaledC64 { mantissa, exponent }
}
