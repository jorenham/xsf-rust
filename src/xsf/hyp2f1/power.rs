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
    // The power at an infinite z is infinite or zero.
    if !z.is_finite() {
        return w.powf(p);
    }
    // Where (p Im z / (1 - Re z))^2 underflows, the power is linear in Im z:
    // (1-z)^p = (1-x)^p (1 - i p y / (1-x)), on the lip of the sign of y if x > 1.
    let y = z.im;
    let tiny_phase = y.abs() * p.abs().max(1.0) < f64::MIN_POSITIVE.sqrt() * w.re.abs();
    let z = if tiny_phase {
        c64(z.re, 0.0_f64.copysign(y))
    } else {
        z
    };
    let (power, exponent) = one_minus_z_pow(z, p);
    // The squarings lose about |p| 2^-104 of the power, or |p z| 2^-53 where 1-z rounds to 1.
    // Beyond the accepted error, only a real or imaginary power, or an underflow, is certain.
    let error = p.abs() * pow2(-104).min(z.l1_norm() * pow2(-53));
    let certain = power.re == 0.0 || power.im == 0.0 || exponent + 2.0 * error < -1075.0;
    if error > HYP2F1_RECUR_TOL && !certain {
        return C64::NAN;
    }
    // A nonzero imaginary part on the cut exceeds the linear term by 2^458.
    let im = if tiny_phase && power.im == 0.0 {
        let [(p, ep), (y, ey), (x, ex)] = [p, y, w.re].map(frexp_f64);
        // (a zero stays positive)
        ldexp_sat(-p * power.re * y / x, exponent + f64::from(ep + ey - ex)) + 0.0
    } else {
        ldexp_sat(power.im, exponent)
    };
    c64(ldexp_sat(power.re, exponent), im)
}

// m 2^e, where any |e| > 4096 under- or overflows.
pub(super) fn ldexp_sat(m: f64, e: f64) -> f64 {
    ldexp_f64(m, e.clamp(-4096.0, 4096.0) as i32)
}

// (1-z)^p as a complex mantissa and a binary exponent, for finite z != 1.
pub(super) fn one_minus_z_pow(z: C64, p: f64) -> (C64, f64) {
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
    (c64(real.hi, imag.hi) * w.powf(p - n), exponent)
}
