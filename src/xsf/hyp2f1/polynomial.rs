use super::hypergeometric_series_generator;
use crate::utils::pow2;
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, frexp_f64, ldexp_f64, two_sum};
use crate::xsf::cephes::hyp2f1::pow_1mx;
use crate::xsf::tools::series_eval_fixed_length;

pub(super) fn evaluate(a: f64, b: f64, c: f64, z: C64, n: u64) -> C64 {
    let other = if a == -(n as f64) { b } else { a };
    // Unlike in xsf, cancel a positive common parameter on the real axis.
    if c > 0.0 && c.is_finite() && other == c && z.is_real() && z.re.is_finite() {
        let (w, w_lo) = two_sum(1.0, -z.re);
        let power = if z.re < 1.0 && w_lo != 0.0 && n > 1 {
            pow_1mx(z.re, n as f64)
        } else {
            w.powf(n as f64)
        };
        if power.is_finite() {
            return c64(power, 0.0);
        }
    }
    let count = n + 1;
    // Unlike in xsf, compare a centered polynomial when the origin sum can cancel.
    let positive_beta = other > 0.0 && other < c && c.is_finite();
    let finite_population = other < 0.0
        && other == other.trunc()
        && n as f64 - 1.0 < -c
        && (-other).max(-c) <= pow2(53);
    let value = if (finite_population || positive_beta)
        && (z.re > 0.0 || z.im != 0.0)
        // Headroom for the squared ordinate in the centered recurrence.
        && z.re.abs().max(z.im.abs()) <= pow2(256)
    {
        population(a, b, c, n, z)
    } else {
        series_eval_fixed_length(hypergeometric_series_generator(a, b, c, z), C64::ZERO, count)
    };
    // Unlike in xsf, recover a finite linear result when forming its
    // coefficient overflows for a subnormal denominator.
    if n == 1 && c.is_subnormal() && !value.is_finite() {
        let coefficient = -other;
        if coefficient.is_finite() && z.is_finite() && (coefficient / c).is_infinite() {
            let (p, ep) = frexp_f64(coefficient);
            let (q, eq) = frexp_f64(c);
            // Multiply mantissas before applying their exponents, and retain
            // the low word for cancellation against the constant term.
            let linear = |z| {
                let (w, ew) = frexp_f64(z);
                (Dd::new(p) * w / q).ldexp(ep + ew - eq)
            };
            let candidate = c64((linear(z.re) + 1.0).hi, 0.0 + linear(z.im).hi);
            if candidate.is_finite() {
                return candidate;
            }
        }
    }
    value
}

// Unlike in xsf, advance the two positive beta moments on the real interval.
fn positive_beta_real(n: u64, b: f64, c: f64, x: f64) -> Option<(C64, f64)> {
    // Keep the products and their double-word residuals away from either exponent limit.
    const MIN: f64 = pow2(-192);
    const MAX: f64 = pow2(384);
    let h = Dd::new(c) - b;
    if h.hi < MIN || b < MIN || c > pow2(256) {
        return None;
    }
    let w = Dd::ONE - x;
    let mut r = h;
    let b = Dd::new(b);
    let mut d = b;
    let c = Dd::new(c);
    let mut divisor = c;
    let mut scale = 0;
    // For T ~ Beta(b, c - b) and L = 1 - xT, r and d are proportional to
    // E[(1 - T)L^(j - 1)] and E[T L^(j - 1)] with the same scale.
    for j in 1..n {
        let j = j as f64;
        let dw = d * w;
        let sum = r + dw;
        (r, d) = (h * sum + r * j, b * sum + dw * j);
        divisor = divisor * (c + j);

        // Rescale the two moments together, independently of the normalizer.
        // Each holds at least min(b,h)/(c+j) of their total.
        let peak = r.hi.max(d.hi);
        if !(MIN..=MAX).contains(&peak) {
            let shift = if peak > MAX { -512 } else { 512 };
            let factor = pow2(shift);
            (r, d) = (r.mul_pwr2(factor), d.mul_pwr2(factor));
            scale -= shift;
        }
        if divisor.hi > MAX {
            divisor = divisor.mul_pwr2(pow2(-512));
            scale -= 512;
        }
    }
    let (m, e) = (r + w * d).frexp();
    let (divisor, divisor_exp) = divisor.frexp();
    let quotient = m / divisor;
    let exponent = scale + e - divisor_exp;
    // Positive sums propagate the largest relative input error without cancellation.
    let relative = (64.0 * n as f64 + 32.0) * pow2(-100);
    // Widen the estimate before projecting its endpoints to allow for DD rounding.
    let radius = 2.0 * quotient.hi * relative;
    let [lower, upper] = [quotient - radius, quotient + radius].map(|v| {
        let value = ldexp_f64(v.hi, exponent);
        if value <= f64::MIN_POSITIVE {
            v.ldexp(exponent + 1074).round_subnormal_units()
        } else {
            value
        }
    });
    (lower == upper && lower.is_finite()).then_some((c64(lower, 0.0), lower * relative))
}

/// Compare the origin and centered sums by their estimated absolute rounding errors.
#[cold]
fn population(a: f64, b: f64, c: f64, n: u64, z: C64) -> C64 {
    let mut terms = hypergeometric_series_generator(a, b, c, z);
    let mut mass = 0.0;
    let origin = series_eval_fixed_length(
        || {
            let term = terms();
            mass += term.l1_norm();
            term
        },
        C64::ZERO,
        n + 1,
    );
    let error = (6.0 * n as f64 + 8.0) * f64::EPSILON * mass;
    let origin_finite = origin.is_finite();
    // Keep the origin sum when the estimate leaves at least half the significand.
    if (c <= 0.0 || ((c + n as f64) * n as f64).is_finite())
        && origin_finite
        && error <= pow2(-26) * origin.re.abs().max(origin.im.abs())
    {
        return origin;
    }

    let other = if a == -(n as f64) { b } else { a };
    let (centered, rounding) = centered(n, other, c, z);
    if !centered.is_finite() {
        return origin;
    }
    if c > 0.0 {
        let select = |old: f64, new: f64, uncertainty: f64| {
            let limit = uncertainty + f64::EPSILON * new.abs() + f64::from_bits(1);
            if (old.is_finite() && (old - new).abs() > 4.0 * limit)
                || (!old.is_finite() && limit < new.abs())
            {
                new
            } else {
                old
            }
        };
        // |Im(F)| <= n |Im(z)| on the closed disk |1-z| <= 1.
        let imag_uncertainty = if z.re >= z.im.abs() && 2.0 - z.re >= z.im.abs() {
            rounding.min(centered.im.abs() + (2.0 * n as f64) * z.im.abs())
        } else {
            rounding
        };
        return c64(
            select(origin.re, centered.re, rounding),
            if z.is_real() {
                0.0
            } else {
                select(origin.im, centered.im, imag_uncertainty)
            },
        );
    }
    let magnitude = centered.re.abs().max(centered.im.abs());
    // Replace a nonfinite origin result only when the centered estimate resolves a component.
    let limit = if origin_finite { error } else { magnitude };
    if rounding + f64::EPSILON * magnitude < limit {
        centered
    } else {
        origin
    }
}

/// Polynomial in 1 - z and its absolute rounding estimate before the final projection to f64.
fn centered(n: u64, b: f64, c: f64, z: C64) -> (C64, f64) {
    // The t = 1 - z coefficients follow from Chu–Vandermonde summation (DLMF 15.4.24).
    let positive_beta = c > 0.0;
    if positive_beta && z.is_real() && (0.0..=1.0).contains(&z.re) {
        if let Some(result) = positive_beta_real(n, b, c, z.re) {
            return result;
        }
    }
    // The beta-binomial coefficients have leading term (b)_n / (c)_n.
    let product = |p: f64| {
        let shift = if p > pow2(384) { frexp_f64(p).1 } else { 0 };
        (0..n).fold((Dd::ONE, 0), |(product, scale), j| {
            let factor = if positive_beta {
                (Dd::new(p) + j as f64).ldexp(-shift)
            } else {
                // The integral falling-factor offset j is exact; p - j is exact when p is integral.
                Dd::new(p - j as f64)
            };
            let product = product * factor;
            let scale = scale + shift;
            if product.hi > pow2(384) {
                (product.ldexp(-512), scale + 512)
            } else if positive_beta && product.hi < pow2(-384) {
                let (mantissa, shift) = product.frexp();
                (mantissa, scale + shift)
            } else {
                (product, scale)
            }
        })
    };
    let (num, num_scale) = product(b.abs());
    let (den, den_scale) = product(c.abs());
    let (leading, mut scale) = (num / den, num_scale - den_scale);
    let t = Dd::ONE - z.re;
    let square = Dd::new(z.im) * z.im;
    let radius = t.hi.hypot(z.im);
    let mut coefficient = Dd::ONE;
    let (mut re, mut im) = (Dd::ONE, Dd::ZERO);
    let mut mass = 1.0;
    let h = Dd::new(c) - b;
    let excess = -h - n as f64;
    // Use one quotient correction when the real centered sum cannot cancel.
    let nonnegative = z.is_real() && z.re <= 1.0 && (positive_beta || excess >= -1.0);
    let divide = |numerator: Dd, denominator: Dd| {
        if nonnegative && denominator.lo == 0.0 {
            let d = denominator.hi;
            let q = numerator.hi / d;
            // Leave 106 exponent bits below q before the normal range ends.
            if q >= pow2(-916) {
                let residual = (-q).mul_add(d, numerator.hi) + numerator.lo;
                return Dd::new(q) + residual / d;
            }
        }
        numerator / denominator
    };
    // Small beta-binomial weights can grow back to normal, so retain their own scale.
    let mut coefficient_scale = scale;
    for j in (1..=n).rev() {
        let j = j as f64;
        let remaining = n as f64 - j;
        let term = if positive_beta {
            let (u, ue) = (h + remaining).frexp();
            let (l, le) = (Dd::new(b) + (j - 1.0)).frexp();
            let (m, e) = divide(coefficient * u * j, l * (remaining + 1.0)).frexp();
            coefficient = m;
            coefficient_scale += e + ue - le;
            let delta = coefficient_scale - scale;
            if delta > 384 {
                [re, im] = [re, im].map(|v| v.ldexp(-delta));
                mass = ldexp_f64(mass, -delta);
                scale = coefficient_scale;
                coefficient
            } else {
                coefficient.ldexp(delta)
            }
        } else {
            coefficient =
                divide(coefficient * (excess + j) * j, Dd::new(-b - j + 1.0) * (remaining + 1.0));
            coefficient
        };
        // Carry Im(P) divided by Im(t) so even a subnormal ordinate survives Horner.
        (re, im) = (re * t - im * square + term, re + im * t);
        mass = mass * radius + term.hi.abs();
        let size = mass.max(im.hi.abs());
        if size != 0.0 && !(pow2(-384)..=pow2(384)).contains(&size) {
            let (_, shift) = frexp_f64(size);
            [re, im] = [re, im].map(|v| v.ldexp(-shift));
            if !positive_beta {
                coefficient = coefficient.ldexp(-shift);
            }
            mass = ldexp_f64(mass, -shift);
            scale += shift;
        }
    }
    let round = |v: Dd, scale: i32| {
        let v = v * leading;
        let x = ldexp_f64(v.hi, scale);
        if x.abs() < f64::MIN_POSITIVE {
            v.ldexp(scale + 1074).round_subnormal_units()
        } else {
            x
        }
    };
    let imag = if z.is_real() {
        0.0
    } else {
        let (u, e) = Dd::new(-z.im).frexp();
        round(im * u, scale + e)
    };
    let rounding = ldexp_f64(mass * leading.hi * ((40.0 * n as f64 + 64.0) * pow2(-104)), scale);
    (c64(round(re, scale), imag), rounding)
}
