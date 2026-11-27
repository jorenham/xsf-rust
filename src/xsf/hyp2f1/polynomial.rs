use super::power::{ldexp_sat, one_minus_z_pow};
use super::{Estimate, HYP2F1_MAX_FIXED_TERMS, hypergeometric_series_generator};
use crate::ext::{Scaled, fold_exp, hyp2f1_poly_dd};
use crate::utils::{DENORM_MIN, MIN_SUBNORMAL_EXP, pow2};
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::dd_real::{DoubleDouble as Dd, frexp_f64, ldexp_f64, two_sum};
use crate::xsf::cephes::hyp2f1::HYP2F1_ETHRESH;
use crate::xsf::tools::series_eval_fixed_length;

/// The double-double polynomial as an `f64` if `err` is small, and the interval of its error
/// bound if that holds more than one rounding of a (nearly) subnormal value
fn dd_candidate(y: Scaled) -> Option<(f64, Option<[f64; 2]>)> {
    let y = fold_exp(y);
    let [lo, v, hi] = [1.0 - y.err, 1.0, 1.0 + y.err].map(|f| ldexp_f64(y.m * f, y.e));
    let open = lo != hi && lo.abs().min(hi.abs()) < f64::MIN_POSITIVE;
    (v.is_finite() && y.err <= HYP2F1_ETHRESH).then_some((v, open.then_some([lo, hi])))
}

/// The series terminated by the negative integer `a` or `b`, or NaN if it is too long.
pub(super) fn terminating(a: f64, b: f64, c: f64, z: C64) -> Estimate {
    let neg_int = |p: f64| p == p.trunc() && p < 0.0;
    let degree = match (neg_int(a), neg_int(b)) {
        (true, true) => -a.max(b),
        (true, false) => -a,
        _ => -b,
    };
    // Unlike in xsf, cap the number of terms in the finite series.
    if degree < HYP2F1_MAX_FIXED_TERMS as f64 {
        evaluate(a, b, c, z, degree as u64)
    } else {
        Estimate::from(C64::NAN)
    }
}

// The polynomial and the estimated absolute error of an origin or Pfaff sum (zero for other
// forms).
fn evaluate(a: f64, b: f64, c: f64, z: C64, n: u64) -> Estimate {
    let other = if a == -(n as f64) { b } else { a };
    // Unlike in xsf, evaluate the uniform beta-binomial polynomial as a geometric sum.
    if other == 1.0 && c == 2.0 && z.is_finite() {
        let mean = uniform_beta_binomial(n + 1, z);
        if mean.is_finite() {
            return Estimate::from(mean);
        }
    }
    // Unlike in xsf, use the double-double polynomial on the real axis (at z = 1 the
    // Chu-Vandermonde product, DLMF 15.4.24), as the terms of the series can cancel, or the
    // series if it lies within an open subnormal rounding of it.
    if z.is_real() {
        match hyp2f1_poly_dd(a, b, c, z.re).and_then(dd_candidate) {
            Some((y, None)) => return Estimate::from(c64(y, 0.0)),
            Some((y, Some([lo, hi]))) => {
                let value = series(a, b, c, z, n, other).value;
                let inside = (lo.min(hi)..=lo.max(hi)).contains(&value.re) && value.im == 0.0;
                return Estimate::from(if inside { value } else { c64(y, 0.0) });
            },
            None => {},
        }
    }
    series(a, b, c, z, n, other)
}

// The series as in `evaluate`, with `other` the upper parameter that is not `-n`.
fn series(a: f64, b: f64, c: f64, z: C64, n: u64, other: f64) -> Estimate {
    // Unlike in xsf, try alternate forms when the origin sum cancels: the centered polynomial for
    // a positive beta or finite population, else the Pfaff transformation and then the centered
    // polynomial.
    let positive_beta = other > 0.0 && other < c && c.is_finite();
    let finite_population = other < 0.0
        && other == other.trunc()
        && n as f64 - 1.0 < -c
        && (-other).max(-c) <= pow2(53);
    let centerable = (z.re > 0.0 || z.im != 0.0)
        // Headroom for the squared ordinate in the centered recurrence.
        && z.linf_norm() <= pow2(256);
    let estimate = if (finite_population || positive_beta) && centerable {
        Estimate::from(population(a, b, c, n, z))
    } else {
        // Keep a linear sum, which keeps exact roots, and real terms of one sign, which cannot
        // cancel.
        let origin = estimated_sum(a, b, c, n, z, 0.0);
        if n == 1 || one_sign(other, c, n, z) || origin.accurate() {
            origin
        } else {
            let pfaff = pfaff(two_sum(c, -other), c, n, z, origin);
            (centerable && !pfaff.accurate())
                .then(|| centered(n, other, c, z))
                .filter(|centered| centered.resolved() && !(centered.error >= pfaff.error))
                .unwrap_or(pfaff)
        }
    };
    // Unlike in xsf, recover a finite linear result when forming its
    // coefficient overflows for a subnormal denominator.
    if n == 1 && c.is_subnormal() && !estimate.value.is_finite() {
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
                return Estimate::from(candidate);
            }
        }
    }
    estimate
}

// Mean of 1, w, ..., w^(count-1), with w = 1-z.
// Store each complex number as a + i*b*z.im to retain subnormal input ordinates.
fn uniform_beta_binomial(count: u64, z: C64) -> C64 {
    let multiply =
        |(a, b): (Dd, Dd), (c, d): (Dd, Dd)| (a * c - (b * z.im) * (d * z.im), a * d + b * c);
    let w = (Dd::ONE - z.re, -Dd::ONE);
    let mut power = w;
    let mut sum = (Dd::ONE, Dd::new(0.0));
    for bit in (0..count.ilog2()).rev() {
        // S(2k) = S(k) (1 + w^k), and S(2k+1) = S(2k) + w^(2k).
        sum = multiply(sum, (power.0 + 1.0, power.1));
        power = multiply(power, power);
        if count >> bit & 1 != 0 {
            sum = (sum.0 + power.0, sum.1 + power.1);
            power = multiply(power, w);
        }
    }
    c64((sum.0 / count as f64).hi, ((sum.1 / count as f64) * z.im).hi)
}

// Unlike in xsf, advance the two positive beta moments on the real interval.
fn positive_beta_real(n: u64, b: f64, c: f64, x: f64) -> Option<Estimate> {
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
            v.ldexp(exponent - MIN_SUBNORMAL_EXP)
                .round_subnormal_units()
        } else {
            value
        }
    });
    (lower == upper && lower.is_finite()).then_some(Estimate {
        value: c64(lower, 0.0),
        error: lower * relative,
    })
}

// Whether the real terms share a sign: the factors of the term ratio
// (k - n)(b + k)z / ((k + 1)(c + k)) keep their signs for 0 <= k < n, and its sign is positive.
fn one_sign(b: f64, c: f64, n: u64, z: C64) -> bool {
    let last = (n - 1) as f64;
    let fixed = |p: f64| p >= 0.0 || p + last <= 0.0;
    z.im == 0.0 && fixed(b) && fixed(c) && ((b < 0.0) == (c < 0.0)) == (z.re < 0.0)
}

// The sum of the terms up to z^n, and its estimated absolute error, for terms with an extra
// relative uncertainty. The estimate is infinite where a coefficient denominator can overflow
// and silently drop a term.
fn estimated_sum(a: f64, b: f64, c: f64, n: u64, z: C64, relative: f64) -> Estimate {
    let mut terms = hypergeometric_series_generator(a, b, c, z);
    let mut mass = 0.0;
    let value = series_eval_fixed_length(
        || {
            let term = terms();
            mass += term.l1_norm();
            term
        },
        C64::ZERO,
        n + 1,
    );
    let n = n as f64;
    let error = if (n * (c.abs() + n)).is_finite() {
        ((6.0 * n + 8.0) * f64::EPSILON + relative) * mass
    } else {
        f64::INFINITY
    };
    Estimate { value, error }
}

/// `F(-n, c - a; c; z)`, given its `estimate` from the rounded `c - a`.
pub(super) fn complement(a: f64, c: f64, n: u64, z: C64, estimate: Estimate) -> Estimate {
    // Where c - a rounds by more than the ulp scale of a, the double-double polynomial on the real
    // axis, or else the Pfaff form, uses the exact c - a instead of the rounded one.
    let (c_minus_a, rounding) = two_sum(c, -a);
    if rounding.abs() <= f64::EPSILON * a.abs() || n >= HYP2F1_MAX_FIXED_TERMS {
        return estimate;
    }
    let degree = n as f64;
    let dd = z
        .is_real()
        .then(|| hyp2f1_poly_dd(-degree, Dd::new(c) - a, c, z.re).and_then(dd_candidate))
        .flatten();
    if let Some((y, _)) = dd {
        // (keeping an estimate that agrees within the accepted double-double error)
        return if (y - estimate.value.re).abs() <= HYP2F1_ETHRESH * y.abs() {
            estimate
        } else {
            Estimate::from(c64(y, 0.0))
        };
    }
    let relative = rounding_change(c_minus_a, rounding, n);
    let rounded = estimated_sum(-degree, c_minus_a, c, n, z, relative);
    let error = estimate.error.max(rounded.error);
    pfaff((a, 0.0), c, n, z, Estimate { error, ..estimate })
}

// F(-n, b; c; z) = (1 - z)^n F(-n, c - b; c; z / (z - 1)) (DLMF 15.8.1), where the tuple is c - b
// rounded and its rounding error
#[cold]
fn pfaff((c_minus_b, rounding): (f64, f64), c: f64, n: u64, z: C64, origin: Estimate) -> Estimate {
    let x = z / (z - 1.0);
    if !x.is_finite() {
        return origin;
    }
    let degree = n as f64;
    let relative = rounding_change(c_minus_b, rounding, n);
    let sum = estimated_sum(-degree, c_minus_b, c, n, x, relative);
    let power = one_minus_z_pow(z, degree);
    let value = (power * sum.value).to_c64();
    // The power and the product each round a component by at most eps |power| |sum|, and the
    // final scaling by less than the spacing of the subnormals.
    let product_error = 2.0 * f64::EPSILON * sum.value.norm();
    let pfaff_error =
        ldexp_sat(power.mantissa.norm() * (sum.error + product_error), power.exponent) + DENORM_MIN;
    // Take each finite component that the Pfaff estimate resolves, and that differs from a
    // less certain origin one by more than its error.
    let select = |new: f64, old: f64| {
        let differs = pfaff_error < origin.error && (new - old).abs() > pfaff_error;
        if new.is_finite() && pfaff_error < new.abs() && (differs || !old.is_finite()) {
            (new, pfaff_error)
        } else {
            (old, origin.error)
        }
    };
    let (re, re_error) = select(value.re, origin.value.re);
    let (im, im_error) = select(value.im, origin.value.im);
    Estimate {
        value: c64(re, im),
        error: re_error.max(im_error),
    }
}

// Rounding p changes each term of F(-n, p; c; z) by a relative n |rounding| / min |p + j| for
// 0 <= j < n to first order, and by at most twice that while it is below one.
fn rounding_change(p: f64, rounding: f64, n: u64) -> f64 {
    let degree = n as f64;
    let change = degree * rounding.abs();
    let distance = (p + (-p).round().clamp(0.0, degree - 1.0)).abs();
    if change == 0.0 {
        0.0
    } else if change <= distance {
        2.0 * change / distance
    } else {
        f64::INFINITY
    }
}

/// Compare the origin and centered sums by their estimated absolute rounding errors.
#[cold]
fn population(a: f64, b: f64, c: f64, n: u64, z: C64) -> C64 {
    let origin = estimated_sum(a, b, c, n, z, 0.0);
    if origin.accurate() {
        return origin.value;
    }
    let other = if a == -(n as f64) { b } else { a };
    let estimate = centered(n, other, c, z);
    let centered = estimate.value;
    if !centered.is_finite() {
        return origin.value;
    }
    if c > 0.0 {
        return select_positive_beta(origin.value, &estimate, n, z);
    }
    let magnitude = centered.linf_norm();
    // Replace a nonfinite origin result only when the centered estimate resolves a component.
    let limit = if origin.value.is_finite() {
        origin.error
    } else {
        magnitude
    };
    if estimate.error + f64::EPSILON * magnitude < limit {
        centered
    } else {
        origin.value
    }
}

fn select_positive_beta(origin: C64, centered: &Estimate, n: u64, z: C64) -> C64 {
    let select = |old: f64, new: f64, uncertainty: f64| {
        let limit = uncertainty + f64::EPSILON * new.abs() + DENORM_MIN;
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
        centered
            .error
            .min(centered.value.im.abs() + (2.0 * n as f64) * z.im.abs())
    } else {
        centered.error
    };
    c64(
        select(origin.re, centered.value.re, centered.error),
        if z.is_real() {
            0.0
        } else {
            select(origin.im, centered.value.im, imag_uncertainty)
        },
    )
}

// The leading coefficient is (b)_n/(c)_n, from the rising factors p + j of each parameter p, or
// else the falling factors -p - j. Carry the products' binary scales separately.
fn leading_coefficient(n: u64, b: f64, c: f64, rising: bool) -> (Dd, i32) {
    let product = |p: f64| {
        let shift = if p.abs() > pow2(384) {
            frexp_f64(p).1
        } else {
            0
        };
        (0..n).fold((Dd::ONE, 0), |(product, scale), j| {
            let factor = if rising {
                (Dd::new(p) + j as f64).ldexp(-shift)
            } else {
                // The integral offset j is exact; -p - j is exact when p is integral.
                Dd::new(-p - j as f64)
            };
            let product = product * factor;
            let scale = scale + shift;
            if product.hi.abs() > pow2(384) {
                (product.ldexp(-512), scale + 512)
            } else if rising && product.hi.abs() < pow2(-384) {
                let (mantissa, shift) = product.frexp();
                (mantissa, scale + shift)
            } else {
                (product, scale)
            }
        })
    };
    let (num, num_scale) = product(b);
    let (den, den_scale) = product(c);
    (num / den, num_scale - den_scale)
}

/// Polynomial in 1 - z and its absolute rounding estimate before the final projection to f64.
fn centered(n: u64, b: f64, c: f64, z: C64) -> Estimate {
    // The t = 1 - z coefficients follow from Chu–Vandermonde summation (DLMF 15.4.24).
    if c > 0.0 && z.is_real() && (0.0..=1.0).contains(&z.re) {
        if let Some(result) = positive_beta_real(n, b, c, z.re) {
            return result;
        }
    }
    // Use falling factors, with exact offsets, for an integral b < 0 and c < 0.
    let rising = c > 0.0 || b > 0.0 || b != b.trunc();
    let (leading, mut scale) = leading_coefficient(n, b, c, rising);
    let t = Dd::ONE - z.re;
    let square = Dd::new(z.im) * z.im;
    let radius = t.hi.hypot(z.im);
    let mut coefficient = Dd::ONE;
    let (mut re, mut im) = (Dd::ONE, Dd::ZERO);
    let mut mass = 1.0;
    let h = Dd::new(c) - b;
    let excess = -h - n as f64;
    // Use one quotient correction when the real centered sum cannot cancel.
    let nonnegative = z.is_real() && z.re <= 1.0 && ((0.0 < b && b < c) || excess >= -1.0);
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
        let term = if rising {
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
            if !rising {
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
            v.ldexp(scale - MIN_SUBNORMAL_EXP).round_subnormal_units()
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
    let rounding =
        ldexp_f64(mass * leading.hi.abs() * ((40.0 * n as f64 + 64.0) * pow2(-104)), scale);
    Estimate {
        value: c64(round(re, scale), imag),
        error: rounding,
    }
}
