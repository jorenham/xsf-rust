use num_complex::Complex;

/// Complex division `a / b` using Smith's algorithm, which avoids the overflow and underflow of
/// `|b|^2` in the naive formula. Tiny operands are first scaled by the same power of 2, so that
/// the intermediate results don't lose precision by becoming subnormal.
///
/// This replaces the complex division in C++, i.e. GCC's `__divdc3`. Like `__divdc3`, it recovers
/// the infinities and zeros that would otherwise be computed as NaN + NaN i (C99 Annex G), but
/// unlike `__divdc3`, it doesn't scale huge operands (e.g. `1 / (MAX + i MAX)` underflows to zero).
///
/// R. L. Smith, "Algorithm 116: Complex division", Communications of the ACM 5(8), 1962.
#[inline]
pub(crate) fn cdiv(mut a: Complex<f64>, mut b: Complex<f64>) -> Complex<f64> {
    const TWO_M400: f64 = f64::from_bits((1023 - 400) << 52); // 2^-400
    const TWO_600: f64 = f64::from_bits((1023 + 600) << 52); // 2^600

    if b.re.abs().max(b.im.abs()) < TWO_M400 {
        a *= TWO_600;
        b *= TWO_600;
    }
    let w = if b.re.abs() >= b.im.abs() {
        let r = b.im / b.re;
        let den = b.re + b.im * r;
        Complex::new((a.re + a.im * r) / den, (a.im - a.re * r) / den)
    } else {
        let r = b.re / b.im;
        let den = b.re * r + b.im;
        Complex::new((a.re * r + a.im) / den, (a.im * r - a.re) / den)
    };

    // Recover infinities and zeros that computed as NaN+iNaN; the only
    // cases are nonzero/zero, infinite/finite, and finite/infinite.
    if w.re.is_nan() && w.im.is_nan() {
        let (Complex { re: a, im: b }, Complex { re: c, im: d }) = (a, b);
        if c == 0.0 && d == 0.0 && (!a.is_nan() || !b.is_nan()) {
            let inf = f64::INFINITY.copysign(c);
            return Complex::new(inf * a, inf * b);
        } else if (a.is_infinite() || b.is_infinite()) && c.is_finite() && d.is_finite() {
            let a = box_inf(a);
            let b = box_inf(b);
            return Complex::new(
                f64::INFINITY * (a * c + b * d),
                f64::INFINITY * (b * c - a * d),
            );
        } else if (c.is_infinite() || d.is_infinite()) && a.is_finite() && b.is_finite() {
            let c = box_inf(c);
            let d = box_inf(d);
            return Complex::new(0.0 * (a * c + b * d), 0.0 * (b * c - a * d));
        }
    }
    w
}

/// Complex multiplication `a * b`, like GCC's `__muldc3` (C99 Annex G)
///
/// This replaces the complex multiplication in C++. For finite results it is the same as the naive
/// formula (i.e. `num_complex`'s `*`), but it recovers the infinities that would otherwise be
/// computed as NaN + NaN i, e.g. `(inf + inf i) * (1 + 0i)`.
#[inline]
pub(crate) fn cmul(a: Complex<f64>, b: Complex<f64>) -> Complex<f64> {
    let (
        Complex {
            re: mut a,
            im: mut b,
        },
        Complex {
            re: mut c,
            im: mut d,
        },
    ) = (a, b);
    let (ac, bd, ad, bc) = (a * c, b * d, a * d, b * c);
    let (mut x, mut y) = (ac - bd, ad + bc);

    if x.is_nan() && y.is_nan() {
        // Recover infinities that computed as NaN+iNaN.
        let mut recalc = false;
        if a.is_infinite() || b.is_infinite() {
            // z is infinite.  "Box" the infinity and change NaNs in
            // the other factor to 0.
            a = box_inf(a);
            b = box_inf(b);
            c = nan_to_zero(c);
            d = nan_to_zero(d);
            recalc = true;
        }
        if c.is_infinite() || d.is_infinite() {
            // w is infinite.  "Box" the infinity and change NaNs in
            // the other factor to 0.
            c = box_inf(c);
            d = box_inf(d);
            a = nan_to_zero(a);
            b = nan_to_zero(b);
            recalc = true;
        }
        if !recalc && (ac.is_infinite() || bd.is_infinite() || ad.is_infinite() || bc.is_infinite())
        {
            // Recover infinities from overflow by changing NaNs to 0.
            a = nan_to_zero(a);
            b = nan_to_zero(b);
            c = nan_to_zero(c);
            d = nan_to_zero(d);
            recalc = true;
        }
        if recalc {
            x = f64::INFINITY * (a * c - b * d);
            y = f64::INFINITY * (a * d + b * c);
        }
    }
    Complex::new(x, y)
}

/// `isnan(x) ? copysign(0, x) : x`, used by [`cmul`]
fn nan_to_zero(x: f64) -> f64 {
    if x.is_nan() { 0.0_f64.copysign(x) } else { x }
}

/// `copysign(isinf(x) ? 1 : 0, x)`, used to "box" infinities in C99 Annex G
fn box_inf(x: f64) -> f64 {
    (if x.is_infinite() { 1.0_f64 } else { 0.0 }).copysign(x)
}

/// `x - z` for real `x` as in C++, where the imaginary part is `-z.im` instead of `0.0 - z.im`
pub(crate) fn rsub(x: f64, z: Complex<f64>) -> Complex<f64> {
    Complex::new(x - z.re, -z.im)
}
