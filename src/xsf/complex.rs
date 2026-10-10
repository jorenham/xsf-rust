use num_complex::Complex;

/// Complex division `a / b` using Smith's algorithm, which avoids the overflow and underflow of
/// `|b|^2` in the naive formula. Tiny operands are first scaled by the same power of 2, so that
/// the intermediate results don't lose precision by becoming subnormal.
///
/// This replaces the complex division in C++, i.e. GCC's `__divdc3`. Unlike `__divdc3`, it
/// doesn't scale huge operands (e.g. `1 / (MAX + i MAX)` underflows to zero), and it doesn't
/// recover infinities from NaN results.
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
    if b.re.abs() >= b.im.abs() {
        let r = b.im / b.re;
        let den = b.re + b.im * r;
        Complex::new((a.re + a.im * r) / den, (a.im - a.re * r) / den)
    } else {
        let r = b.re / b.im;
        let den = b.re * r + b.im;
        Complex::new((a.re * r + a.im) / den, (a.im * r - a.re) / den)
    }
}
