/// Generalized exponential integral $E_n(x)$
///
/// Corresponds to [`scipy.special.expn`][expn]
///
/// [expn]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expn.html
///
/// # Definition
/// For integer $n \geq 0$ and real $x \geq 0$ the
/// generalized exponential integral is defined as:
///
/// $$ E_n(x) = x^{n-1} \int_x^\infty \frac{e^{-t}}{t^n} \dd t $$
///
/// # Panics
/// - If n cannot be converted to a [`c_int`](core::ffi::c_int)
///
/// # See also
/// - [`exp1`](fn.exp1.html): Special case of $E_n$ for $n = 1$
/// - [`expi`](fn.expi.html): Related to $E_n$ when $n = 1$
#[must_use]
#[inline]
pub fn expn(n: u32, x: f64) -> f64 {
    unsafe { crate::ffi::xsf::expn(n.try_into().unwrap(), x) }
}

#[cfg(test)]
mod tests {
    use num_traits::ToPrimitive;

    // The xsref tables also contain negative `n`, which in scipy returns NaN.

    #[test]
    fn test_expn_u32() {
        xsref::test("expn", "p_d-d", |x| {
            if x[0] < 0.0 {
                f64::NAN
            } else {
                crate::expn(x[0].to_u32().unwrap(), x[1])
            }
        });
    }

    #[test]
    fn test_expn_f64() {
        // not sure why this table exists; but we might as well use it 🤷🏻
        xsref::test("expn", "d_d-d", |x| {
            if x[0] < 0.0 {
                f64::NAN
            } else {
                crate::expn(x[0].to_u32().unwrap(), x[1])
            }
        });
    }
}
