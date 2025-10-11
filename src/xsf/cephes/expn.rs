use crate::ffi;

/// Generalized exponential integral Eₙ(x)
///
/// For integer n ≥ 0 and real x ≥ 0 the
/// generalized exponential integral is defined as:
///
/// Eₙ(x) = xⁿ⁻¹ ∫[x,∞) e⁻ᵗ / tⁿ dt
///
/// ## See also:
/// - [`exp1`](fn.exp1.html): Special case of Eₙ for n = 1
/// - [`expi`](fn.expi.html): Related to Eₙ when n = 1
pub fn expn(n: u32, x: f64) -> f64 {
    unsafe { ffi::xsf::expn(n as core::ffi::c_int, x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_expn_u32() {
        xsref::test::<f64, _>("expn", "p_d-d", |x: &[f64]| expn(x[0] as u32, x[1]));
    }

    #[test]
    fn test_expn_f64() {
        // not sure why this table exists; but we might as well use it 🤷🏻
        xsref::test::<f64, _>("expn", "d_d-d", |x: &[f64]| expn(x[0] as u32, x[1]));
    }
}
