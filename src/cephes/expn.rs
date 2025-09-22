use crate::bindings;

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
    unsafe { bindings::expn(n as core::ffi::c_int, x) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing;

    #[test]
    fn test_expn_u32() {
        testing::test::<f64, _>("expn", "p_d-d", |x: &[f64]| expn(x[0] as u32, x[1]));
    }

    #[test]
    fn test_expn_f64() {
        // not sure why this table exists; but we might as well use it 🤷🏻
        testing::test::<f64, _>("expn", "d_d-d", |x: &[f64]| expn(x[0] as u32, x[1]));
    }
}
