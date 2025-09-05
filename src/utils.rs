use crate::bindings;

#[inline(always)]
pub(crate) fn c_complex64_nan() -> bindings::root::std::complex<f64> {
    num_complex::Complex::new(f64::NAN, f64::NAN).into()
}
