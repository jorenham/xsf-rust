use crate::bindings;
use num_complex::Complex;
use std::os::raw::c_int;

mod sealed {
    use num_complex::Complex;

    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for Complex<f64> {}
}

pub trait AiryArg: sealed::Sealed {
    type Output;

    fn airy(self) -> (Self::Output, Self::Output, Self::Output, Self::Output);
    fn airye(self) -> (Self::Output, Self::Output, Self::Output, Self::Output);
}

fn c_c64_nan() -> bindings::root::std::complex<f64> {
    Complex::new(f64::NAN, f64::NAN).into()
}

impl AiryArg for f64 {
    type Output = f64;

    fn airy(self) -> (Self::Output, Self::Output, Self::Output, Self::Output) {
        let mut ai = f64::NAN;
        let mut aip = f64::NAN;
        let mut bi = f64::NAN;
        let mut bip = f64::NAN;

        unsafe {
            bindings::airy(self, &mut ai, &mut aip, &mut bi, &mut bip);
        }
        (ai, aip, bi, bip)
    }

    fn airye(self) -> (Self::Output, Self::Output, Self::Output, Self::Output) {
        let mut ai = f64::NAN;
        let mut aip = f64::NAN;
        let mut bi = f64::NAN;
        let mut bip = f64::NAN;

        unsafe {
            bindings::airye(self, &mut ai, &mut aip, &mut bi, &mut bip);
        }
        (ai, aip, bi, bip)
    }
}

impl AiryArg for Complex<f64> {
    type Output = Complex<f64>;

    fn airy(self) -> (Self::Output, Self::Output, Self::Output, Self::Output) {
        let mut ai = c_c64_nan();
        let mut bi = c_c64_nan();
        let mut ad = c_c64_nan();
        let mut bd = c_c64_nan();

        unsafe {
            bindings::airy_1(self.into(), &mut ai, &mut bi, &mut ad, &mut bd);
        }
        (ai.into(), bi.into(), ad.into(), bd.into())
    }

    fn airye(self) -> (Self::Output, Self::Output, Self::Output, Self::Output) {
        let mut ai = c_c64_nan();
        let mut bi = c_c64_nan();
        let mut ad = c_c64_nan();
        let mut bd = c_c64_nan();

        unsafe {
            bindings::airye_1(self.into(), &mut ai, &mut bi, &mut ad, &mut bd);
        }
        (ai.into(), bi.into(), ad.into(), bd.into())
    }
}

/// Airy functions and their derivatives `Ai(z)`, `Ai'(z)`, `Bi(z)`, `Bi'(z)`
pub fn airy<T: AiryArg>(z: T) -> (T::Output, T::Output, T::Output, T::Output) {
    z.airy()
}

/// Exponentially scaled Airy functions and their derivatives `eAi(z)`, `eAi'(z)`, `eBi(z)`,
/// `eBi'(z)`
pub fn airye<T: AiryArg>(z: T) -> (T::Output, T::Output, T::Output, T::Output) {
    z.airye()
}

/// Compute the integrals of Airy functions with respect to t from 0 and x (x ≥ 0)
pub fn itairy(x: f64) -> (f64, f64, f64, f64) {
    let mut apt = f64::NAN;
    let mut bpt = f64::NAN;
    let mut ant = f64::NAN;
    let mut bnt = f64::NAN;

    unsafe {
        bindings::itairy(x, &mut apt, &mut bpt, &mut ant, &mut bnt);
    }
    (apt, bpt, ant, bnt)
}

/// Compute Airy functions and their derivatives
pub fn airyb(x: f64) -> (f64, f64, f64, f64) {
    let mut ai = f64::NAN;
    let mut bi = f64::NAN;
    let mut ad = f64::NAN;
    let mut bd = f64::NAN;

    unsafe {
        bindings::airyb(x, &mut ai, &mut bi, &mut ad, &mut bd);
    }
    (ai, bi, ad, bd)
}

/// Compute the first NT zeros of Airy functions Ai(x) and Ai'(x), a and a', and the associated
/// values of Ai(a') and Ai'(a); and the first NT zeros of Airy functions Bi(x) and Bi'(x), b and
/// b', and the associated values of Bi(b') and Bi'(b)
pub fn airyzo(nt: i32, kf: i32) -> (f64, f64, f64, f64) {
    let mut xa = f64::NAN;
    let mut xb = f64::NAN;
    let mut xc = f64::NAN;
    let mut xd = f64::NAN;

    unsafe {
        bindings::airyzo(nt as c_int, kf as c_int, &mut xa, &mut xb, &mut xc, &mut xd);
    }
    (xa, xb, xc, xd)
}
