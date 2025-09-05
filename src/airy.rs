use crate::bindings;
use crate::utils::c_complex64_nan;
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

impl AiryArg for f64 {
    type Output = f64;

    #[inline(always)]
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

    #[inline(always)]
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

    #[inline(always)]
    fn airy(self) -> (Self::Output, Self::Output, Self::Output, Self::Output) {
        let mut ai = c_complex64_nan();
        let mut bi = c_complex64_nan();
        let mut ad = c_complex64_nan();
        let mut bd = c_complex64_nan();

        unsafe {
            bindings::airy_1(self.into(), &mut ai, &mut bi, &mut ad, &mut bd);
        }
        (ai.into(), bi.into(), ad.into(), bd.into())
    }

    #[inline(always)]
    fn airye(self) -> (Self::Output, Self::Output, Self::Output, Self::Output) {
        let mut ai = c_complex64_nan();
        let mut bi = c_complex64_nan();
        let mut ad = c_complex64_nan();
        let mut bd = c_complex64_nan();

        unsafe {
            bindings::airye_1(self.into(), &mut ai, &mut bi, &mut ad, &mut bd);
        }
        (ai.into(), bi.into(), ad.into(), bd.into())
    }
}

/// Airy functions and their derivatives.
///
/// # Arguments
///
/// - `z` - real (`f64`) or complex (`num_complex::Complex<f64>`) argument
///
/// # Returns
///
/// A tuple `(Ai, Aip, Bi, Bip)` where:
/// - `Ai` - Ai(z)
/// - `Aip` - Ai'(z)
/// - `Bi` - Bi(z)
/// - `Bip` - Bi'(z)
pub fn airy<T: AiryArg>(z: T) -> (T::Output, T::Output, T::Output, T::Output) {
    z.airy()
}

/// Exponentially scaled Airy functions and their derivatives.
///
/// Scaling:
///
/// ```plain
/// eAi(z)  = Ai(z)  * exp(2/3 * z * sqrt(z))
/// eAi'(z) = Ai'(z) * exp(2/3 * z * sqrt(z))
/// eBi(z)  = Bi(z)  * exp(-abs(2/3 * (z * sqrt(z)).real))
/// eBi'(z) = Bi'(z) * exp(-abs(2/3 * (z * sqrt(z)).real))
/// ```
///
/// # Arguments
///
/// - `z` - real (`f64`) or complex (`num_complex::Complex<f64>`) argument
///
/// # Returns
///
/// A tuple `(eAi, eAip, eBi, eBip)` where:
/// - `eAi` - eAi(z)
/// - `eAip` - eAi'(z)
/// - `eBi` - eBi(z)
/// - `eBip` - eBi'(z)
pub fn airye<T: AiryArg>(z: T) -> (T::Output, T::Output, T::Output, T::Output) {
    z.airye()
}

/// Integrals of Airy functions
///
/// Calculates the integrals of Airy functions from 0 to `x`.
///
/// # Arguments
///
/// - `x` - Upper limit of the integral (x ≥ 0)
///
/// # Returns
///
/// A tuple `(Apt, Bpt, Ant, Bnt)` where:
/// - `Apt` - Integral of Ai(t) from 0 to x
/// - `Bpt` - Integral of Bi(t) from 0 to x
/// - `Ant` - Integral of Ai(-t) from 0 to x
/// - `Bnt` - Integral of Bi(-t) from 0 to x
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

/// Airy functions and their derivatives.
///
/// # Arguments
///
/// - `x` - Real argument
///
/// # Returns
///
/// A tuple `(Ai, Bi, Aip, Bip)` where:
/// - `Ai` - Ai(x)
/// - `Bi` - Bi(x)
/// - `Aip` - Ai'(x)
/// - `Bip` - Bi'(x)
pub fn airyb(x: f64) -> (f64, f64, f64, f64) {
    let mut ai = f64::NAN;
    let mut bi = f64::NAN;
    let mut aip = f64::NAN;
    let mut bip = f64::NAN;

    unsafe {
        bindings::airyb(x, &mut ai, &mut bi, &mut aip, &mut bip);
    }
    (ai, bi, aip, bip)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AiryKind {
    Ai = 1,
    Bi = 2,
}

/// Zeros of Airy functions and their associated values.
///
/// This function computes the first `nt` zeros of Airy functions Ai(x) and Ai'(x), a and a',
/// and the associated values of Ai(a') and Ai'(a); and the first `nt` zeros of Airy functions
/// Bi(x) and Bi'(x), b and b', and the associated values of Bi(b') and Bi'(b).
///
/// # Arguments
///
/// - `nt` - Total number of zeros to compute
/// - `kf` - Function code:
///   - `1` for Ai(x) and Ai'(x)
///   - `2` for Bi(x) and Bi'(x)
///
/// # Returns
///
/// A tuple `(xa, xb, xc, xd)` where:
/// - `xa` - The m-th zero *a* of Ai(x) or the m-th zero *b* of Bi(x)
/// - `xb` - The m-th zero *a'* of Ai'(x) or the m-th zero *b'* of Bi'(x)
/// - `xc` - Ai(*a'*) or Bi(*b'*)
/// - `xd` - Ai'(*a*) or Bi'(*b*)
///
/// where m is the serial number of zeros.
pub fn airyzo(nt: u32, kf: AiryKind) -> (f64, f64, f64, f64) {
    assert!(nt > 0);

    let mut xa = f64::NAN;
    let mut xb = f64::NAN;
    let mut xc = f64::NAN;
    let mut xd = f64::NAN;

    unsafe {
        bindings::airyzo(nt as c_int, kf as c_int, &mut xa, &mut xb, &mut xc, &mut xd);
    }
    (xa, xb, xc, xd)
}
