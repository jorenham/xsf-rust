use alloc::{vec, vec::Vec};

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait AiryArg: sealed::Sealed + Sized {
    fn airy(self) -> (Self, Self, Self, Self);
    fn airye(self) -> (Self, Self, Self, Self);
}

impl AiryArg for f64 {
    #[inline(always)]
    fn airy(self) -> (Self, Self, Self, Self) {
        let mut ai = f64::NAN;
        let mut aip = f64::NAN;
        let mut bi = f64::NAN;
        let mut bip = f64::NAN;

        unsafe {
            crate::ffi::xsf::airy(self, &mut ai, &mut aip, &mut bi, &mut bip);
        }
        (ai, aip, bi, bip)
    }

    #[inline(always)]
    fn airye(self) -> (Self, Self, Self, Self) {
        let mut ai = f64::NAN;
        let mut aip = f64::NAN;
        let mut bi = f64::NAN;
        let mut bip = f64::NAN;

        unsafe {
            crate::ffi::xsf::airye(self, &mut ai, &mut aip, &mut bi, &mut bip);
        }
        (ai, aip, bi, bip)
    }
}

impl AiryArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn airy(self) -> (Self, Self, Self, Self) {
        let mut ai = f64::NAN.into();
        let mut bi = f64::NAN.into();
        let mut ad = f64::NAN.into();
        let mut bd = f64::NAN.into();

        unsafe {
            crate::ffi::xsf::airy_1(self.into(), &mut ai, &mut bi, &mut ad, &mut bd);
        }
        (ai.into(), bi.into(), ad.into(), bd.into())
    }

    #[inline(always)]
    fn airye(self) -> (Self, Self, Self, Self) {
        let mut ai = f64::NAN.into();
        let mut bi = f64::NAN.into();
        let mut ad = f64::NAN.into();
        let mut bd = f64::NAN.into();

        unsafe {
            crate::ffi::xsf::airye_1(self.into(), &mut ai, &mut bi, &mut ad, &mut bd);
        }
        (ai.into(), bi.into(), ad.into(), bd.into())
    }
}

/// Airy functions and their derivatives.
///
/// # Arguments
/// - `z` - real (`f64`) or complex (`num_complex::Complex<f64>`) argument
///
/// # Returns
/// - `Ai`: Ai(z)
/// - `Aip`: Ai'(z)
/// - `Bi`: Bi(z)
/// - `Bip`: Bi'(z)
pub fn airy<T: AiryArg>(z: T) -> (T, T, T, T) {
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
/// - `z` - real (`f64`) or complex (`num_complex::Complex<f64>`) argument
///
/// # Returns
/// - `eAi`: eAi(z)
/// - `eAip`: eAi'(z)
/// - `eBi`: eBi(z)
/// - `eBip`: eBi'(z)
pub fn airye<T: AiryArg>(z: T) -> (T, T, T, T) {
    z.airye()
}

/// Integrals of Airy functions
///
/// Calculates the integrals of Airy functions from 0 to `x`.
///
/// # Arguments
/// - `x` - Upper limit of the integral (x ≥ 0)
///
/// # Returns
/// - `Apt`: Integral of Ai(t) from 0 to x
/// - `Bpt`: Integral of Bi(t) from 0 to x
/// - `Ant`: Integral of Ai(-t) from 0 to x
/// - `Bnt`: Integral of Bi(-t) from 0 to x
pub fn itairy(x: f64) -> (f64, f64, f64, f64) {
    let mut apt = f64::NAN;
    let mut bpt = f64::NAN;
    let mut ant = f64::NAN;
    let mut bnt = f64::NAN;

    unsafe {
        crate::ffi::xsf::itairy(x, &mut apt, &mut bpt, &mut ant, &mut bnt);
    }
    (apt, bpt, ant, bnt)
}

/// Compute `nt` zeros and values of the Airy function Ai and its derivative
///
/// Computes the first `nt` zeros, `a`, of the Airy function Ai(x);
/// first `nt` zeros, `ap`, of the derivative of the Airy function Ai'(x);
/// the corresponding values Ai(a'); and the corresponding values Ai'(a).
///
/// Corresponds to `scipy.special.ai_zeros` in SciPy, and `xsf::airyzo` in the XSF C++ library.
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - `a`: First `nt` zeros of Ai(x)
/// - `ap`: First `nt` zeros of Ai'(x)
/// - `ai`: Values of Ai(x) evaluated at first `nt` zeros of Ai'(x)
/// - `aip`: Values of Ai'(x) evaluated at first `nt` zeros of Ai(x)
#[doc(alias = "airyzo")]
pub fn ai_zeros(nt: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    assert!(nt > 0);

    let mut a = vec![f64::NAN; nt];
    let mut ap = vec![f64::NAN; nt];
    let mut ai = vec![f64::NAN; nt];
    let mut aip = vec![f64::NAN; nt];

    unsafe {
        crate::ffi::xsf::airyzo(
            nt,
            1,
            a.as_mut_ptr(),
            ap.as_mut_ptr(),
            ai.as_mut_ptr(),
            aip.as_mut_ptr(),
        );
    }
    (a, ap, ai, aip)
}

/// Compute `nt` zeros and values of the Airy function Bi and its derivative
///
/// Computes the first `nt` zeros, b, of the Airy function Bi(x);
/// first `nt` zeros, b', of the derivative of the Airy function Bi'(x);
/// the corresponding values Bi(b'); and the corresponding values Bi'(b).
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - `b`: First `nt` zeros of Bi(x)
/// - `bp`: First `nt` zeros of Bi'(x)
/// - `bi`: Values of Bi(x) evaluated at first `nt` zeros of Bi'(x)
/// - `bip`: Values of Bi'(x) evaluated at first `nt` zeros of Bi(x)
pub fn bi_zeros(nt: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    assert!(nt > 0);

    let mut xa = vec![f64::NAN; nt];
    let mut xb = vec![f64::NAN; nt];
    let mut xc = vec![f64::NAN; nt];
    let mut xd = vec![f64::NAN; nt];

    unsafe {
        crate::ffi::xsf::airyzo(
            nt,
            2,
            xa.as_mut_ptr(),
            xb.as_mut_ptr(),
            xc.as_mut_ptr(),
            xd.as_mut_ptr(),
        );
    }
    (xa, xb, xc, xd)
}

#[cfg(test)]
mod tests {
    use crate::testing::np_assert_allclose;
    use num_complex::c64;

    #[test]
    fn test_airy_f64() {
        crate::xsref::test("airy", "d-d_d_d_d", |x| crate::airy(x[0]));
    }

    #[test]
    fn test_airy_c64() {
        crate::xsref::test("airy", "cd-cd_cd_cd_cd", |x| crate::airy(c64(x[0], x[1])));
    }

    #[test]
    fn test_airye_f64() {
        crate::xsref::test("airye", "d-d_d_d_d", |x| crate::airye(x[0]));
    }

    #[test]
    fn test_airye_c64() {
        crate::xsref::test("airye", "cd-cd_cd_cd_cd", |x| crate::airye(c64(x[0], x[1])));
    }

    #[test]
    fn test_itairy() {
        crate::xsref::test("itairy", "d-d_d_d_d", |x| crate::itairy(x[0]));
    }

    #[test]
    fn test_ai_zeros() {
        // >>> from scipy import special
        // >>> a, ap, ai, aip = special.ai_zeros(3)
        let (a, ap, ai, aip) = crate::ai_zeros(3);
        // >>> a
        // array([-2.33810741, -4.08794944, -5.52055983])
        np_assert_allclose(&a, &[-2.33810741, -4.08794944, -5.52055983], 0.0, 1e-8);
        // >>> ap
        // array([-1.01879297, -3.24819758, -4.82009921])
        np_assert_allclose(&ap, &[-1.01879297, -3.24819758, -4.82009921], 0.0, 1e-8);
        // >>> ai
        // array([ 0.53565666, -0.41901548,  0.38040647])
        np_assert_allclose(&ai, &[0.53565666, -0.41901548, 0.38040647], 0.0, 1e-8);
        // >>> aip
        // array([ 0.70121082, -0.80311137,  0.86520403])
        np_assert_allclose(&aip, &[0.70121082, -0.80311137, 0.86520403], 0.0, 1e-8);
    }

    #[test]
    fn test_bi_zeros() {
        // >>> from scipy import special
        // >>> b, bp, bi, bip = special.bi_zeros(3)
        let (b, bp, bi, bip) = crate::bi_zeros(3);
        // >>> b
        // array([-1.17371322, -3.2710933 , -4.83073784])
        np_assert_allclose(&b, &[-1.17371322, -3.2710933, -4.83073784], 0.0, 1e-8);
        // >>> bp
        // array([-2.29443968, -4.07315509, -5.51239573])
        np_assert_allclose(&bp, &[-2.29443968, -4.07315509, -5.51239573], 0.0, 1e-8);
        // >>> bi
        // array([-0.45494438,  0.39652284, -0.36796916])
        np_assert_allclose(&bi, &[-0.45494438, 0.39652284, -0.36796916], 0.0, 1e-8);
        // >>> bip
        // array([ 0.60195789, -0.76031014,  0.83699101])
        np_assert_allclose(&bip, &[0.60195789, -0.76031014, 0.83699101], 0.0, 1e-8);
    }
}
