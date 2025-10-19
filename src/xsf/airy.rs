use num_complex::Complex;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for num_complex::Complex<f64> {}
}

pub trait AiryArg: sealed::Sealed + Sized {
    fn airy(self) -> [Self; 4];
    fn airye(self) -> [Self; 4];
}

const CNAN: Complex<f64> = Complex::new(f64::NAN, f64::NAN);

impl AiryArg for f64 {
    #[inline(always)]
    fn airy(self) -> [Self; 4] {
        let [mut ai, mut aip, mut bi, mut bip] = [f64::NAN; 4];
        unsafe {
            crate::ffi::xsf::airy(self, &mut ai, &mut aip, &mut bi, &mut bip);
        }
        [ai, aip, bi, bip]
    }

    #[inline(always)]
    fn airye(self) -> [Self; 4] {
        let [mut ai, mut aip, mut bi, mut bip] = [f64::NAN; 4];
        unsafe {
            crate::ffi::xsf::airye(self, &mut ai, &mut aip, &mut bi, &mut bip);
        }
        [ai, aip, bi, bip]
    }
}

impl AiryArg for num_complex::Complex<f64> {
    #[inline(always)]
    fn airy(self) -> [Self; 4] {
        let [mut ai, mut aip, mut bi, mut bip] = [CNAN; 4];
        unsafe {
            crate::ffi::xsf::airy_1(self, &mut ai, &mut bi, &mut aip, &mut bip);
        }
        [ai, aip, bi, bip]
    }

    #[inline(always)]
    fn airye(self) -> [Self; 4] {
        let [mut ai, mut bi, mut ad, mut bd] = [CNAN; 4];
        unsafe {
            crate::ffi::xsf::airye_1(self, &mut ai, &mut bi, &mut ad, &mut bd);
        }
        [ai, bi, ad, bd]
    }
}

/// Airy functions and derivatives
///
/// The Airy functions Ai(*z*) and Bi(*z*) are linearly independent solutions to the differential
/// equation *y(z)'' - z y(z) = 0*,  Airy equation or Stokes equation.
///
/// This corresponds to [`scipy.special.airy`][airy] in SciPy.
///
/// [airy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.airy.html
///
/// # Arguments
///
/// - `z: T`, where `T` is either
///   - [`f64`], or
///   - [`Complex<f64>`](num_complex::Complex).
///
/// # Returns
/// - Ai(*z*)
/// - Ai'(*z*)
/// - Bi(*z*)
/// - Bi'(*z*)
///
/// # See also
/// - [`airy_scaled`]: Exponentially scaled Airy functions and derivatives
pub fn airy<T: AiryArg>(z: T) -> [T; 4] {
    z.airy()
}

/// Exponentially scaled Airy functions and derivatives
///
/// This scales Ai and Ai' by a factor exp(2/3 z<sup>3/2</sup>), and Bi and Bi' by
/// exp(-2/3 |Re[z<sup>3/2</sup>]|).
///
/// This corresponds to [`scipy.special.airye`][airye] in SciPy.
///
/// [airye]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.airye.html
///
/// # Arguments
/// - `z`: Real ([`f64`]) or complex ([`Complex<f64>`](num_complex::Complex)) argument
///
/// # Returns
/// - exp(2/3 z<sup>3/2</sup>) Ai(z)
/// - exp(2/3 z<sup>3/2</sup>) Ai'(z)
/// - exp(-2/3 |Re[z<sup>3/2</sup>]|) Bi(z)
/// - exp(-2/3 |Re[z<sup>3/2</sup>]|) Bi'(z)
///
/// # See also
/// - [`airy`]
#[doc(alias = "airye")]
pub fn airy_scaled<T: AiryArg>(z: T) -> [T; 4] {
    z.airye()
}

/// Integrals of Airy functions
///
/// Calculates the integrals of Ai(*t*), Bi(*t*), Ai(-*t*), and Bi(-*t*) for *t from 0 to *x*.
///
/// If `x` is negative, the integrals are evaluated from 0 to `x` (i.e., from 0 down to `x`), following the same definitions for the Airy functions. The sign and order of the results reflect the direction of integration.
///
///
/// This corresponds to [`scipy.special.itairy`][itairy] in SciPy.
///
/// [itairy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.itairy.html
///
/// # Arguments
/// - `x` - Upper limit of the integral
///
/// # Returns
/// - ∫<sub>[0, x]</sub> Ai(*t*) *dt*
/// - ∫<sub>[0, x]</sub> Bi(*t*) *dt*
/// - ∫<sub>[0, x]</sub> Ai(*-t*) *dt*
/// - ∫<sub>[0, x]</sub> Bi(*-t*) *dt*
#[doc(alias = "itairy")]
pub fn airy_integrals(x: f64) -> [f64; 4] {
    let [mut itaip, mut itbip, mut itain, mut itbin] = [f64::NAN; 4];
    unsafe {
        crate::ffi::xsf::itairy(x, &mut itaip, &mut itbip, &mut itain, &mut itbin);
    }
    [itaip, itbip, itain, itbin]
}

#[inline(always)]
fn airyzo<const NT: usize, const KF: core::ffi::c_int>() -> [[f64; NT]; 4] {
    let [mut x, mut xp, mut ai_xp, mut aip_x] = [
        [f64::NAN; NT],
        [f64::NAN; NT],
        [f64::NAN; NT],
        [f64::NAN; NT],
    ];
    if NT > 0 {
        unsafe {
            crate::ffi::xsf::airyzo(
                NT,
                KF,
                x.as_mut_ptr(),
                xp.as_mut_ptr(),
                ai_xp.as_mut_ptr(),
                aip_x.as_mut_ptr(),
            );
        }
    }
    [x, xp, ai_xp, aip_x]
}

/// Zeros and values of the Airy function Ai and its derivative
///
/// This corresponds to [`scipy.special.ai_zeros`][ai_zeros] in SciPy.
///
/// [ai_zeros]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ai_zeros.html
///
/// # Returns
///
/// - Roots x<sub>k</sub> s.t. Ai(x<sub>k</sub>) = 0
/// - Roots x'<sub>k</sub> s.t. Ai'(x'<sub>k</sub>) = 0
/// - Critical points Ai(x'<sub>k</sub>), i.e. the values of Ai where its Ai' is zero
/// - Slopes Ai'(x<sub>k</sub>) at roots x<sub>k</sub>, i.e. the values of Ai' where Ai is zero
///
/// # See also
/// - [`airy_bi_zeros`]: Zeros of Bi and Bi'
#[doc(alias = "airyzo")]
#[doc(alias = "ai_zeros")]
pub fn airy_ai_zeros<const NT: usize>() -> [[f64; NT]; 4] {
    airyzo::<NT, 1>()
}

/// Zeros and values of the Airy function Bi and its derivative
///
/// This corresponds to [`scipy.special.bi_zeros`][bi_zeros] in SciPy.
///
/// [bi_zeros]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bi_zeros.html
///
/// # Returns
///
/// - Roots x<sub>k</sub> s.t. Bi(x<sub>k</sub>) = 0
/// - Roots x'<sub>k</sub> s.t. Bi'(x'<sub>k</sub>) = 0
/// - Critical points Bi(x'<sub>k</sub>), i.e. the values of Bi where its Bi' is zero
/// - Slopes Bi'(x<sub>k</sub>) at roots x<sub>k</sub>, i.e. the values of Bi' where Bi is zero
///
/// # See also
/// - [`airy_ai_zeros`]: Zeros of Ai and Ai'
#[doc(alias = "airyzo")]
#[doc(alias = "bi_zeros")]
pub fn airy_bi_zeros<const NT: usize>() -> [[f64; NT]; 4] {
    airyzo::<NT, 2>()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    #[inline(always)]
    fn as_tuple<T: Clone>(x: [T; 4]) -> (T, T, T, T) {
        (x[0].clone(), x[1].clone(), x[2].clone(), x[3].clone())
    }

    #[test]
    fn test_airy_f64() {
        crate::xsref::test("airy", "d-d_d_d_d", |x| as_tuple(crate::airy(x[0])));
    }

    #[test]
    fn test_airy_c64() {
        crate::xsref::test("airy", "cd-cd_cd_cd_cd", |x| {
            as_tuple(crate::airy(c64(x[0], x[1])))
        });
    }

    #[test]
    fn test_airye_f64() {
        crate::xsref::test("airye", "d-d_d_d_d", |x| as_tuple(crate::airy_scaled(x[0])));
    }

    #[test]
    fn test_airye_c64() {
        crate::xsref::test("airye", "cd-cd_cd_cd_cd", |x| {
            as_tuple(crate::airy_scaled(c64(x[0], x[1])))
        });
    }

    #[test]
    fn test_itairy() {
        crate::xsref::test("itairy", "d-d_d_d_d", |x| {
            as_tuple(crate::airy_integrals(x[0]))
        });
    }

    #[test]
    fn test_ai_zeros_0() {
        assert!(crate::airy_ai_zeros::<0>().iter().all(|a| a.is_empty()));
    }

    #[test]
    fn test_ai_zeros() {
        // >>> from scipy import special
        // >>> a, ap, ai, aip = special.ai_zeros(3)
        let [a, ap, ai, aip] = crate::airy_ai_zeros::<3>();
        // >>> a
        // array([-2.33810741, -4.08794944, -5.52055983])
        crate::np_assert_allclose!(&a, &[-2.33810741, -4.08794944, -5.52055983]);
        // >>> ap
        // array([-1.01879297, -3.24819758, -4.82009921])
        crate::np_assert_allclose!(&ap, &[-1.01879297, -3.24819758, -4.82009921]);
        // >>> ai
        // array([ 0.53565666, -0.41901548,  0.38040647])
        crate::np_assert_allclose!(&ai, &[0.53565666, -0.41901548, 0.38040647]);
        // >>> aip
        // array([ 0.70121082, -0.80311137,  0.86520403])
        crate::np_assert_allclose!(&aip, &[0.70121082, -0.80311137, 0.86520403]);
    }

    #[test]
    fn test_bi_zeros_0() {
        assert!(crate::airy_bi_zeros::<0>().iter().all(|a| a.is_empty()));
    }

    #[test]
    fn test_bi_zeros() {
        // >>> from scipy import special
        // >>> b, bp, bi, bip = special.bi_zeros(3)
        let [b, bp, bi, bip] = crate::airy_bi_zeros::<3>();
        // >>> b
        // array([-1.17371322, -3.2710933 , -4.83073784])
        crate::np_assert_allclose!(&b, &[-1.17371322, -3.2710933, -4.83073784]);
        // >>> bp
        // array([-2.29443968, -4.07315509, -5.51239573])
        crate::np_assert_allclose!(&bp, &[-2.29443968, -4.07315509, -5.51239573]);
        // >>> bi
        // array([-0.45494438,  0.39652284, -0.36796916])
        crate::np_assert_allclose!(&bi, &[-0.45494438, 0.39652284, -0.36796916]);
        // >>> bip
        // array([ 0.60195789, -0.76031014,  0.83699101])
        crate::np_assert_allclose!(&bip, &[0.60195789, -0.76031014, 0.83699101]);
    }
}
