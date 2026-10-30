use core::f64::consts::{PI, SQRT_2};

use num_complex::{Complex64, c64};

use crate::xsf::amos::{self, ierr_to_sferr, set_error_and_nan};

pub trait AiryArg: crate::sealed::Sealed + Sized {
    fn airy(self) -> [Self; 4];
    fn airye(self) -> [Self; 4];
}

impl AiryArg for f64 {
    #[inline]
    fn airy(self) -> [Self; 4] {
        airy_f64(self)
    }

    #[inline]
    fn airye(self) -> [Self; 4] {
        airye_f64(self)
    }
}

impl AiryArg for Complex64 {
    #[inline]
    fn airy(self) -> [Self; 4] {
        airy_c64(self)
    }

    #[inline]
    fn airye(self) -> [Self; 4] {
        airye_c64(self)
    }
}

// Translated from `xsf/airy.h` (xsf v0.2.2); returns `[ai, aip, bi, bip]`

fn airy_c64(z: Complex64) -> [Complex64; 4] {
    let mut id = 0;
    let mut ierr = 0;
    let kode = 1;
    let mut nz = 0;

    let mut ai = amos::airy(z, id, kode, &mut nz, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut ai);

    nz = 0;
    let mut bi = amos::biry(z, id, kode, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut bi);

    id = 1;
    let mut aip = amos::airy(z, id, kode, &mut nz, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut aip);

    nz = 0;
    let mut bip = amos::biry(z, id, kode, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut bip);

    [ai, aip, bi, bip]
}

fn airye_c64(z: Complex64) -> [Complex64; 4] {
    let mut id = 0;
    let kode = 2; /* Exponential scaling */
    let mut nz = 0;
    let mut ierr = 0;

    let mut ai = amos::airy(z, id, kode, &mut nz, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut ai);

    nz = 0;
    let mut bi = amos::biry(z, id, kode, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut bi);

    id = 1;
    let mut aip = amos::airy(z, id, kode, &mut nz, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut aip);

    nz = 0;
    let mut bip = amos::biry(z, id, kode, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut bip);

    [ai, aip, bi, bip]
}

fn airye_f64(z: f64) -> [f64; 4] {
    let mut id = 0;
    let kode = 2; /* Exponential scaling */
    let mut nz = 0;
    let mut ierr = 0;
    let (ai, aip, bi, bip);

    if z < 0.0 {
        ai = f64::NAN;
    } else {
        let mut cai = amos::airy(c64(z, 0.0), id, kode, &mut nz, &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cai);
        ai = cai.re;
    }

    nz = 0;
    let mut cbi = amos::biry(c64(z, 0.0), id, kode, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cbi);
    bi = cbi.re;

    id = 1;
    if z < 0.0 {
        aip = f64::NAN;
    } else {
        let mut caip = amos::airy(c64(z, 0.0), id, kode, &mut nz, &mut ierr);
        set_error_and_nan(ierr_to_sferr(nz, ierr), &mut caip);
        aip = caip.re;
    }

    nz = 0;
    let mut cbip = amos::biry(c64(z, 0.0), id, kode, &mut ierr);
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut cbip);
    bip = cbip.re;

    [ai, aip, bi, bip]
}

#[allow(clippy::manual_range_contains)]
fn airy_f64(x: f64) -> [f64; 4] {
    /* For small arguments, use Cephes as it's slightly faster.
     * For large arguments, use AMOS as it's more accurate.
     */
    if x < -10.0 || x > 10.0 {
        let [zai, zaip, zbi, zbip] = airy_c64(c64(x, 0.0));
        [zai.re, zaip.re, zbi.re, zbip.re]
    } else {
        let (_, ai, aip, bi, bip) = crate::xsf::cephes::airy(x);
        [ai, aip, bi, bip]
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
#[inline]
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
#[inline]
pub fn airy_scaled<T: AiryArg>(z: T) -> [T; 4] {
    z.airye()
}

#[allow(
    clippy::excessive_precision,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::unreadable_literal
)]
fn itairy_impl(mut x: f64) -> [f64; 4] {
    // ======================================================
    // Purpose: Compute the integrals of Airy fnctions with
    //          respect to t from 0 and x ( x ≥ 0 )
    // Input  : x   --- Upper limit of the integral
    // Output : APT --- Integration of Ai(t) from 0 and x
    //          BPT --- Integration of Bi(t) from 0 and x
    //          ANT --- Integration of Ai(-t) from 0 and x
    //          BNT --- Integration of Bi(-t) from 0 and x
    // ======================================================

    const A: [f64; 16] = [
        0.569444444444444,
        0.891300154320988,
        0.226624344493027e+01,
        0.798950124766861e+01,
        0.360688546785343e+02,
        0.198670292131169e+03,
        0.129223456582211e+04,
        0.969483869669600e+04,
        0.824184704952483e+05,
        0.783031092490225e+06,
        0.822210493622814e+07,
        0.945557399360556e+08,
        0.118195595640730e+10,
        0.159564653040121e+11,
        0.231369166433050e+12,
        0.358622522796969e+13,
    ];
    let c1 = 0.355028053887817;
    let c2 = 0.258819403792807;
    let sr3 = 1.732050807568877;
    let q0 = 1.0 / 3.0;
    let q1 = 2.0 / 3.0;
    let q2 = SQRT_2;
    let eps = 1e-5;

    let mut apt = 0.0;
    let mut bpt = 0.0;
    let mut ant = 0.0;
    let mut bnt = 0.0;

    if x == 0.0 {
        return [apt, bpt, ant, bnt];
    }
    if x.abs() <= 9.25 {
        for l in 0..2 {
            x *= (-1.0_f64).powf(f64::from(l));
            let mut fx = x;
            let mut r = x;
            for k in 1..41 {
                let k = f64::from(k);
                r = r * (3.0 * k - 2.0) / (3.0 * k + 1.0) * x / (3.0 * k) * x / (3.0 * k - 1.0) * x;
                fx += r;
                if r.abs() < fx.abs() * eps {
                    break;
                }
            }
            let mut gx = 0.5 * x * x;
            r = gx;
            for k in 1..41 {
                let k = f64::from(k);
                r = r * (3.0 * k - 1.0) / (3.0 * k + 2.0) * x / (3.0 * k) * x / (3.0 * k + 1.0) * x;
                gx += r;
                if r.abs() < gx.abs() * eps {
                    break;
                }
            }
            ant = c1 * fx - c2 * gx;
            bnt = sr3 * (c1 * fx + c2 * gx);
            if l == 0 {
                apt = ant;
                bpt = bnt;
            } else {
                ant = -ant;
                bnt = -bnt;
                x = -x;
            }
        }
    } else {
        let xe = x * x.sqrt() / 1.5;
        let xp6 = 1.0 / (6.0 * PI * xe).sqrt();
        let mut su1 = 1.0;
        let mut r = 1.0;
        let xr1 = 1.0 / xe;
        for k in 1..17 {
            r = -r * xr1;
            su1 += A[k - 1] * r;
        }
        let mut su2 = 1.0;
        r = 1.0;
        for k in 1..17 {
            r *= xr1;
            su2 += A[k - 1] * r;
        }
        apt = q0 - (-xe).exp() * xp6 * su1;
        bpt = 2.0 * xe.exp() * xp6 * su2;
        let mut su3 = 1.0;
        r = 1.0;
        let xr2 = 1.0 / (xe * xe);
        for k in 1..9 {
            r = -r * xr2;
            su3 += A[2 * k - 1] * r;
        }
        let mut su4 = A[0] * xr1;
        r = xr1;
        for k in 1..8 {
            r = -r * xr2;
            su4 += A[2 * k] * r;
        }
        let su5 = su3 + su4;
        let su6 = su3 - su4;
        ant = q1 - q2 * xp6 * (su5 * xe.cos() - su6 * xe.sin());
        bnt = q2 * xp6 * (su5 * xe.sin() + su6 * xe.cos());
    }
    [apt, bpt, ant, bnt]
}

/// Integrals of Airy functions
///
/// Calculates the integrals of Ai(*t*), Bi(*t*), Ai(-*t*), and Bi(-*t*) for *t* from 0 to *x*.
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
#[must_use]
#[inline]
pub fn airy_integrals(mut x: f64) -> [f64; 4] {
    let x_signbit = x.is_sign_negative();
    if x_signbit {
        x = -x;
    }

    let [mut apt, mut bpt, mut ant, mut bnt] = itairy_impl(x);
    if x_signbit {
        /* negative limit -- switch signs and roles */
        let mut tmp = apt;
        apt = -ant;
        ant = -tmp;

        tmp = bpt;
        bpt = -bnt;
        bnt = -tmp;
    }
    [apt, bpt, ant, bnt]
}

#[inline]
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
#[doc(alias = "airyzo", alias = "ai_zeros")]
#[must_use]
#[inline]
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
#[doc(alias = "airyzo", alias = "bi_zeros")]
#[must_use]
#[inline]
pub fn airy_bi_zeros<const NT: usize>() -> [[f64; NT]; 4] {
    airyzo::<NT, 2>()
}

#[cfg(test)]
mod tests {
    use num_complex::c64;

    fn as_tuple<T: Clone>(x: &[T; 4]) -> (T, T, T, T) {
        (x[0].clone(), x[1].clone(), x[2].clone(), x[3].clone())
    }

    #[test]
    fn test_airy_f64() {
        xsref::test("airy", "d-d_d_d_d", |x| as_tuple(&crate::airy(x[0])));
    }

    #[test]
    fn test_airy_c64() {
        xsref::test("airy", "cd-cd_cd_cd_cd", |x| {
            as_tuple(&crate::airy(c64(x[0], x[1])))
        });
    }

    #[test]
    fn test_airye_f64() {
        xsref::test("airye", "d-d_d_d_d", |x| {
            as_tuple(&crate::airy_scaled(x[0]))
        });
    }

    #[test]
    fn test_airye_c64() {
        xsref::test("airye", "cd-cd_cd_cd_cd", |x| {
            as_tuple(&crate::airy_scaled(c64(x[0], x[1])))
        });
    }

    #[test]
    fn test_itairy() {
        xsref::test("itairy", "d-d_d_d_d", |x| {
            as_tuple(&crate::airy_integrals(x[0]))
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
        crate::np_assert_allclose!(&a, &[-2.338_107_41, -4.087_949_44, -5.520_559_83]);
        // >>> ap
        // array([-1.01879297, -3.24819758, -4.82009921])
        crate::np_assert_allclose!(&ap, &[-1.018_792_97, -3.248_197_58, -4.820_099_21]);
        // >>> ai
        // array([ 0.53565666, -0.41901548,  0.38040647])
        crate::np_assert_allclose!(&ai, &[0.535_656_66, -0.419_015_48, 0.380_406_47]);
        // >>> aip
        // array([ 0.70121082, -0.80311137,  0.86520403])
        crate::np_assert_allclose!(&aip, &[0.701_210_82, -0.803_111_37, 0.865_204_03]);
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
        crate::np_assert_allclose!(&b, &[-1.173_713_22, -3.271_093_3, -4.830_737_84]);
        // >>> bp
        // array([-2.29443968, -4.07315509, -5.51239573])
        crate::np_assert_allclose!(&bp, &[-2.294_439_68, -4.073_155_09, -5.512_395_73]);
        // >>> bi
        // array([-0.45494438,  0.39652284, -0.36796916])
        crate::np_assert_allclose!(&bi, &[-0.454_944_38, 0.396_522_84, -0.367_969_16]);
        // >>> bip
        // array([ 0.60195789, -0.76031014,  0.83699101])
        crate::np_assert_allclose!(&bip, &[0.601_957_89, -0.760_310_14, 0.836_991_01]);
    }
}
