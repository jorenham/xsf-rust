//! Translated into pure Rust from `xsf/fresnel.h` (xsf v0.2.2).
//!
//! Unlike in xsf (scipy/xsf#367), for `x < 0`, the imaginary part of F±(x) in `ffk` is computed from its own value
//! (as in the original Fortran `FI=(-1)**KS*PP2-FI`), instead of from the already updated real part.

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use core::ffi::c_int;
use num_complex::{Complex, c64};

use crate::utils::{SQRT_TAU, SQRT1OPI, SQRT2OPI};
use crate::xsf::complex::cdiv;
use crate::xsf::trig::{ccos, csin};

#[allow(clippy::many_single_char_names, clippy::float_cmp)]
fn cfc(z: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    // =========================================================
    // Purpose: Compute complex Fresnel integral C(z) and C'(z)
    // Input :  z --- Argument of C(z)
    // Output:  ZF --- C(z)
    //          ZD --- C'(z)
    // =========================================================

    const EPS: f64 = 1.0e-14;

    let w0 = z.norm();
    let zp = FRAC_PI_2 * z * z;
    let zp2 = zp * zp;
    let z0 = c64(0.0, 0.0);

    let c = if z == z0 {
        z0
    } else if w0 <= 2.5 {
        let mut cr = z;
        let mut c = cr;
        let mut wa0 = 0.0;
        for k in 1..=80 {
            let k = f64::from(k);
            cr = -0.5 * cr * (4.0 * k - 3.0) / k / (2.0 * k - 1.0) / (4.0 * k + 1.0) * zp2;
            c += cr;
            let wa = c.norm();
            if (((wa - wa0) / wa).abs() < EPS) && (k > 10.0) {
                return (c, ccos(FRAC_PI_2 * z * z));
            }
            wa0 = wa;
        }
        c
    } else if (w0 > 2.5) && (w0 < 4.5) {
        let m = 85;
        let mut c = z0;
        let mut cf = z0;
        let mut cf1 = z0;
        let mut cf0 = c64(1.0e-100, 0.0);
        for k in (0..=m).rev() {
            cf = cdiv((2.0 * f64::from(k) + 3.0) * cf0, zp) - cf1;
            if k % 2 == 0 {
                c += cf;
            }
            cf1 = cf0;
            cf0 = cf;
        }
        c *= cdiv(cdiv(2.0.into(), PI * z) * csin(zp), cf);
        c
    } else {
        // See comment at CFS(), use C(z) = iC(-iz)
        let d = if (z.im > -z.re) && (z.im <= z.re) {
            // right quadrant
            c64(0.5, 0.0)
        } else if (z.im > z.re) && (z.im >= -z.re) {
            // upper quadrant
            c64(0.0, 0.5)
        } else if (z.im < -z.re) && (z.im >= z.re) {
            // left quadrant
            c64(-0.5, 0.0)
        } else {
            c64(0.0, -0.5)
        };
        let mut cr = c64(1.0, 0.0);
        let mut cf = c64(1.0, 0.0);
        for k in 1..=20 {
            let k = f64::from(k);
            cr = cdiv(-0.25 * cr * (4.0 * k - 1.0) * (4.0 * k - 3.0), zp2);
            cf += cr;
        }
        cr = cdiv(1.0.into(), PI * z * z);
        let mut cg = cr;
        for k in 1..=12 {
            let k = f64::from(k);
            cr = cdiv(-0.25 * cr * (4.0 * k + 1.0) * (4.0 * k - 1.0), zp2);
            cg += cr;
        }
        d + cdiv(cf * csin(zp) - cg * ccos(zp), PI * z)
    };
    (c, ccos(FRAC_PI_2 * z * z))
}

#[allow(clippy::many_single_char_names, clippy::float_cmp)]
fn cfs(z: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    // =========================================================
    // Purpose: Compute complex Fresnel Integral S(z) and S'(z)
    // Input :  z  --- Argument of S(z)
    // Output:  ZF --- S(z)
    //          ZD --- S'(z)
    // =========================================================

    const EPS: f64 = 1.0e-14;

    let w0 = z.norm();
    let zp = FRAC_PI_2 * z * z;
    let zp2 = zp * zp;
    let z0 = c64(0.0, 0.0);

    let s = if z == z0 {
        z0
    } else if w0 <= 2.5 {
        let mut s = z * zp / 3.0;
        let mut cr = s;
        let mut wb0 = 0.0;
        for k in 1..=80 {
            let k = f64::from(k);
            cr = -0.5 * cr * (4.0 * k - 1.0) / k / (2.0 * k + 1.0) / (4.0 * k + 3.0) * zp2;
            s += cr;
            let wb = s.norm();
            if ((wb - wb0).abs() < EPS) && (k > 10.0) {
                return (s, csin(FRAC_PI_2 * z * z));
            }
            wb0 = wb;
        }
        s
    } else if (w0 > 2.5) && (w0 < 4.5) {
        let m = 85;
        let mut s = z0;
        let mut cf = z0;
        let mut cf1 = z0;
        let mut cf0 = c64(1.0e-100, 0.0);
        for k in (0..=m).rev() {
            cf = cdiv((2.0 * f64::from(k) + 3.0) * cf0, zp) - cf1;
            if k % 2 == 1 {
                s += cf;
            }
            cf1 = cf0;
            cf0 = cf;
        }
        cdiv(cdiv(2.0.into(), PI * z) * csin(zp), cf) * s
    } else {
        // Auxiliary functions f(z) and g(z) can be computed using an
        // asymptotic expansion in the right quadrant |arg(z)| <= pi/4, not pi/2
        // as sometimes suggested. Use the symmetry S(z) = -iS(-iz).
        // Interestingly, most of the expansion code is the same across
        // the quadrants. (The forth power in Z is the equalizer here.)
        // Only one constant has to be adapted.
        let d = if (z.im > -z.re) && (z.im <= z.re) {
            // right quadrant
            c64(0.5, 0.0)
        } else if (z.im > z.re) && (z.im >= -z.re) {
            // upper quadrant
            c64(0.0, -0.5)
        } else if (z.im < -z.re) && (z.im >= z.re) {
            // left quadrant
            c64(-0.5, 0.0)
        } else {
            c64(0.0, 0.5)
        };
        let mut cr = c64(1.0, 0.0);
        let mut cf = c64(1.0, 0.0);
        for k in 1..=20 {
            let k = f64::from(k);
            cr = cdiv(-0.25 * cr * (4.0 * k - 1.0) * (4.0 * k - 3.0), zp2);
            cf += cr;
        }
        cr = c64(1.0, 0.0);
        let mut cg = c64(1.0, 0.0);
        for k in 1..=12 {
            let k = f64::from(k);
            cr = cdiv(-0.25 * cr * (4.0 * k + 1.0) * (4.0 * k - 1.0), zp2);
            cg += cr;
        }
        cg = cdiv(cg, PI * z * z);
        d - cdiv(cf * ccos(zp) + cg * csin(zp), PI * z)
    };
    (s, csin(FRAC_PI_2 * z * z))
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names
)]
fn ffk(ks: i32, x: f64) -> (Complex<f64>, Complex<f64>) {
    // =======================================================
    // Purpose: Compute modified Fresnel integrals F±(x)
    //          and K±(x)
    // Input :  x   --- Argument of F±(x) and K±(x)
    //          KS  --- Sign code
    //                  KS=0 for calculating F+(x) and K+(x)
    //                  KS=1 for calculating F_(x) and K_(x)
    // Output:  FR  --- Re[F±(x)]
    //          FI  --- Im[F±(x)]
    //          GR  --- Re[K±(x)]
    //          GI  --- Im[K±(x)]
    // ======================================================

    const EPS: f64 = 1.0e-15;
    const PP2: f64 = 1.253_314_137_315_5;

    let mut f = c64(0.0, 0.0);
    let mut g = c64(0.0, 0.0);

    let xa = x.abs();
    let x2 = x * x;
    let x4 = x2 * x2;

    if x == 0.0 {
        f.re = 0.5 * FRAC_PI_2.sqrt();
        f.im = (-1.0_f64).powi(ks) * f.re;
        g = c64(0.5, 0.0);
    } else {
        let (c1, s1);
        if xa <= 2.5 {
            let mut xr = SQRT2OPI * xa;
            let mut c = xr;

            for k in 1..=50 {
                let k = f64::from(k);
                xr = -0.5 * xr * (4.0 * k - 3.0) / k / (2.0 * k - 1.0) / (4.0 * k + 1.0) * x4;
                c += xr;
                if (xr / c).abs() < EPS {
                    break;
                }
            }

            let mut s = SQRT2OPI * xa * xa * xa / 3.0;
            xr = s;

            for k in 1..=50 {
                let k = f64::from(k);
                xr = -0.5 * xr * (4.0 * k - 1.0) / k / (2.0 * k + 1.0) / (4.0 * k + 3.0) * x4;
                s += xr;
                if (xr / s).abs() < EPS {
                    break;
                }
            }

            (c1, s1) = (c, s);
            f.re = PP2 * (0.5 - c1);
            let fi0 = PP2 * (0.5 - s1);
            f.im = (-1.0_f64).powi(ks) * fi0;
        } else if xa < 5.5 {
            let m = (42.0 + 1.75 * x2) as i32;
            let mut xsu = 0.0;
            let mut xc = 0.0;
            let mut xs = 0.0;
            let mut xf1 = 0.0;
            let mut xf0 = 1.0e-100;

            for k in (0..=m).rev() {
                let xf = (2.0 * f64::from(k) + 3.0) * xf0 / x2 - xf1;
                if k % 2 == 0 {
                    xc += xf;
                } else {
                    xs += xf;
                }
                xsu += (2.0 * f64::from(k) + 1.0) * xf * xf;
                xf1 = xf0;
                xf0 = xf;
            }

            let xq = xsu.sqrt();
            let xw = SQRT2OPI * xa / xq;
            c1 = xc * xw;
            s1 = xs * xw;
        } else {
            let mut xr = 1.0;
            let mut xf = 1.0;

            for k in 1..=12 {
                let k = f64::from(k);
                xr = -0.25 * xr * (4.0 * k - 1.0) * (4.0 * k - 3.0) / x4;
                xf += xr;
            }

            xr = 1.0 / (2.0 * xa * xa);
            let mut xg = xr;

            for k in 1..=12 {
                let k = f64::from(k);
                xr = -0.25 * xr * (4.0 * k + 1.0) * (4.0 * k - 1.0) / x4;
                xg += xr;
            }

            c1 = 0.5 + (xf * x2.sin() - xg * x2.cos()) / SQRT_TAU / xa;
            s1 = 0.5 - (xf * x2.cos() + xg * x2.sin()) / SQRT_TAU / xa;
        }

        f.re = PP2 * (0.5 - c1);
        let fi0 = PP2 * (0.5 - s1);
        f.im = (-1.0_f64).powi(ks) * fi0;

        let xp = x2 + FRAC_PI_4;
        let cs = xp.cos();
        let ss = xp.sin();
        let xq2 = SQRT1OPI;

        g.re = xq2 * (f.re * cs + fi0 * ss);
        g.im = (-1.0_f64).powi(ks) * xq2 * (fi0 * cs - f.re * ss);

        if x < 0.0 {
            f.re = PP2 - f.re;
            f.im = (-1.0_f64).powi(ks) * PP2 - f.im;
            g.re = x2.cos() - g.re;
            g.im = -(-1.0_f64).powi(ks) * x2.sin() - g.im;
        }
    }
    (f, g)
}

#[allow(clippy::cast_precision_loss, clippy::many_single_char_names)]
fn fcszo(kf: i32, nt: usize, zo: &mut [Complex<f64>]) {
    // ===============================================================
    // Purpose: Compute the complex zeros of Fresnel integral C(z)
    //          or S(z) using modified Newton's iteration method
    // Input :  KF  --- Function code
    //                  KF=1 for C(z) or KF=2 for S(z)
    //          NT  --- Total number of zeros
    // Output:  ZO(L) --- L-th zero of C(z) or S(z)
    // Routines called:
    //      (1) CFC for computing Fresnel integral C(z)
    //      (2) CFS for computing Fresnel integral S(z)
    // ==============================================================

    let mut psq = 0.0;
    let mut w = 0.0;

    for nr in 1..=nt {
        if kf == 1 {
            psq = (4.0 * nr as f64 - 1.0).sqrt();
        }
        if kf == 2 {
            psq = 2.0 * (nr as f64).sqrt();
        }

        let px = psq - (PI * psq).ln() / (PI * PI * psq * psq * psq);
        let py = (PI * psq).ln() / (PI * psq);
        let mut z = c64(px, py);

        if kf == 2 {
            if nr == 2 {
                z = c64(2.8334, 0.2443);
            }
            if nr == 3 {
                z = c64(3.4674, 0.2185);
            }
            if nr == 4 {
                z = c64(4.0025, 0.2008);
            }
        }

        let mut it = 0;
        loop {
            it += 1;
            let (zf, zd) = if kf == 1 { cfc(z) } else { cfs(z) };

            let mut zp = c64(1.0, 0.0);
            for i in 1..nr {
                zp *= z - zo[i - 1];
            }

            let zfd = cdiv(zf, zp);
            let mut zq = c64(0.0, 0.0);
            for i in 1..nr {
                let mut zw = c64(1.0, 0.0);
                for j in 1..nr {
                    if j == i {
                        continue;
                    }
                    zw *= z - zo[j - 1];
                }
                zq += zw;
            }
            let zgd = cdiv(zd - zq * zfd, zp);
            z -= cdiv(zfd, zgd);
            let w0 = w;
            w = z.norm();
            if !((it <= 50) && (((w - w0) / w).abs() > 1.0e-12)) {
                break;
            }
        }
        zo[nr - 1] = z;
    }
}

pub trait FresnelArg: crate::sealed::Sealed + Sized {
    fn fresnel(self) -> (Self, Self);
}

impl FresnelArg for f64 {
    #[inline]
    fn fresnel(self) -> (Self, Self) {
        crate::xsf::cephes::fresnl(self)
    }
}

impl FresnelArg for Complex<f64> {
    #[inline]
    fn fresnel(self) -> (Self, Self) {
        (cfs(self).0, cfc(self).0)
    }
}

/// Fresnel integrals S(z) and C(z) for real or complex argument
///
/// The Fresnel integrals are defined as:
/// - S(z) = ∫₀ᶻ sin(π t² / 2) dt
/// - C(z) = ∫₀ᶻ cos(π t² / 2) dt
///
/// where the integrals are taken from 0 to z.
///
/// Corresponds to [`scipy.special.fresnel`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fresnel.html
///
/// # See also
/// - [`fresnel_zeros`] - Zeros of S(z) and C(z)
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals
#[must_use]
#[inline]
pub fn fresnel<T: FresnelArg>(z: T) -> (T, T) {
    z.fresnel()
}

/// Modified Fresnel positive integrals
///
/// Computes the modified Fresnel integrals:
/// - F₊(x) = ∫ₓ∞ exp(i t²) dt
/// - K₊(x) = F₊(x) × exp(-i (x² + π/4)) / √π
///
/// Corresponds to [`scipy.special.modfresnelp`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modfresnelp.html
///
/// # See also
/// - [`fresnel`] - Standard Fresnel integrals
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals
#[doc(alias = "modfresnelp")]
#[must_use]
#[inline]
pub fn modified_fresnel_plus(x: f64) -> (Complex<f64>, Complex<f64>) {
    ffk(0, x)
}

/// Modified Fresnel negative integrals
///
/// Computes the modified Fresnel integrals:
/// - F₋(x) = ∫ₓ∞ exp(-i t²) dt
/// - K₋(x) = F₋(x) × exp(i (x² + π/4)) / √π
///
/// Corresponds to [`scipy.special.modfresnelm`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.modfresnelm.html
///
/// # See also
/// - [`fresnel`] - Standard Fresnel integrals S(z) and C(z)
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals
#[doc(alias = "modfresnelm")]
#[must_use]
#[inline]
pub fn modified_fresnel_minus(x: f64) -> (Complex<f64>, Complex<f64>) {
    ffk(1, x)
}

/// Zeros of Fresnel integrals S(z) and C(z)
///
/// Compute `nt` complex zeros of the Fresnel integrals S(z) and C(z).
///
/// Corresponds to [`scipy.special.fresnel_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.fresnel_zeros.html
///
/// # Panics
/// - if `nt` is greater than `[c_int::MAX`](core::ffi::c_int::MAX)
///
/// # See also
/// - [`fresnel`] - Fresnel integrals S(z) and C(z)
/// - [`modified_fresnel_plus`] - Modified Fresnel positive integrals F₊(x) and K₊(x)
/// - [`modified_fresnel_minus`] - Modified Fresnel negative integrals F₋(x) and K₋(x)
#[doc(alias = "fresnelc_zeros", alias = "fresnels_zeros")]
#[must_use]
#[inline]
pub fn fresnel_zeros(nt: usize) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
    assert!(nt <= c_int::MAX as usize);

    let mut szo = vec![c64(f64::NAN, 0.0); nt];
    let mut czo = vec![c64(f64::NAN, 0.0); nt];
    fcszo(2, nt, &mut szo);
    fcszo(1, nt, &mut czo);
    (szo, czo)
}
