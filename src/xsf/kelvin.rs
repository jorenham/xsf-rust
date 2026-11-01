//! Translated into pure Rust from `xsf/kelvin.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM only does
//! when optimizing.

use core::f64::consts::{PI, SQRT_2};

use num_complex::{Complex64, c64};

/// `SPECFUN_ZCONVINF` (only the real part is checked)
#[allow(clippy::float_cmp)]
fn specfun_zconvinf(z: &mut Complex64) {
    if z.re == 1.0e300 {
        z.re = f64::INFINITY;
    }
    if z.re == -1.0e300 {
        z.re = f64::NEG_INFINITY;
    }
}

struct Klvna {
    ber: f64,
    bei: f64,
    ger: f64,
    gei: f64,
    der: f64,
    dei: f64,
    her: f64,
    hei: f64,
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::unreadable_literal
)]
fn klvna(x: f64) -> Klvna {
    // ======================================================
    // Purpose: Compute Kelvin functions ber x, bei x, ker x
    //          and kei x, and their derivatives  ( x > 0 )
    // Input :  x   --- Argument of Kelvin functions
    // Output:  BER --- ber x
    //          BEI --- bei x
    //          GER --- ker x
    //          GEI --- kei x
    //          DER --- ber'x
    //          DEI --- bei'x
    //          HER --- ker'x
    //          HEI --- kei'x
    // ================================================

    let el = 0.5772156649015329;
    let eps = 1.0e-15;

    let (mut ber, mut bei, mut ger, mut gei, mut der, mut dei, mut her, mut hei);
    let mut r;
    let mut gs;

    if x == 0.0 {
        return Klvna {
            ber: 1.0,
            bei: 0.0,
            ger: 1.0e+300,
            gei: -0.25 * PI,
            der: 0.0,
            dei: 0.0,
            her: -1.0e+300,
            hei: 0.0,
        };
    }

    let x2 = 0.25 * x * x;
    let x4 = x2 * x2;

    if x.abs() < 10.0 {
        ber = 1.0;
        r = 1.0;
        for m in 1..=60 {
            let t = 2.0 * f64::from(m) - 1.0;
            r = -0.25 * r / f64::from(m * m) / (t * t) * x4;
            ber += r;
            if r.abs() < ber.abs() * eps {
                break;
            }
        }

        bei = x2;
        r = x2;
        for m in 1..=60 {
            let t = 2.0 * f64::from(m) + 1.0;
            r = -0.25 * r / f64::from(m * m) / (t * t) * x4;
            bei += r;
            if r.abs() < bei.abs() * eps {
                break;
            }
        }

        ger = -((x / 2.0).ln() + el) * ber + 0.25 * PI * bei;

        r = 1.0;
        gs = 0.0;
        for m in 1..=60 {
            let t = 2.0 * f64::from(m) - 1.0;
            r = -0.25 * r / f64::from(m * m) / (t * t) * x4;
            gs = gs + 1.0 / (2.0 * f64::from(m) - 1.0) + 1.0 / (2.0 * f64::from(m));
            ger += r * gs;
            if (r * gs).abs() < ger.abs() * eps {
                break;
            }
        }

        gei = x2 - ((x / 2.0).ln() + el) * bei - 0.25 * PI * ber;

        r = x2;
        gs = 1.0;
        for m in 1..=60 {
            let t = 2.0 * f64::from(m) + 1.0;
            r = -0.25 * r / f64::from(m * m) / (t * t) * x4;
            gs = gs + 1.0 / (2.0 * f64::from(m)) + 1.0 / (2.0 * f64::from(m) + 1.0);
            gei += r * gs;
            if (r * gs).abs() < gei.abs() * eps {
                break;
            }
        }

        der = -0.25 * x * x2;
        r = der;
        for m in 1..=60 {
            let t = 2.0 * f64::from(m) + 1.0;
            r = -0.25 * r / f64::from(m) / (f64::from(m) + 1.0) / (t * t) * x4;
            der += r;
            if r.abs() < der.abs() * eps {
                break;
            }
        }

        dei = 0.5 * x;
        r = dei;
        for m in 1..=60 {
            r = -0.25 * r
                / f64::from(m * m)
                / (2.0 * f64::from(m) - 1.0)
                / (2.0 * f64::from(m) + 1.0)
                * x4;
            dei += r;
            if r.abs() < dei.abs() * eps {
                break;
            }
        }

        r = -0.25 * x * x2;
        gs = 1.5;
        her = 1.5 * r - ber / x - ((x / 2.0).ln() + el) * der + 0.25 * PI * dei;
        for m in 1..=60 {
            let t = 2.0 * f64::from(m) + 1.0;
            r = -0.25 * r / f64::from(m) / (f64::from(m) + 1.0) / (t * t) * x4;
            gs = gs + 1.0 / (2.0 * f64::from(m) + 1.0) + 1.0 / (f64::from(2 * m) + 2.0);
            her += r * gs;
            if (r * gs).abs() < her.abs() * eps {
                break;
            }
        }

        r = 0.5 * x;
        gs = 1.0;
        hei = 0.5 * x - bei / x - ((x / 2.0).ln() + el) * dei - 0.25 * PI * der;
        for m in 1..=60 {
            r = -0.25 * r / f64::from(m * m) / (f64::from(2 * m) - 1.0) / (f64::from(2 * m) + 1.0)
                * x4;
            gs = gs + 1.0 / (2.0 * f64::from(m)) + 1.0 / (f64::from(2 * m) + 1.0);
            hei += r * gs;
            if (r * gs).abs() < hei.abs() * eps {
                break;
            }
        }
    } else {
        let mut pp0 = 1.0;
        let mut pn0 = 1.0;
        let mut qp0 = 0.0;
        let mut qn0 = 0.0;
        let mut r0 = 1.0;
        let mut km = 18;
        if x.abs() >= 40.0 {
            km = 10;
        }
        let mut fac = 1.0;
        for k in 1..=km {
            let kf = f64::from(k);
            fac = -fac;
            let xt = 0.25 * kf * PI - (0.125 * kf).trunc() * 2.0 * PI;
            let cs = xt.cos();
            let ss = xt.sin();
            let t = 2.0 * kf - 1.0;
            r0 = 0.125 * r0 * (t * t) / kf / x;
            let rc = r0 * cs;
            let rs = r0 * ss;
            pp0 += rc;
            pn0 += fac * rc;
            qp0 += rs;
            qn0 += fac * rs;
        }

        let xd = x / SQRT_2;
        let xe1 = xd.exp();
        let xe2 = (-xd).exp();
        let xc1 = 1.0 / (2.0 * PI * x).sqrt();
        let xc2 = (0.5 * PI / x).sqrt();
        let cp0 = (xd + 0.125 * PI).cos();
        let cn0 = (xd - 0.125 * PI).cos();
        let sp0 = (xd + 0.125 * PI).sin();
        let sn0 = (xd - 0.125 * PI).sin();

        ger = xc2 * xe2 * (pn0 * cp0 - qn0 * sp0);
        gei = xc2 * xe2 * (-pn0 * sp0 - qn0 * cp0);
        ber = xc1 * xe1 * (pp0 * cn0 + qp0 * sn0) - gei / PI;
        bei = xc1 * xe1 * (pp0 * sn0 - qp0 * cn0) + ger / PI;

        let mut pp1 = 1.0;
        let mut pn1 = 1.0;
        let mut qp1 = 0.0;
        let mut qn1 = 0.0;
        let mut r1 = 1.0;
        fac = 1.0;
        for k in 1..=km {
            let kf = f64::from(k);
            fac = -fac;
            let xt = 0.25 * kf * PI - f64::from((0.125 * kf) as i32) * 2.0 * PI;
            let cs = xt.cos();
            let ss = xt.sin();
            let t = 2.0 * kf - 1.0;
            r1 = 0.125 * r1 * (4.0 - t * t) / (kf * x);
            let rc = r1 * cs;
            let rs = r1 * ss;
            pp1 += fac * rc;
            pn1 += rc;
            qp1 += fac * rs;
            qn1 += rs;
        }
        her = xc2 * xe2 * (-pn1 * cn0 + qn1 * sn0);
        hei = xc2 * xe2 * (pn1 * sn0 + qn1 * cn0);
        der = xc1 * xe1 * (pp1 * cp0 + qp1 * sp0) - hei / PI;
        dei = xc1 * xe1 * (pp1 * sp0 - qp1 * cp0) + her / PI;
    }
    Klvna {
        ber,
        bei,
        ger,
        gei,
        der,
        dei,
        her,
        hei,
    }
}

/// Kelvin function *ber*
///
/// Corresponds to [`scipy.special.ber`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ber.html
///
/// # See also
/// - [`bei`] - Kelvin function *bei*
/// - [`ker`] - Kelvin function *ker*
/// - [`kei`] - Kelvin function *kei*
#[must_use]
#[inline]
pub fn ber(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }

    let k = klvna(x);
    let mut be = c64(k.ber, k.bei);
    specfun_zconvinf(&mut be);
    be.re
}

/// Kelvin function *bei*
///
/// Corresponds to [`scipy.special.bei`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bei.html
///
/// # See also
/// - [`ber`] - Kelvin function *ber*
/// - [`ker`] - Kelvin function *ker*
/// - [`kei`] - Kelvin function *kei*
#[must_use]
#[inline]
pub fn bei(mut x: f64) -> f64 {
    if x < 0.0 {
        x = -x;
    }

    let k = klvna(x);
    let mut be = c64(k.ber, k.bei);
    specfun_zconvinf(&mut be);
    be.im
}

/// Kelvin function *ker*
///
/// Corresponds to [`scipy.special.ker`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ker.html
///
/// # See also
/// - [`ber`]
/// - [`bei`]
/// - [`kei`]
#[must_use]
#[inline]
pub fn ker(x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    let k = klvna(x);
    let mut ke = c64(k.ger, k.gei);
    specfun_zconvinf(&mut ke);
    ke.re
}

/// Kelvin function *kei*
///
/// Corresponds to [`scipy.special.kei`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kei.html
///
/// # See also
/// - [`ber`]
/// - [`bei`]
/// - [`ker`]
#[must_use]
#[inline]
pub fn kei(x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }

    let k = klvna(x);
    let mut ke = c64(k.ger, k.gei);
    specfun_zconvinf(&mut ke);
    ke.im
}

/// Derivative of the Kelvin function [`ber`]
///
/// Corresponds to [`scipy.special.berp`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.berp.html
///
/// # See also
/// - [`ber`]
/// - [`beip`]
/// - [`kelvin`]
#[must_use]
#[inline]
#[doc(alias = "ber_prime")]
pub fn berp(mut x: f64) -> f64 {
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }

    let k = klvna(x);
    let mut bep = c64(k.der, k.dei);
    specfun_zconvinf(&mut bep);
    if flag {
        return -bep.re;
    }
    bep.re
}

/// Derivative of the Kelvin function [`bei`]
///
/// Corresponds to [`scipy.special.beip`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.beip.html
///
/// # See also
/// - [`bei`]
/// - [`berp`]
/// - [`kelvin`]
#[must_use]
#[inline]
#[doc(alias = "bei_prime")]
pub fn beip(mut x: f64) -> f64 {
    let mut flag = false;

    if x < 0.0 {
        x = -x;
        flag = true;
    }
    let k = klvna(x);
    let mut bep = c64(k.der, k.dei);
    specfun_zconvinf(&mut bep);
    if flag {
        return -bep.im;
    }
    bep.im
}

/// Derivative of the Kelvin function [`ker`]
///
/// Corresponds to [`scipy.special.kerp`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kerp.html
///
/// # See also
/// - [`ker`]
/// - [`keip`]
/// - [`kelvin`]
#[must_use]
#[inline]
#[doc(alias = "ker_prime")]
pub fn kerp(x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    let k = klvna(x);
    let mut kep = c64(k.her, k.hei);
    specfun_zconvinf(&mut kep);
    kep.re
}

/// Derivative of the Kelvin function [`kei`]
///
/// Corresponds to [`scipy.special.keip`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.keip.html
///
/// # See also
/// - [`kei`]
/// - [`kerp`]
/// - [`kelvin`]
#[must_use]
#[inline]
#[doc(alias = "kei_prime")]
pub fn keip(x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    let k = klvna(x);
    let mut kep = c64(k.her, k.hei);
    specfun_zconvinf(&mut kep);
    kep.im
}

/// Kelvin functions as complex numbers
///
/// Corresponds to [`scipy.special.kelvin`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kelvin.html
///
/// # Returns
/// - *Be*: [`ber`] + *i* [`bei`]
/// - *Ke*: [`ker`] + *i* [`kei`]
/// - *Be*': [`berp`] + *i* [`beip`]
/// - *Ke*': [`kerp`] + *i* [`keip`]
///
/// # See also
/// - [`ber`]
/// - [`bei`]
/// - [`ker`]
/// - [`kei`]
/// - [`kelvin_zeros`]
///
/// Returned in a length-4 array of complex numbers:
/// - *Be*: Value of [`ber`] + *i* [`bei`]
/// - *Ke*: Value of [`ker`] + *i* [`kei`]
/// - *Be*': Derivative of [`berp`] + *i* [`beip`]
/// - *Ke*': Derivative of [`kerp`] + *i* [`keip`]
#[must_use]
#[inline]
pub fn kelvin(mut x: f64) -> [Complex64; 4] {
    let mut flag = false;
    if x < 0.0 {
        x = -x;
        flag = true;
    }

    let k = klvna(x);
    let mut be = c64(k.ber, k.bei);
    let mut ke = c64(k.ger, k.gei);
    let mut bep = c64(k.der, k.dei);
    let mut kep = c64(k.her, k.hei);

    specfun_zconvinf(&mut be);
    specfun_zconvinf(&mut ke);
    specfun_zconvinf(&mut bep);
    specfun_zconvinf(&mut kep);
    if flag {
        bep.re = -bep.re;
        bep.im = -bep.im;
        ke.re = f64::NAN;
        ke.im = f64::NAN;
        kep.re = f64::NAN;
        kep.im = f64::NAN;
    }
    [be, ke, bep, kep]
}

#[derive(Clone, Copy)]
enum KelvinFunction {
    Ber = 1,
    Bei = 2,
    Ker = 3,
    Kei = 4,
    Berp = 5,
    Beip = 6,
    Kerp = 7,
    Keip = 8,
}

#[allow(clippy::unreadable_literal)]
fn klvnzo(nt: usize, kd: KelvinFunction) -> Vec<f64> {
    // ====================================================
    // Purpose: Compute the zeros of Kelvin functions
    // Input :  NT  --- Total number of zeros
    //          KD  --- Function code
    //          KD=1 to 8 for ber x, bei x, ker x, kei x,
    //                    ber'x, bei'x, ker'x and kei'x,
    //                    respectively.
    // Output:  ZO(M) --- the M-th zero of Kelvin function
    //                    for code KD
    // Routine called:
    //          KLVNA for computing Kelvin functions and
    //          their derivatives
    // ====================================================

    let mut zo = vec![f64::NAN; nt];
    let mut rt0 = [
        0.0, 2.84891, 5.02622, 1.71854, 3.91467, 6.03871, 3.77268, 2.66584, 4.93181,
    ];
    let kd = kd as usize;
    let mut rt = rt0[kd];

    for m in 1..=nt {
        loop {
            let Klvna {
                ber,
                bei,
                ger,
                gei,
                der,
                dei,
                her,
                hei,
            } = klvna(rt);
            if kd == 1 {
                rt -= ber / der;
            } else if kd == 2 {
                rt -= bei / dei;
            } else if kd == 3 {
                rt -= ger / her;
            } else if kd == 4 {
                rt -= gei / hei;
            } else if kd == 5 {
                rt -= der / (-bei - der / rt);
            } else if kd == 6 {
                rt -= dei / (ber - dei / rt);
            } else if kd == 7 {
                rt -= her / (-gei - her / rt);
            } else {
                rt -= hei / (ger - hei / rt);
            }

            if (rt - rt0[kd]).abs() <= 5e-10 {
                break;
            }
            rt0[kd] = rt;
        }
        zo[m - 1] = rt;
        rt += 4.44;
    }
    zo
}

/// First `nt` zeros of Kelvin function [`ber`]
///
/// Corresponds to [`scipy.special.ber_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ber_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`ber`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`bei_zeros`]
#[must_use]
#[inline]
pub fn ber_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Ber)
}

/// First `nt` zeros of Kelvin function [`bei`]
///
/// Corresponds to [`scipy.special.bei_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.bei_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`bei`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`ber_zeros`]
#[must_use]
#[inline]
pub fn bei_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Bei)
}

/// First `nt` zeros of Kelvin function [`ker`]
///
/// Corresponds to [`scipy.special.ker_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ker_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`ker`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`kei_zeros`]
#[must_use]
#[inline]
pub fn ker_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Ker)
}

/// First `nt` zeros of Kelvin function [`kei`]
///
/// Corresponds to [`scipy.special.kei_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kei_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`kei`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`ker_zeros`]
#[must_use]
#[inline]
pub fn kei_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Kei)
}

/// First `nt` zeros of Kelvin function derivative [`berp`]
///
/// Corresponds to [`scipy.special.berp_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.berp_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`berp`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`beip_zeros`]
#[must_use]
#[inline]
pub fn berp_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Berp)
}

/// First `nt` zeros of Kelvin function derivative [`beip`]
///
/// Corresponds to [`scipy.special.beip_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.beip_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`beip`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`berp_zeros`]
#[must_use]
#[inline]
pub fn beip_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Beip)
}

/// First `nt` zeros of Kelvin function derivative [`kerp`]
///
/// Corresponds to [`scipy.special.kerp_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kerp_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`kerp`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`keip_zeros`]
#[must_use]
#[inline]
pub fn kerp_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Kerp)
}

/// First `nt` zeros of Kelvin function derivative [`keip`]
///
/// Corresponds to [`scipy.special.keip_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.keip_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute
///
/// # Returns
/// - Vector of length `nt` containing the smallest zeros of [`keip`]
///
/// # See also
/// - [`kelvin_zeros`]
/// - [`kerp_zeros`]
#[must_use]
#[inline]
pub fn keip_zeros(nt: usize) -> Vec<f64> {
    klvnzo(nt, KelvinFunction::Keip)
}

/// First `nt` zeros of all Kelvin functions and their derivatives
///
/// Corresponds to [`scipy.special.kelvin_zeros`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kelvin_zeros.html
///
/// # Arguments
/// - `nt` - Number of zeros to compute for each function
///
/// # Returns
/// - Length-8 array containing vectors of the zeros of [`ber`], [`bei`], [`ker`], [`kei`],
///   [`berp`], [`beip`], [`kerp`], and [`keip`], respectively
#[must_use]
#[inline]
pub fn kelvin_zeros(nt: usize) -> [Vec<f64>; 8] {
    [
        klvnzo(nt, KelvinFunction::Ber),
        klvnzo(nt, KelvinFunction::Bei),
        klvnzo(nt, KelvinFunction::Ker),
        klvnzo(nt, KelvinFunction::Kei),
        klvnzo(nt, KelvinFunction::Berp),
        klvnzo(nt, KelvinFunction::Beip),
        klvnzo(nt, KelvinFunction::Kerp),
        klvnzo(nt, KelvinFunction::Keip),
    ]
}
