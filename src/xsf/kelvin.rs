//! Translated into pure Rust from `xsf/kelvin.h` (scipy/xsf@b1eb270).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the scalar functions (`ber`, ..., `keip`) share `kelvin` (xsf duplicates it).
//!
//! `pow(y, 2)` is written as `y * y`, which is what GCC compiles it to, and what LLVM only does
//! when optimizing.

use core::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4, PI, SQRT_2, TAU};

use num_complex::Complex64;

use crate::xsf::amos;
use crate::xsf::c64::{C64, c64};
use crate::xsf::specfun::{SPECFUN_INF, specfun_zconvinf};

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

/// The relative tolerance of the power series
const EPS: f64 = 1.0e-15;

/// The offsets `[i, j, k]` of the ber x, bei x, ber'x and bei'x series (see `ber_term`)
const BER_OFFSETS: [[f64; 3]; 4] = [
    [0.0, -1.0, -1.0],
    [0.0, 1.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, -1.0, 1.0],
];

/// The term `r_m = -r_(m-1) x4 / (4 m (m + i) (2m + j) (2m + k))` of a series with the offsets
/// `[i, j, k]` of `BER_OFFSETS`, where `x4 = (x/2)^4`
fn ber_term(r: f64, m: f64, x4: f64, [i, j, k]: [f64; 3]) -> f64 {
    // The products of these small integers are exact
    -0.25 * r / (m * (m + i)) / ((2.0 * m + j) * (2.0 * m + k)) * x4
}

/// ber x, bei x, ber'x and bei'x from their power series
fn ber_series(x: f64) -> (f64, f64, f64, f64) {
    let x2 = 0.25 * x * x;
    let x4 = x2 * x2;
    // `r_0 + r_1 + ...` with the `r_m` of `ber_term`
    let sum = |mut r: f64, d: [f64; 3]| {
        let mut s = r;
        for m in 1..=60 {
            r = ber_term(r, f64::from(m), x4, d);
            s += r;
            if r.abs() < s.abs() * EPS {
                break;
            }
        }
        s
    };
    (
        sum(1.0, BER_OFFSETS[0]),
        sum(x2, BER_OFFSETS[1]),
        sum(-0.25 * x * x2, BER_OFFSETS[2]),
        sum(0.5 * x, BER_OFFSETS[3]),
    )
}

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

    let (ber, bei, ger, gei, der, dei, her, hei);

    if x == 0.0 {
        return Klvna {
            ber: 1.0,
            bei: 0.0,
            ger: SPECFUN_INF,
            gei: -FRAC_PI_4,
            der: 0.0,
            dei: 0.0,
            her: -SPECFUN_INF,
            hei: 0.0,
        };
    }

    let x2 = 0.25 * x * x;
    let x4 = x2 * x2;

    if x.abs() < 3.0 {
        (ber, bei, der, dei) = ber_series(x);

        // `s + r_1 g_1 + r_2 g_2 + ...` with the `r_m` of `ber_term` and
        // `g_m = g_(m-1) + 1 / (2m + k - 1) + 1 / (2m + k)`
        let sum = |mut s: f64, mut r: f64, mut gs: f64, k: f64, d: [f64; 3]| {
            for m in 1..=60 {
                let m = f64::from(m);
                r = ber_term(r, m, x4, d);
                gs = gs + 1.0 / (2.0 * m + k - 1.0) + 1.0 / (2.0 * m + k);
                s += r * gs;
                if (r * gs).abs() < s.abs() * EPS {
                    break;
                }
            }
            s
        };
        let l = (x / 2.0).ln() + el;
        ger = sum(-l * ber + FRAC_PI_4 * bei, 1.0, 0.0, 0.0, BER_OFFSETS[0]);
        gei = sum(x2 - l * bei - FRAC_PI_4 * ber, x2, 1.0, 1.0, BER_OFFSETS[1]);
        let r = -0.25 * x * x2;
        let s = 1.5 * r - ber / x - l * der + FRAC_PI_4 * dei;
        her = sum(s, r, 1.5, 2.0, BER_OFFSETS[2]);
        let s = 0.5 * x - bei / x - l * dei - FRAC_PI_4 * der;
        hei = sum(s, 0.5 * x, 1.0, 1.0, BER_OFFSETS[3]);
    } else if x.abs() < 18.0 {
        // Unlike in xsf (scipy/xsf#232):
        // - ber, bei, ber' and bei' use the series up to x = 18, where the asymptotic expansion
        //   becomes the more accurate one;
        // - ker x + i kei x = K_0(z) and ker'x + i kei'x = -w K_1(z), with z = w x, w = e^(i pi/4)
        //   (DLMF 10.61.2, 10.29.3), since their series cancels from x = 3 on. AMOS can't
        //   overflow or underflow here.
        (ber, bei, der, dei) = ber_series(x);

        let w = c64(FRAC_1_SQRT_2, FRAC_1_SQRT_2);
        let mut k = [C64::ZERO; 2];
        amos::besk(w * x, 0.0, 1, 2, &mut k, &mut 0);
        let (ke, kep) = (k[0], -w * k[1]);
        (ger, gei, her, hei) = (ke.re, ke.im, kep.re, kep.im);
    } else {
        let km = if x.abs() >= 40.0 { 10 } else { 18 };
        // `(sum r_k c_k, sum (-1)^k r_k c_k, sum r_k s_k, sum (-1)^k r_k s_k)` over k = 0..=km,
        // with `(c_k, s_k) = (cos, sin)(k pi / 4)`, `r_0 = 1` and `r_k = next(r_(k-1), k, x)`
        let sums = |next: fn(f64, f64, f64) -> f64| {
            let (mut pp, mut pn, mut qp, mut qn) = (1.0, 1.0, 0.0, 0.0);
            let mut r = 1.0;
            let mut fac = 1.0;
            for k in 1..=km {
                let kf = f64::from(k);
                fac = -fac;
                let xt = 0.25 * kf * PI - (0.125 * kf).trunc() * TAU;
                let cs = xt.cos();
                let ss = xt.sin();
                r = next(r, kf, x);
                let rc = r * cs;
                let rs = r * ss;
                pp += rc;
                pn += fac * rc;
                qp += rs;
                qn += fac * rs;
            }
            (pp, pn, qp, qn)
        };
        let (pp0, pn0, qp0, qn0) = sums(|r0, k, x| 0.125 * r0 * (2.0 * k - 1.0).powi(2) / k / x);

        let xd = x / SQRT_2;
        let xe1 = xd.exp();
        let xe2 = (-xd).exp();
        let xc1 = 1.0 / (TAU * x).sqrt();
        let xc2 = (FRAC_PI_2 / x).sqrt();
        let cp0 = (xd + 0.125 * PI).cos();
        let cn0 = (xd - 0.125 * PI).cos();
        let sp0 = (xd + 0.125 * PI).sin();
        let sn0 = (xd - 0.125 * PI).sin();

        ger = xc2 * xe2 * (pn0 * cp0 - qn0 * sp0);
        gei = xc2 * xe2 * (-pn0 * sp0 - qn0 * cp0);
        ber = xc1 * xe1 * (pp0 * cn0 + qp0 * sn0) - gei / PI;
        bei = xc1 * xe1 * (pp0 * sn0 - qp0 * cn0) + ger / PI;

        let (pn1, pp1, qn1, qp1) =
            sums(|r1, k, x| 0.125 * r1 * (4.0 - (2.0 * k - 1.0).powi(2)) / (k * x));
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
pub fn ber(x: f64) -> f64 {
    kelvin(x)[0].re
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
pub fn bei(x: f64) -> f64 {
    kelvin(x)[0].im
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
    kelvin(x)[1].re
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
    kelvin(x)[1].im
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
pub fn berp(x: f64) -> f64 {
    kelvin(x)[2].re
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
pub fn beip(x: f64) -> f64 {
    kelvin(x)[2].im
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
    kelvin(x)[3].re
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
    kelvin(x)[3].im
}

/// Kelvin functions as complex numbers
///
/// For $x \ge 0$, these are defined as [^DLMF]:
///
/// $$
/// \ber x + i \bei x = J_0\left(x e^{3\pi i/4}\right) , \quad
/// \ker x + i \kei x = K_0\left(x e^{\pi i/4}\right) ,
/// $$
///
/// with $J_0$ the Bessel function of the first kind and $K_0$ the modified Bessel function of the
/// second kind.
///
/// Corresponds to [`scipy.special.kelvin`][scipy]
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kelvin.html
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/10.61.E1> and
///   <https://dlmf.nist.gov/10.61.E2>.
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
#[must_use]
#[inline]
pub fn kelvin(mut x: f64) -> [Complex64; 4] {
    let mut flag = false;
    if x < 0.0 {
        x = -x;
        flag = true;
    }

    let k = klvna(x);
    let mut be = Complex64::new(k.ber, k.bei);
    let mut ke = Complex64::new(k.ger, k.gei);
    let mut bep = Complex64::new(k.der, k.dei);
    let mut kep = Complex64::new(k.her, k.hei);

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
