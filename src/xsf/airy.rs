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
        // Unlike in xsf (scipy/xsf#310), fall back to Cephes if AMOS overflows (NaN Bi), e.g.
        // `airy(415.6)` was [0, 0, NaN, NaN] instead of [0, 0, inf, inf].
        if x < 0.0 || !zbi.re.is_nan() {
            return [zai.re, zaip.re, zbi.re, zbip.re];
        }
    }

    let (_, ai, aip, bi, bip) = crate::xsf::cephes::airy(x);
    [ai, aip, bi, bip]
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
    let eps = 1e-15;

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
/// For negative `x`, the integrals run from 0 down to `x`.
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

#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::excessive_precision,
    clippy::unreadable_literal
)]
fn airyb(x: f64, ai: &mut f64, bi: &mut f64, ad: &mut f64, bd: &mut f64) {
    // =======================================================
    // Purpose: Compute Airy functions and their derivatives
    // Input:   x  --- Argument of Airy function
    // Output:  AI --- Ai(x)
    //          BI --- Bi(x)
    //          AD --- Ai'(x)
    //          BD --- Bi'(x)
    // =======================================================

    const EPS: f64 = 1.0e-15;
    const C1: f64 = 0.355028053887817;
    const C2: f64 = 0.258819403792807;
    const SR3: f64 = 1.732050807568877;

    let mut km: i32;
    let mut km2: i32;
    let kmax: i32;
    let mut ck = [0.0_f64; 52];
    let mut dk = [0.0_f64; 52];
    let mut r: f64;

    km2 = 0;
    let xa = x.abs();
    let xq = xa.sqrt();
    let mut xm = 8.0;
    if x > 0.0 {
        xm = 5.0;
    }

    if x == 0.0 {
        *ai = C1;
        *bi = SR3 * C1;
        *ad = -C2;
        *bd = SR3 * C2;
        return;
    }

    if xa <= xm {
        let mut fx = 1.0;
        r = 1.0;
        for k in 1..=40 {
            r = r * x / (3.0 * f64::from(k)) * x / (3.0 * f64::from(k) - 1.0) * x;
            fx += r;
            if r.abs() < fx.abs() * EPS {
                break;
            }
        }

        let mut gx = x;
        r = x;
        for k in 1..=40 {
            r = r * x / (3.0 * f64::from(k)) * x / (3.0 * f64::from(k) + 1.0) * x;
            gx += r;
            if r.abs() < gx.abs() * EPS {
                break;
            }
        }

        *ai = C1 * fx - C2 * gx;
        *bi = SR3 * (C1 * fx + C2 * gx);

        let mut df = 0.5 * x * x;
        r = df;
        for k in 1..=40 {
            r = r * x / (3.0 * f64::from(k)) * x / (3.0 * f64::from(k) + 2.0) * x;
            df += r;
            if r.abs() < df.abs() * EPS {
                break;
            }
        }

        let mut dg = 1.0;
        r = 1.0;
        for k in 1..=40 {
            r = r * x / (3.0 * f64::from(k)) * x / (3.0 * f64::from(k) - 2.0) * x;
            dg += r;
            if r.abs() < dg.abs() * EPS {
                break;
            }
        }

        *ad = C1 * df - C2 * dg;
        *bd = SR3 * (C1 * df + C2 * dg);
    } else {
        // NOTE: `(int)(24.5 - xa)` is UB in C++ if `xa` is NaN or `xa > 24.5 + 2^31`:
        // - huge `xa`: `km` is overwritten below (`xa > 15`)
        // - NaN `x`: `km` is 0 here (INT_MIN in C++ on x86-64); all outputs are NaN either way
        km = (24.5 - xa) as i32;
        if xa < 6.0 {
            km = 14;
        }
        if xa > 15.0 {
            km = 10;
        }

        if x > 0.0 {
            kmax = km;
        } else {
            // Choose cutoffs so that the remainder term in asymptotic
            // expansion is epsilon size. The X<0 branch needs to be fast
            // in order to make AIRYZO efficient
            if xa > 70.0 {
                km = 3;
            }
            if xa > 500.0 {
                km = 2;
            }
            if xa > 1000.0 {
                km = 1;
            }

            km2 = km;
            if xa > 150.0 {
                km2 = 1;
            }
            if xa > 3000.0 {
                km2 = 0;
            }
            kmax = 2 * km + 1;
        }
        let xe = xa * xq / 1.5;
        let xr1 = 1.0 / xe;
        let xar = 1.0 / xq;
        let xf = xar.sqrt();
        let rp = 0.5641895835477563;
        r = 1.0;
        for k in 1..=kmax {
            r = r * (6.0 * f64::from(k) - 1.0) / 216.0 * (6.0 * f64::from(k) - 3.0) / f64::from(k)
                * (6.0 * f64::from(k) - 5.0)
                / (2.0 * f64::from(k) - 1.0);
            ck[(k - 1) as usize] = r;
            dk[(k - 1) as usize] = -(6.0 * f64::from(k) + 1.0) / (6.0 * f64::from(k) - 1.0) * r;
        }

        if x > 0.0 {
            let mut sai = 1.0;
            let mut sad = 1.0;
            r = 1.0;
            for k in 1..=km {
                r *= -xr1;
                sai += ck[(k - 1) as usize] * r;
                sad += dk[(k - 1) as usize] * r;
            }
            let mut sbi = 1.0;
            let mut sbd = 1.0;
            r = 1.0;
            for k in 1..=km {
                r *= xr1;
                sbi += ck[(k - 1) as usize] * r;
                sbd += dk[(k - 1) as usize] * r;
            }
            let xp1 = (-xe).exp();
            *ai = 0.5 * rp * xf * xp1 * sai;
            *bi = rp * xf / xp1 * sbi;
            *ad = -0.5 * rp / xf * xp1 * sad;
            *bd = rp / xf / xp1 * sbd;
        } else {
            let xcs = (xe + PI / 4.0).cos();
            let xss = (xe + PI / 4.0).sin();
            let mut ssa = 1.0;
            let mut sda = 1.0;
            r = 1.0;
            let xr2 = 1.0 / (xe * xe);
            for k in 1..=km {
                r *= -xr2;
                ssa += ck[(2 * k - 1) as usize] * r;
                sda += dk[(2 * k - 1) as usize] * r;
            }
            let mut ssb = ck[0] * xr1;
            let mut sdb = dk[0] * xr1;
            r = xr1;
            for k in 1..=km2 {
                r *= -xr2;
                ssb += ck[(2 * k) as usize] * r;
                sdb += dk[(2 * k) as usize] * r;
            }

            *ai = rp * xf * (xss * ssa - xcs * ssb);
            *bi = rp * xf * (xcs * ssa + xss * ssb);
            *ad = -rp / xf * (xcs * sda + xss * sdb);
            *bd = rp / xf * (xss * sda - xcs * sdb);
        }
    }
}

#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::cast_sign_loss,
    clippy::excessive_precision,
    clippy::unreadable_literal
)]
fn airyzo(nt: i32, kf: i32, xa: &mut [f64], xb: &mut [f64], xc: &mut [f64], xd: &mut [f64]) {
    // ========================================================
    // Purpose: Compute the first NT zeros of Airy functions
    //          Ai(x) and Ai'(x), a and a', and the associated
    //          values of Ai(a') and Ai'(a); and the first NT
    //          zeros of Airy functions Bi(x) and Bi'(x), b and
    //          b', and the associated values of Bi(b') and
    //          Bi'(b)
    // Input :  NT    --- Total number of zeros
    //          KF    --- Function code
    //                    KF=1 for Ai(x) and Ai'(x)
    //                    KF=2 for Bi(x) and Bi'(x)
    // Output:  XA(m) --- a, the m-th zero of Ai(x) or
    //                    b, the m-th zero of Bi(x)
    //          XB(m) --- a', the m-th zero of Ai'(x) or
    //                    b', the m-th zero of Bi'(x)
    //          XC(m) --- Ai(a') or Bi(b')
    //          XD(m) --- Ai'(a) or Bi'(b)
    //                    ( m --- Serial number of zeros )
    // Routine called: AIRYB for computing Airy functions and
    //                 their derivatives
    // =======================================================

    // NOTE: for `kf` other than 1 or 2 (not reachable from the public API), `rt` stays 0, so
    // `err` is NaN and the Newton loops never terminate (in C++ as well).
    let mut rt = 0.0;
    let mut rt0: f64;
    let mut u = 0.0;
    let mut u1 = 0.0;
    let mut x: f64;
    let (mut ai, mut bi, mut ad, mut bd) = (0.0, 0.0, 0.0, 0.0);
    let mut err: f64;

    for i in 1..=nt {
        rt0 = 0.0;
        if kf == 1 {
            u = 3.0 * PI * (4.0 * f64::from(i) - 1.0) / 8.0;
            u1 = 1.0 / (u * u);
        } else if kf == 2 {
            if i == 1 {
                rt0 = -1.17371;
            } else {
                u = 3.0 * PI * (4.0 * f64::from(i) - 3.0) / 8.0;
                u1 = 1.0 / (u * u);
            }
        }

        if rt0 == 0.0 {
            // DLMF 9.9.18
            rt0 = -(u * u).powf(1.0 / 3.0)
                * (1.0
                    + u1 * (5.0 / 48.0
                        + u1 * (-5.0 / 36.0
                            + u1 * (77125.0 / 82944.0 + u1 * (-108056875.0 / 6967296.0)))));
        }

        loop {
            x = rt0;
            airyb(x, &mut ai, &mut bi, &mut ad, &mut bd);

            if kf == 1 {
                rt = rt0 - ai / ad;
            } else if kf == 2 {
                rt = rt0 - bi / bd;
            }

            err = ((rt - rt0) / rt).abs();
            if err <= 1.0e-12 {
                break;
            }
            rt0 = rt;
        }

        xa[(i - 1) as usize] = rt;
        if err > 1.0e-14 {
            airyb(rt, &mut ai, &mut bi, &mut ad, &mut bd);
        }

        if kf == 1 {
            xd[(i - 1) as usize] = ad;
        } else if kf == 2 {
            xd[(i - 1) as usize] = bd;
        }
    }

    for i in 1..=nt {
        rt0 = 0.0;

        if kf == 1 {
            if i == 1 {
                rt0 = -1.01879;
            } else {
                u = 3.0 * PI * (4.0 * f64::from(i) - 3.0) / 8.0;
                u1 = 1.0 / (u * u);
            }
        } else if kf == 2 {
            if i == 1 {
                rt0 = -2.29444;
            } else {
                u = 3.0 * PI * (4.0 * f64::from(i) - 1.0) / 8.0;
                u1 = 1.0 / (u * u);
            }
        }

        if rt0 == 0.0 {
            // DLMF 9.9.19
            rt0 = -(u * u).powf(1.0 / 3.0)
                * (1.0
                    + u1 * (-7.0 / 48.0
                        + u1 * (35.0 / 288.0
                            + u1 * (-181223.0 / 207360.0 + u1 * (18683371.0 / 1244160.0)))));
        }

        loop {
            x = rt0;
            airyb(x, &mut ai, &mut bi, &mut ad, &mut bd);

            if kf == 1 {
                rt = rt0 - ad / (ai * x);
            } else if kf == 2 {
                rt = rt0 - bd / (bi * x);
            }

            err = ((rt - rt0) / rt).abs();
            if err <= 1.0e-12 {
                break;
            }
            rt0 = rt;
        }
        xb[(i - 1) as usize] = rt;

        if err > 1.0e-14 {
            airyb(rt, &mut ai, &mut bi, &mut ad, &mut bd);
        }

        if kf == 1 {
            xc[(i - 1) as usize] = ai;
        } else if kf == 2 {
            xc[(i - 1) as usize] = bi;
        }
    }
}

#[inline]
fn airy_zeros<const NT: usize, const KF: i32>() -> [[f64; NT]; 4] {
    let [mut x, mut xp, mut ai_xp, mut aip_x] = [
        [f64::NAN; NT],
        [f64::NAN; NT],
        [f64::NAN; NT],
        [f64::NAN; NT],
    ];
    if NT > 0 {
        airyzo(
            NT.try_into().unwrap(),
            KF,
            &mut x,
            &mut xp,
            &mut ai_xp,
            &mut aip_x,
        );
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
    airy_zeros::<NT, 1>()
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
    airy_zeros::<NT, 2>()
}
