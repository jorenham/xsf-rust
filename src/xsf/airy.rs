use core::f64::consts::{PI, SQRT_2};

use num_complex::Complex64;

use crate::utils::mul_exp;
use crate::xsf::amos::{self, ierr_to_sferr, set_error_and_nan};
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes;

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
        airy_c64(self.into(), 1).map(Into::into)
    }

    #[inline]
    fn airye(self) -> [Self; 4] {
        airy_c64(self.into(), 2).map(Into::into)
    }
}

// Translated from `xsf/airy.h` (scipy/xsf@76d55c4); returns `[ai, aip, bi, bip]`

/// Unlike in xsf, `airy` and `airye` share this implementation (xsf duplicates it), with
/// `kode = 2` for exponential scaling.
fn airy_c64(z: C64, kode: i32) -> [C64; 4] {
    let ai = amos_airy_or_biry(false, z, 0, kode);
    let bi = amos_airy_or_biry(true, z, 0, kode);
    let aip = amos_airy_or_biry(false, z, 1, kode);
    let bip = amos_airy_or_biry(true, z, 1, kode);
    [ai, aip, bi, bip]
}

/// `amos::airy` (Ai) or, for `bi`, `amos::biry` (Bi), or their derivative for `id = 1`, with its
/// error reported and NaN where it fails
fn amos_airy_or_biry(bi: bool, z: C64, id: i32, kode: i32) -> C64 {
    let (mut nz, mut ierr) = (0, 0);
    let mut out = if bi {
        amos::biry(z, id, kode, &mut ierr)
    } else {
        amos::airy(z, id, kode, &mut nz, &mut ierr)
    };
    set_error_and_nan(ierr_to_sferr(nz, ierr), &mut out);
    out
}

fn airye_f64(z: f64) -> [f64; 4] {
    // Unlike in xsf, Cephes is used for z >= 0: it is faster than AMOS, whose series cancel for
    // Ai and Ai' at z < 2.09, and it has no upper limit.
    if z >= 0.0 {
        return cephes::airye(z);
    }
    // exponentially scaled
    let bi = |id| amos_airy_or_biry(true, c64(z, 0.0), id, 2).re;
    [f64::NAN, f64::NAN, bi(0), bi(1)]
}

fn airy_f64(x: f64) -> [f64; 4] {
    let cephes = || {
        let (_, ai, aip, bi, bip) = crate::xsf::cephes::airy(x);
        [ai, aip, bi, bip]
    };
    /* For small arguments, use Cephes as it's slightly faster.
     * For large arguments, use AMOS as it's more accurate.
     */
    // Unlike in xsf (scipy/xsf#532), Cephes is also used for x < -10, where its double-double
    // phase makes it more accurate than AMOS.
    if x > 10.0 {
        let [ai, aip, bi, bip] = airy_c64(c64(x, 0.0), 1).map(|z| z.re);
        // Unlike in xsf (scipy/xsf#310), Cephes or `airye` takes over where AMOS over- or
        // underflows.
        if x > 103.892 && x < 108.0 {
            // Cephes returns (0, 0, inf, inf) for x > 103.892, while Bi(x) only overflows at
            // x = 104.44 and Ai(x) underflows at x = 107.47, so `airye` is unscaled in between.
            let [ai, aip, bi, bip] = airye_f64(x);
            let zeta = 2.0 * x * x.sqrt() / 3.0;
            [
                mul_exp(ai, -zeta),
                mul_exp(aip, -zeta),
                mul_exp(bi, zeta),
                mul_exp(bip, zeta),
            ]
        } else if bip.is_nan() {
            // Bi' overflows first
            cephes()
        } else if ai == 0.0 {
            // AMOS returns 0 for Ai(x) below about 1e-305
            let [ai, aip, ..] = cephes();
            [ai, aip, bi, bip]
        } else {
            [ai, aip, bi, bip]
        }
    } else {
        cephes()
    }
}

/// Airy functions and derivatives
///
/// The Airy functions $\mathrm{Ai}(z)$ and $\mathrm{Bi}(z)$ are linearly independent solutions to
/// the differential equation $y\'\'(z) - z y(z) = 0$,  Airy equation or Stokes equation.
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
/// - $\mathrm{Ai}(z)$
/// - $\mathrm{Ai}\'(z)$
/// - $\mathrm{Bi}(z)$
/// - $\mathrm{Bi}\'(z)$
///
/// # See also
/// - [`airy_scaled`]: Exponentially scaled Airy functions and derivatives
#[inline]
pub fn airy<T: AiryArg>(z: T) -> [T; 4] {
    z.airy()
}

/// Exponentially scaled Airy functions and derivatives
///
/// This scales $\mathrm{Ai}$ and $\mathrm{Ai}\'$ by a factor $\exp(\tfrac{2}{3} z^{3/2})$, and
/// $\mathrm{Bi}$ and $\mathrm{Bi}\'$ by $\exp(-\tfrac{2}{3} \abs{\mathrm{Re}(z^{3/2})})$.
///
/// This corresponds to [`scipy.special.airye`][airye] in SciPy.
///
/// [airye]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.airye.html
///
/// # Arguments
/// - `z`: Real ([`f64`]) or complex ([`Complex<f64>`](num_complex::Complex)) argument
///
/// # Returns
/// - $\exp(\tfrac{2}{3} z^{3/2}) \mathrm{Ai}(z)$
/// - $\exp(\tfrac{2}{3} z^{3/2}) \mathrm{Ai}\'(z)$
/// - $\exp(-\tfrac{2}{3} \abs{\mathrm{Re}(z^{3/2})}) \mathrm{Bi}(z)$
/// - $\exp(-\tfrac{2}{3} \abs{\mathrm{Re}(z^{3/2})}) \mathrm{Bi}\'(z)$
///
/// # See also
/// - [`airy`]
#[doc(alias = "airye")]
#[inline]
pub fn airy_scaled<T: AiryArg>(z: T) -> [T; 4] {
    z.airye()
}

/// `int_z^inf Ai(t) dt` for `|ph z| <= pi/3`, given `zeta = 2/3 z^(3/2)`
///
/// By DLMF 9.6.2 and 10.32.9, this is
/// `1/(pi sqrt(3)) int_0^inf e^(-zeta cosh u) cosh(u/3) / cosh(u) du`. On the path
/// `cosh u = 1 + v^2 / zeta`, i.e. `u = 2 asinh(v / sqrt(2 zeta))`, the exponent is `-zeta - v^2`,
/// so that the trapezoidal rule in `v` converges exponentially:
///
/// - the step `H` keeps the error from the poles at `v = +-i sqrt(zeta)`, i.e.
///   `exp(-pi sqrt(2 |zeta|) / H)` for `ph zeta = pi/2`, below `2^-53` for `|zeta| > 1.88`
///   (`x > 2`)
/// - the last node `v = 41 H` has `e^(-v^2) < 2^-54`
fn ai_tail(zeta: C64) -> C64 {
    const H: f64 = 0.15;
    let r = 1.0 / (2.0 * zeta).sqrt();
    let sum = (1..=41)
        .map(|k| {
            let v = H * f64::from(k);
            let s = v * r; // sinh(u/2)
            let c = (1.0 + s * s).sqrt(); // cosh(u/2)
            let u3 = (s + c).ln() * (2.0 / 3.0); // u/3
            (-v * v).exp() * u3.cosh() / ((1.0 + 2.0 * s * s) * c)
        })
        .fold(C64::ZERO, |sum, t| sum + t);
    // (+ 0.5 for the node v = 0, whose integrand is 1, with the trapezoidal weight 1/2)
    (sum + 0.5) * (-zeta).exp() * (2.0 * H / PI) / (6.0 * zeta).sqrt()
}

// Translated from `xsf/cpu/airy.h` (scipy/xsf@306433e)
fn itairy_impl(x: f64) -> [f64; 4] {
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
    // the power series of int_0^x of f and g (DLMF 9.4.1-9.4.3), whose term ratios are
    // (3k + o0) / (3k + o1) x^3 / (3k (3k + o2)) with the offsets `o`
    let series = |x: f64| {
        let sum = |mut r: f64, o: [f64; 3]| {
            let mut sum = r;
            // Unlike in xsf (scipy/xsf#448), allow 80 terms instead of 40, which cover x <= 20
            for k in 1..81 {
                let k = 3.0 * f64::from(k);
                r = r * (k + o[0]) / (k + o[1]) * x / k * x / (k + o[2]) * x;
                sum += r;
                if r.abs() < sum.abs() * eps {
                    break;
                }
            }
            sum
        };
        let (fx, gx) = (
            sum(x, [-2.0, 1.0, -1.0]),
            sum(0.5 * x * x, [-1.0, 2.0, 1.0]),
        );
        (c1 * fx - c2 * gx, sr3 * (c1 * fx + c2 * gx))
    };

    let xe = x * x.sqrt() / 1.5;
    if x <= 2.0 {
        (apt, bpt) = series(x);
        (ant, bnt) = series(-x);
        ant = -ant;
        bnt = -bnt;
    } else if x <= 20.0 {
        // Unlike in xsf (scipy/xsf#448), use the series of BPT (with positive terms) and the tail
        // integral of Ai for 2 < x <= 20:
        // - APT = 1/3 - int_x^inf Ai(t) dt
        // - ANT + i BNT = 2/3 - 2 int_w^inf Ai(t) dt with w = x e^(i pi/3), by DLMF 9.2.10 and
        //   9.2.12
        bpt = series(x).1;
        apt = q0 - ai_tail(c64(xe, 0.0)).re;
        let t = ai_tail(c64(0.0, xe));
        ant = q1 - 2.0 * t.re;
        bnt = -2.0 * t.im;
    } else {
        let xp6 = 1.0 / (6.0 * PI * xe).sqrt();
        if xp6 == 0.0 {
            // 6 pi xe overflows: the limits as x -> inf (DLMF 9.10.11, 9.7.7, 9.10.12)
            return [q0, f64::INFINITY, q1, 0.0];
        }
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
        if bpt.is_infinite() {
            // exp(xe) may overflow before bpt does
            let g = (0.5 * xe).exp();
            bpt = 2.0 * g * xp6 * su2 * g;
        }
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
/// Calculates the integrals of $\mathrm{Ai}(t)$, $\mathrm{Bi}(t)$, $\mathrm{Ai}(-t)$, and
/// $\mathrm{Bi}(-t)$ for $t$ from 0 to $x$.
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
/// - $\int_0^x \mathrm{Ai}(t) \dd t$
/// - $\int_0^x \mathrm{Bi}(t) \dd t$
/// - $\int_0^x \mathrm{Ai}(-t) \dd t$
/// - $\int_0^x \mathrm{Bi}(-t) \dd t$
///
/// # See also
/// - [`airy`]: Airy functions and derivatives
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
    // Routine called: CEPHES::AIRY for computing Airy functions and
    //                 their derivatives
    // =======================================================

    // NOTE: for `kf` other than 1 or 2 (not reachable from the public API), `rt` stays 0, so
    // `err` is NaN and the Newton loops never terminate (in C++ as well).
    let mut rt = 0.0;
    let mut rt0: f64;
    let mut u = 0.0;
    let mut u1 = 0.0;
    let mut x: f64;
    let (mut ai, mut bi, mut ad, mut bd): (f64, f64, f64, f64);
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
            (_, ai, ad, bi, bd) = cephes::airy(x);

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
            (_, _, ad, _, bd) = cephes::airy(rt);
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
            (_, ai, ad, bi, bd) = cephes::airy(x);

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
            (_, ai, _, bi, _) = cephes::airy(rt);
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

/// The first `nt` zeros of Ai, for a runtime `nt` (see [`airy_ai_zeros`])
pub(crate) fn airy_ai_zeros_vec(nt: usize) -> Vec<f64> {
    let mut x = vec![f64::NAN; nt];
    let [mut xp, mut ai_xp, mut aip_x] = [x.clone(), x.clone(), x.clone()];
    if nt > 0 {
        airyzo(
            nt.try_into().unwrap(),
            1,
            &mut x,
            &mut xp,
            &mut ai_xp,
            &mut aip_x,
        );
    }
    x
}

/// Zeros and values of the Airy function Ai and its derivative
///
/// This corresponds to [`scipy.special.ai_zeros`][ai_zeros] in SciPy.
///
/// [ai_zeros]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ai_zeros.html
///
/// # Returns
///
/// - Roots $x_k$ s.t. $\mathrm{Ai}(x_k) = 0$
/// - Roots $x\'_k$ s.t. $\mathrm{Ai}\'(x\'_k) = 0$
/// - Critical points $\mathrm{Ai}(x\'_k)$, i.e. the values of $\mathrm{Ai}$ where its
///   $\mathrm{Ai}\'$ is zero
/// - Slopes $\mathrm{Ai}\'(x_k)$ at roots $x_k$, i.e. the values of $\mathrm{Ai}\'$ where
///   $\mathrm{Ai}$ is zero
///
/// # See also
/// - [`airy_bi_zeros`]: Zeros of $\mathrm{Bi}$ and $\mathrm{Bi}\'$
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
/// - Roots $x_k$ s.t. $\mathrm{Bi}(x_k) = 0$
/// - Roots $x\'_k$ s.t. $\mathrm{Bi}\'(x\'_k) = 0$
/// - Critical points $\mathrm{Bi}(x\'_k)$, i.e. the values of $\mathrm{Bi}$ where its
///   $\mathrm{Bi}\'$ is zero
/// - Slopes $\mathrm{Bi}\'(x_k)$ at roots $x_k$, i.e. the values of $\mathrm{Bi}\'$ where
///   $\mathrm{Bi}$ is zero
///
/// # See also
/// - [`airy_ai_zeros`]: Zeros of $\mathrm{Ai}$ and $\mathrm{Ai}\'$
#[doc(alias = "airyzo", alias = "bi_zeros")]
#[must_use]
#[inline]
pub fn airy_bi_zeros<const NT: usize>() -> [[f64; NT]; 4] {
    airy_zeros::<NT, 2>()
}
