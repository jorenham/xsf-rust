use core::f64::consts::{PI, SQRT_2};

use num_complex::Complex64;

use crate::utils::mul_exp;
use crate::xsf::amos::{self, ierr_to_sferr, set_error_and_nan};
use crate::xsf::c64::{C64, c64};
use crate::xsf::cephes::dd_real::DoubleDouble;
use crate::xsf::cephes::{self, TWO_THIRDS};

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
    // as AMOS, with -0 as +0 for real `z`
    let z = c64(z.re, z.im + 0.0);
    let zeta = zeta_err(z);
    // whether Ai decays: AMOS evaluates it as sqrt(z) K(zta) for |ph z| <= pi/3, and continues K
    // from -zta beyond that
    let decays = zeta.is_some_and(|(zta, _)| z.re > 0.0 && zta.re >= 0.0);
    // Unlike in xsf, use a Taylor expansion where the series for K(1/3, zeta)
    // and K(2/3, zeta) lose digits by cancellation.
    let taylor = zeta.filter(|(zta, _)| decays && z.re >= 1.0 && zta.norm() <= 2.0);
    let (ai, aip) = if let Some((zta, d)) = taylor {
        let (value, derivative) = ai_taylor(z);
        let (ai, aip) = if kode == 2 {
            let scale = zta.exp() * (1.0 + d);
            (value * scale, derivative * scale)
        } else {
            (value, derivative)
        };
        if z.is_real() {
            (c64(ai.re, 0.0), c64(aip.re, -0.0))
        } else {
            (ai, aip)
        }
    } else {
        // Unlike in xsf, unscale a decaying Ai from `airye` below. AMOS loses up to |zeta| eps
        // in e^-zeta and flushes Ai to 0 below ~1e-305.
        let ai_kode = if decays { 2 } else { kode };
        (amos_airy_or_biry(false, z, 0, ai_kode), amos_airy_or_biry(false, z, 1, ai_kode))
    };
    let bi = amos_airy_or_biry(true, z, 0, kode);
    let bip = amos_airy_or_biry(true, z, 1, kode);
    let Some((zta, d)) = zeta else {
        return [ai, aip, bi, bip];
    };
    // Unlike in xsf, the rounding error d of zeta = zta + d is undone. AMOS evaluates
    // w = sqrt(z) F(zta) and w' = z G(zta) for w = Ai and w = Bi. With d zeta / dz = sqrt(z) and
    // w'' = z w, to first order in d, w += d (w' - w / (2 z)) / sqrt(z) and
    // w' += d (z w - w' / z) / sqrt(z). These are applied as factors of q = w' / w, as `z w` can
    // overflow. If w' does, q ~ +-sqrt(z), as w ~ e^(+-zeta) there (+ where Ai decays).
    //
    // A w that AMOS scaled by e^zta is also off by e^d, which adds c = sqrt(z) to q.
    let s = z.sqrt();
    let q_inf = if decays { s } else { -s };
    let undo = |w: C64, wp: C64, c: C64| {
        let q = Some(wp / w)
            .filter(|q| q.is_finite() && *q != C64::ZERO)
            .unwrap_or(q_inf);
        (w * (1.0 + d * (q - 0.5 / z + c) / s), wp * (1.0 + d * (z / q - 1.0 / z + c) / s))
    };
    let (bi, bip) = undo(bi, bip, C64::ZERO);
    // (a decaying Ai is its `airye` here, which varies slowly)
    let c = if kode == 2 { s } else { C64::ZERO };
    let (ai, aip) = if decays { (ai, aip) } else { undo(ai, aip, c) };
    if kode == 2 {
        // The scaling of Bi by AMOS is off as well, by e^-(|Re zeta| - |Re zta|), where Re zeta
        // can be just across 0 from Re zta.
        let (a, e) = (zta.re.abs(), d.re * zta.re.signum());
        let fb = 1.0 - if e < -a { -2.0 * a - e } else { e };
        return [ai, aip, bi * fb, bip * fb];
    }
    if !decays || taylor.is_some() {
        return [ai, aip, bi, bip];
    }
    // Ai is unscaled from its `airye` by e^-zeta = e^-zta (1 - d).
    let (fa, ez) = ((C64::ONE - d) * C64::cis(-zta.im), (-zta.re).exp());
    // by component if `e^-zta.re` underflows, as `w fa e^-zta.re` may not
    let unscale = |w: C64| {
        let w = w * fa;
        if ez.is_normal() {
            w * ez
        } else {
            c64(mul_exp(w.re, -zta.re), mul_exp(w.im, -zta.re))
        }
    };
    [unscale(ai), unscale(aip), bi, bip]
}

/// Ai and Ai' from the Taylor series at 2, where Ai''(z) = z Ai(z) gives
/// n(n - 1) c_n = 2 c_(n - 2) + c_(n - 3).
fn ai_taylor(z: C64) -> (C64, C64) {
    // Re z >= 1 gives |z - 2| <= |z|. For |zeta| <= 2, |z| < 2.081, and
    // the omitted terms for Ai and Ai' are below 4e-18 in absolute value.
    const COEFFS: [f64; 36] = {
        // Ai(2) and Ai'(2), rounded from Arb.
        let mut c = [0.0; 36];
        c[0] = 0.034_924_130_423_274_38;
        c[1] = -0.053_090_384_433_653_63;
        c[2] = c[0];
        let mut n = 3;
        while n < c.len() {
            c[n] = (2.0 * c[n - 2] + c[n - 3]) / ((n * (n - 1)) as f64);
            n += 1;
        }
        c
    };

    let h = z - 2.0;
    COEFFS
        .iter()
        .rfold((C64::ZERO, C64::ZERO), |(ai, aip), &coefficient| {
            let derivative = aip * h + ai;
            (ai * h + coefficient, derivative)
        })
}

/// AMOS' `zta = 2/3 z sqrt(z)` in double, and `zeta - zta` for the exact `zeta = 2/3 z^(3/2)`
///
/// This is `None` where AMOS doesn't evaluate Ai and Bi as functions of `zta`, i.e. for
/// `|z| <= 1`.
fn zeta_err(z: C64) -> Option<(C64, C64)> {
    if z.norm() <= 1.0 {
        return None;
    }
    let s = z.sqrt();
    // (for real z < 0, AMOS also drops the sign of Re zta = -0)
    let zta = z * s * (2.0 / 3.0);
    let dd = DoubleDouble::new;
    let (zr, zi, sr, si) = (dd(z.re), dd(z.im), dd(s.re), dd(s.im));
    // sqrt(z) = s + (z - s^2) / (2 s)
    let e = c64((zr - sr * sr + si * si).hi, (zi - sr * si * 2.0).hi) / (2.0 * s);
    let ze = z * e;
    let zeta_re = (zr * sr - zi * si + ze.re) * TWO_THIRDS;
    let zeta_im = (zr * si + zi * sr + ze.im) * TWO_THIRDS;
    Some((zta, c64((zeta_re - zta.re).hi, (zeta_im - zta.im).hi)))
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
    // Unlike in xsf, Cephes is used instead of AMOS for x > 10 (and x < -10, scipy/xsf#532), as
    // its double-double zeta makes it more accurate.
    if x > 103.892 && x < 108.0 {
        // Unlike in xsf (scipy/xsf#614), `airye` is unscaled here, as Cephes returns
        // (0, 0, inf, inf) for x > 103.892, while Bi(x) only overflows at x = 104.44 and Ai(x)
        // underflows at x = 107.47.
        let [ai, aip, bi, bip] = airye_f64(x);
        let zeta = cephes::zeta_dd(x, x.sqrt());
        // e^(+-zeta.lo) = 1 +- zeta.lo
        let (ai, aip) = (ai * (1.0 - zeta.lo), aip * (1.0 - zeta.lo));
        let (bi, bip) = (bi * (1.0 + zeta.lo), bip * (1.0 + zeta.lo));
        [
            mul_exp(ai, -zeta.hi),
            mul_exp(aip, -zeta.hi),
            mul_exp(bi, zeta.hi),
            mul_exp(bip, zeta.hi),
        ]
    } else {
        let (_, ai, aip, bi, bip) = cephes::airy(x);
        [ai, aip, bi, bip]
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
        let (fx, gx) = (sum(x, [-2.0, 1.0, -1.0]), sum(0.5 * x * x, [-1.0, 2.0, 1.0]));
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
    let [mut x, mut xp, mut ai_xp, mut aip_x] =
        [[f64::NAN; NT], [f64::NAN; NT], [f64::NAN; NT], [f64::NAN; NT]];
    if NT > 0 {
        airyzo(NT.try_into().unwrap(), KF, &mut x, &mut xp, &mut ai_xp, &mut aip_x);
    }
    [x, xp, ai_xp, aip_x]
}

/// The first `nt` zeros of Ai, for a runtime `nt` (see [`airy_ai_zeros`])
pub(crate) fn airy_ai_zeros_vec(nt: usize) -> Vec<f64> {
    let mut x = vec![f64::NAN; nt];
    let [mut xp, mut ai_xp, mut aip_x] = [x.clone(), x.clone(), x.clone()];
    if nt > 0 {
        airyzo(nt.try_into().unwrap(), 1, &mut x, &mut xp, &mut ai_xp, &mut aip_x);
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
