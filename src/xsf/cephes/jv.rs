//! Translated into pure Rust from `xsf/cephes/jv.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `cbrt` is the one from the Rust standard library, which is more accurate than
//! the Cephes one (see [`crate::cbrt`]), and the `double` -> `int` conversions saturate instead of
//! being undefined behavior (taking the same branches as the x86 conversion in C++).
//!
//! Also unlike in xsf (scipy/xsf#341), the power series returns `-inf` instead of `inf` when it
//! overflows for a negative result, e.g. `jv(-123.57388919180612, 0.028605192190944706)`
//! (`-4.24e431`).

#![allow(clippy::excessive_precision, clippy::unreadable_literal)]

use core::f64::consts::PI;

use crate::utils::{MACHEP, MAXGAM, MAXLOG};
use crate::xsf::cephes::dd_real::frexp_f64;
use crate::xsf::cephes::polevl::polevl;
use crate::xsf::cephes::{airy, j0, j1, lgam_sgn, rgamma};

const JV_BIG: f64 = 1.44115188075855872E+17;

/* Reduce the order by backward recurrence.
 * AMS55 #9.1.27 and 9.1.73.
 */

#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::if_not_else,
    clippy::neg_cmp_op_on_partial_ord
)]
fn jv_recur(n: &mut f64, x: f64, newn: &mut f64, cancel: bool) -> f64 {
    let mut pkm2;
    let mut pkm1;
    let mut pk;
    let mut qkm2;
    let mut qkm1;
    let mut qk;

    /* double pkp1; */
    let mut ans;
    let mut r;
    let mut t;
    let big = JV_BIG;

    /* Continued fraction for Jn(x)/Jn-1(x)
     * AMS 9.1.73
     *
     *    x       -x^2      -x^2
     * ------  ---------  ---------   ...
     * 2 n +   2(n+1) +   2(n+2) +
     *
     * Compute it with the simplest possible algorithm.
     *
     * This continued fraction starts to converge when (|n| + m) > |x|.
     * Hence, at least |x|-|n| iterations are necessary before convergence is
     * achieved. There is a hard limit set below, m <= 30000, which is chosen
     * so that no branch in `jv` requires more iterations to converge.
     * The exact maximum number is (500/3.6)^2 - 500 ~ 19000
     */

    let maxiter = 22000;
    let mut miniter = (x.abs() - n.abs()) as i32;
    if miniter < 1 {
        miniter = 1;
    }

    let mut nflag = i32::from(*n < 0.0);

    // fstart:
    loop {
        pkm2 = 0.0;
        qkm2 = 1.0;
        pkm1 = x;
        qkm1 = *n + *n;
        let xk = -x * x;
        let mut yk = qkm1;
        ans = 0.0; /* ans=0.0 ensures that t=1.0 in the first iteration */
        let mut ctr = 0;
        loop {
            yk += 2.0;
            pk = pkm1 * yk + pkm2 * xk;
            qk = qkm1 * yk + qkm2 * xk;
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;

            /* check convergence */
            if qk != 0.0 && ctr > miniter {
                r = pk / qk;
            } else {
                r = 0.0;
            }

            if r != 0.0 {
                t = ((ans - r) / r).abs();
                ans = r;
            } else {
                t = 1.0;
            }

            ctr += 1;
            if ctr > maxiter {
                break; // goto done
            }
            if t < MACHEP {
                break; // goto done
            }

            /* renormalize coefficients */
            if pk.abs() > big {
                pkm2 /= big;
                pkm1 /= big;
                qkm2 /= big;
                qkm1 /= big;
            }
            if !(t > MACHEP) {
                break;
            }
        }

        // done:
        if ans == 0.0 {
            ans = 1.0;
        }

        /* Change n to n-1 if n < 0 and the continued fraction is small */
        if nflag > 0 && ans.abs() < 0.125 {
            nflag = -1;
            *n -= 1.0;
            continue; // goto fstart
        }
        break;
    }

    let kf = *newn;

    /* backward recurrence
     *              2k
     *  J   (x)  =  --- J (x)  -  J   (x)
     *   k-1         x   k         k+1
     */

    pk = 1.0;
    pkm1 = 1.0 / ans;
    let mut k = *n - 1.0;
    r = 2.0 * k;
    loop {
        pkm2 = (pkm1 * r - pk * x) / x;
        /*      pkp1 = pk; */
        pk = pkm1;
        pkm1 = pkm2;
        r -= 2.0;
        /*
         * t = fabs(pkp1) + fabs(pk);
         * if( (k > (kf + 2.5)) && (fabs(pkm1) < 0.25*t) )
         * {
         * k -= 1.0;
         * t = x*x;
         * pkm2 = ( (r*(r+2.0)-t)*pk - r*x*pkp1 )/t;
         * pkp1 = pk;
         * pk = pkm1;
         * pkm1 = pkm2;
         * r -= 2.0;
         * }
         */
        k -= 1.0;
        if !(k > (kf + 0.5)) {
            break;
        }
    }

    /* Take the larger of the last two iterates
     * on the theory that it may have less cancellation error.
     */

    if cancel && (kf >= 0.0) && (pk.abs() > pkm1.abs()) {
        k += 1.0;
        pkm2 = pk;
    }
    *newn = k;
    pkm2
}

/* Ascending power series for Jv(x).
 * AMS55 #9.1.10.
 */

#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn jv_jvs(n: f64, x: f64) -> f64 {
    let z = -x * x / 4.0;
    let mut u = 1.0;
    let mut y = u;
    let mut k = 1.0;
    let mut t: f64 = 1.0;

    while t > MACHEP {
        u *= z / (k * (n + k));
        y += u;
        k += 1.0;
        if y != 0.0 {
            t = (u / y).abs();
        }
    }
    let (_, mut ex) = frexp_f64(0.5 * x);
    ex = (f64::from(ex) * n) as i32;
    if (ex > -1023) && (ex < 1023) && (n > 0.0) && (n < (MAXGAM - 1.0)) {
        t = (0.5 * x).powf(n) * rgamma(n + 1.0);
        y *= t;
    } else {
        let (lgam, mut sgngam) = lgam_sgn(n + 1.0);
        t = n * (0.5 * x).ln() - lgam;
        if y < 0.0 {
            sgngam = -sgngam;
            y = -y;
        }
        t += y.ln();
        if t < -MAXLOG {
            return 0.0;
        }
        if t > MAXLOG {
            // (unlike in xsf, scipy/xsf#341) keep the sign
            return f64::from(sgngam) * f64::INFINITY;
        }
        y = f64::from(sgngam) * t.exp();
    }
    y
}

/* Hankel's asymptotic expansion
 * for large x.
 * AMS55 #9.2.5.
 */

fn jv_hankel(n: f64, x: f64) -> f64 {
    let m = 4.0 * n * n;
    let mut j = 1.0;
    let z = 8.0 * x;
    let mut k = 1.0;
    let mut p = 1.0;
    let mut u = (m - 1.0) / z;
    let mut q = u;
    let mut sign = 1.0;
    let mut conv = 1.0;
    let mut flag = false;
    let mut t: f64 = 1.0;
    let mut pp = 1.0e38;
    let mut qq = 1.0e38;

    while t > MACHEP {
        k += 2.0;
        j += 1.0;
        sign = -sign;
        u *= (m - k * k) / (j * z);
        p += sign * u;
        k += 2.0;
        j += 1.0;
        u *= (m - k * k) / (j * z);
        q += sign * u;
        t = (u / p).abs();
        if t < conv {
            conv = t;
            qq = q;
            pp = p;
            flag = true;
        }
        /* stop if the terms start getting larger */
        if flag && (t > conv) {
            break; // goto hank1
        }
    }

    // hank1:
    u = x - (0.5 * n + 0.25) * PI;
    (2.0 / (PI * x)).sqrt() * (pp * u.cos() - qq * u.sin())
}

/* Asymptotic expansion for transition region,
 * n large and x close to n.
 * AMS55 #9.3.23.
 */

const JV_PF2: [f64; 2] = [-9.0000000000000000000e-2, 8.5714285714285714286e-2];

const JV_PF3: [f64; 3] = [
    1.3671428571428571429e-1,
    -5.4920634920634920635e-2,
    -4.4444444444444444444e-3,
];

const JV_PF4: [f64; 4] = [
    1.3500000000000000000e-3,
    -1.6036054421768707483e-1,
    4.2590187590187590188e-2,
    2.7330447330447330447e-3,
];

const JV_PG1: [f64; 2] = [-2.4285714285714285714e-1, 1.4285714285714285714e-2];

const JV_PG2: [f64; 3] = [
    -9.0000000000000000000e-3,
    1.9396825396825396825e-1,
    -1.1746031746031746032e-2,
];

const JV_PG3: [f64; 3] = [
    1.9607142857142857143e-2,
    -1.5983694083694083694e-1,
    6.3838383838383838384e-3,
];

#[allow(non_snake_case)]
fn jv_jnt(n: f64, x: f64) -> f64 {
    let cbn = n.cbrt();
    let z = (x - n) / cbn;
    let cbtwo = 2.0_f64.cbrt();

    /* Airy function */
    let mut zz = -cbtwo * z;
    let (_, ai, aip, _, _) = airy(zz);

    /* polynomials in expansion */
    zz = z * z;
    let z3 = zz * z;
    let F = [
        1.0,
        -z / 5.0,
        polevl(z3, &JV_PF2) * zz,
        polevl(z3, &JV_PF3),
        polevl(z3, &JV_PF4) * z,
    ];
    let G = [
        0.3 * zz,
        polevl(z3, &JV_PG1),
        polevl(z3, &JV_PG2) * z,
        polevl(z3, &JV_PG3) * zz,
    ];

    let mut pp = 0.0;
    let mut qq = 0.0;
    let mut nk = 1.0;
    let n23 = (n * n).cbrt();

    for k in 0..=4 {
        let fk = F[k] * nk;
        pp += fk;
        if k != 4 {
            let gk = G[k] * nk;
            qq += gk;
        }
        nk /= n23;
    }

    cbtwo * ai * pp / cbn + 4.0_f64.cbrt() * aip * qq / n
}

/* Asymptotic expansion for large n.
 * AMS55 #9.3.35.
 */

const JV_LAMBDA: [f64; 11] = [
    1.0,
    1.041666666666666666666667E-1,
    8.355034722222222222222222E-2,
    1.282265745563271604938272E-1,
    2.918490264641404642489712E-1,
    8.816272674437576524187671E-1,
    3.321408281862767544702647E+0,
    1.499576298686255465867237E+1,
    7.892301301158651813848139E+1,
    4.744515388682643231611949E+2,
    3.207490090890661934704328E+3,
];

const JV_MU: [f64; 11] = [
    1.0,
    -1.458333333333333333333333E-1,
    -9.874131944444444444444444E-2,
    -1.433120539158950617283951E-1,
    -3.172272026784135480967078E-1,
    -9.424291479571202491373028E-1,
    -3.511203040826354261542798E+0,
    -1.572726362036804512982712E+1,
    -8.228143909718594444224656E+1,
    -4.923553705236705240352022E+2,
    -3.316218568547972508762102E+3,
];

const JV_P1: [f64; 2] = [
    -2.083333333333333333333333E-1,
    1.250000000000000000000000E-1,
];

const JV_P2: [f64; 3] = [
    3.342013888888888888888889E-1,
    -4.010416666666666666666667E-1,
    7.031250000000000000000000E-2,
];

const JV_P3: [f64; 4] = [
    -1.025812596450617283950617E+0,
    1.846462673611111111111111E+0,
    -8.912109375000000000000000E-1,
    7.324218750000000000000000E-2,
];

const JV_P4: [f64; 5] = [
    4.669584423426247427983539E+0,
    -1.120700261622299382716049E+1,
    8.789123535156250000000000E+0,
    -2.364086914062500000000000E+0,
    1.121520996093750000000000E-1,
];

const JV_P5: [f64; 6] = [
    -2.8212072558200244877E1,
    8.4636217674600734632E1,
    -9.1818241543240017361E1,
    4.2534998745388454861E1,
    -7.3687943594796316964E0,
    2.27108001708984375E-1,
];

const JV_P6: [f64; 7] = [
    2.1257013003921712286E2,
    -7.6525246814118164230E2,
    1.0599904525279998779E3,
    -6.9957962737613254123E2,
    2.1819051174421159048E2,
    -2.6491430486951555525E1,
    5.7250142097473144531E-1,
];

const JV_P7: [f64; 8] = [
    -1.9194576623184069963E3,
    8.0617221817373093845E3,
    -1.3586550006434137439E4,
    1.1655393336864533248E4,
    -5.3056469786134031084E3,
    1.2009029132163524628E3,
    -1.0809091978839465550E2,
    1.7277275025844573975E0,
];

#[allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::needless_late_init,
    clippy::similar_names
)]
fn jv_jnx(n: f64, x: f64) -> f64 {
    let mut u = [0.0; 8];

    /* Test for x very close to n. Use expansion for transition region if so. */
    let cbn = n.cbrt();
    let mut z = (x - n) / cbn;
    if z.abs() <= 0.7 {
        return jv_jnt(n, x);
    }

    z = x / n;
    let zz = 1.0 - z * z;
    if zz == 0.0 {
        return 0.0;
    }

    let sz;
    let mut t;
    let zeta;
    let nflg;
    if zz > 0.0 {
        sz = zz.sqrt();
        t = 1.5 * (((1.0 + sz) / z).ln() - sz); /* zeta ** 3/2          */
        zeta = (t * t).cbrt();
        nflg = 1;
    } else {
        sz = (-zz).sqrt();
        t = 1.5 * (sz - (1.0 / z).acos());
        zeta = -(t * t).cbrt();
        nflg = -1;
    }
    let z32i = (1.0 / t).abs();
    let sqz = t.cbrt();

    /* Airy function */
    let n23 = (n * n).cbrt();
    t = n23 * zeta;

    let (_, ai, aip, _, _) = airy(t);

    /* polynomials in expansion */
    u[0] = 1.0;
    let zzi = 1.0 / zz;
    u[1] = polevl(zzi, &JV_P1) / sz;
    u[2] = polevl(zzi, &JV_P2) / zz;
    u[3] = polevl(zzi, &JV_P3) / (sz * zz);
    let mut pp = zz * zz;
    u[4] = polevl(zzi, &JV_P4) / pp;
    u[5] = polevl(zzi, &JV_P5) / (pp * sz);
    pp *= zz;
    u[6] = polevl(zzi, &JV_P6) / pp;
    u[7] = polevl(zzi, &JV_P7) / (pp * sz);

    pp = 0.0;
    let mut qq = 0.0;
    let mut np = 1.0;
    /* flags to stop when terms get larger */
    let mut doa = true;
    let mut dob = true;
    let mut akl = f64::INFINITY;
    let mut bkl = f64::INFINITY;

    for k in 0..=3 {
        let tk = 2 * k;
        let tkp1 = tk + 1;
        let mut zp = 1.0;
        let mut ak = 0.0;
        let mut bk = 0.0;
        for s in 0..=tk {
            if doa {
                let sign = if (s & 3) > 1 { nflg } else { 1 };
                ak += f64::from(sign) * JV_MU[s] * zp * u[tk - s];
            }

            if dob {
                let m = tkp1 - s;
                let sign = if ((m + 1) & 3) > 1 { nflg } else { 1 };
                bk += f64::from(sign) * JV_LAMBDA[s] * zp * u[m];
            }
            zp *= z32i;
        }

        if doa {
            ak *= np;
            t = ak.abs();
            if t < akl {
                akl = t;
                pp += ak;
            } else {
                doa = false;
            }
        }

        if dob {
            bk += JV_LAMBDA[tkp1] * zp * u[0];
            bk *= -np / sqz;
            t = bk.abs();
            if t < bkl {
                bkl = t;
                qq += bk;
            } else {
                dob = false;
            }
        }
        if np < MACHEP {
            break;
        }
        np /= n * n;
    }

    /* normalizing factor ( 4*zeta/(1 - z**2) )**1/4    */
    t = 4.0 * zeta / zz;
    t = t.sqrt().sqrt();

    t *= ai * pp / n.cbrt() + aip * qq / (n23 * n);
    t
}

/// Bessel function of noninteger order
///
/// Returns Bessel function of order v of the argument,
/// where v is real.  Negative x is allowed if v is an integer.
///
/// Several expansions are included: the ascending power
/// series, the Hankel expansion, and two transitional
/// expansions for large v.  If v is not too large, it
/// is reduced by recurrence to a region of best accuracy.
/// The transitional expansions give 12D accuracy for v > 500.
///
/// # Accuracy
///
/// Results for integer v are indicated by *, where x and v
/// both vary from -125 to +125.  Otherwise,
/// x ranges from 0 to 125, v ranges as indicated by "domain."
/// Error criterion is absolute, except relative when |jv()| > 1.
///
/// ```text
/// arithmetic  v domain  x domain    # trials      peak       rms
///    IEEE      0,125     0,125       20000      3.1e-15    2.5e-16
///    IEEE   -125,0       0,125       20000      1.1e-10    8.5e-13
///    IEEE      0,500     0,500       20000      3.4e-15    4.5e-16
/// Integer v:
///    IEEE   -125,125   -125,125      20000      2.9e-15*   2.2e-16*
/// ```
#[allow(
    clippy::cast_possible_truncation,
    clippy::float_cmp,
    clippy::manual_range_contains,
    clippy::too_many_lines
)]
pub(crate) fn jv(mut n: f64, mut x: f64) -> f64 {
    let mut k;
    let mut q;
    let mut t;
    let mut y;

    let mut nint = false; /* Flag for integer n */
    let mut sign = 1; /* Flag for sign inversion */
    let an = n.abs();
    y = an.floor();
    if y == an {
        nint = true;
        let i = (an - 16384.0 * (an / 16384.0).floor()) as i32;
        if n < 0.0 {
            if i & 1 != 0 {
                sign = -sign;
            }
            n = an;
        }
        if x < 0.0 {
            if i & 1 != 0 {
                sign = -sign;
            }
            x = -x;
        }
        if n == 0.0 {
            return j0(x);
        }
        if n == 1.0 {
            return f64::from(sign) * j1(x);
        }
    }

    let y = 'done: {
        if (x < 0.0) && (y != an) {
            break 'done f64::NAN;
        }

        if x == 0.0 && n < 0.0 && !nint {
            return f64::INFINITY * rgamma(n + 1.0);
        }

        y = x.abs();

        if y * y < (n + 1.0).abs() * MACHEP {
            // Unlike in xsf, keep the sign (of odd integer n < 0 or x < 0), e.g.
            // `jv(-5.0, 1e-60)` was 2.6e-304 instead of -2.6e-304.
            return f64::from(sign) * (0.5 * x).powf(n) * rgamma(n + 1.0);
        }

        k = 3.6 * y.sqrt();
        t = 3.6 * an.sqrt();
        if (y < t) && (an > 21.0) {
            return f64::from(sign) * jv_jvs(n, x);
        }
        if (an < k) && (y > 21.0) {
            return f64::from(sign) * jv_hankel(n, x);
        }

        if an < 500.0 {
            /* Note: if x is too large, the continued fraction will fail; but then the
             * Hankel expansion can be used. */
            if nint {
                k = 0.0;
                q = jv_recur(&mut n, x, &mut k, true);
                if k == 0.0 {
                    break 'done j0(x) / q;
                }
                if k == 1.0 {
                    break 'done j1(x) / q;
                }
            }

            if an > 2.0 * y || ((n >= 0.0) && (n < 20.0) && (y > 6.0) && (y < 20.0)) {
                /* Recur backwards from a larger value of n */
                // rlarger:
                k = n;

                y = y + an + 1.0;
                if y < 30.0 {
                    y = 30.0;
                }
                y = n + (y - n).floor();
                q = jv_recur(&mut y, x, &mut k, false);
                break 'done jv_jvs(y, x) * q;
            }

            if k <= 30.0 {
                k = 2.0;
            } else if k < 90.0 {
                k = (3.0 * k) / 4.0;
            }
            if an > (k + 3.0) {
                if n < 0.0 {
                    k = -k;
                }
                q = n - n.floor();
                k = k.floor() + q;
                if n > 0.0 {
                    q = jv_recur(&mut n, x, &mut k, true);
                } else {
                    t = k;
                    k = n;
                    q = jv_recur(&mut t, x, &mut k, true);
                    k = t;
                }
                if q == 0.0 {
                    break 'done 0.0;
                }
            } else {
                k = n;
                q = 1.0;
            }

            /* boundary between convergence of
             * power series and Hankel expansion
             */
            y = k.abs();
            if y < 26.0 {
                t = (0.0083 * y + 0.09) * y + 12.9;
            } else {
                t = 0.9 * y;
            }

            if x > t {
                y = jv_hankel(k, x);
            } else {
                y = jv_jvs(k, x);
            }
            if n > 0.0 {
                y /= q;
            } else {
                y *= q;
            }
        } else {
            /* For large n, use the uniform expansion or the transitional expansion.
             * But if x is of the order of n**2, these may blow up, whereas the
             * Hankel expansion will then work.
             */
            if n < 0.0 {
                break 'done f64::NAN;
            }
            t = x / n;
            t /= n;
            if t > 0.3 {
                y = jv_hankel(n, x);
            } else {
                y = jv_jnx(n, x);
            }
        }
        y
    };

    // done:
    f64::from(sign) * y
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_jv_overflow_sign() {
        assert_eq!(
            super::jv(-123.57388919180612, 0.028605192190944706),
            f64::NEG_INFINITY
        );
        assert_eq!(
            super::jv(-124.68311766689978, 0.032989216721374115),
            f64::INFINITY
        );
        assert_eq!(
            crate::bessel_j(-111.29415504271454, 0.06094218693261655),
            f64::NEG_INFINITY
        );
    }
}
