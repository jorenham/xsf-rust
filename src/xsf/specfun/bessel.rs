//! Bessel functions of integer order, translated into pure Rust from `xsf/specfun/specfun.h`
//! (scipy/xsf@d54dd52).

#![allow(clippy::approx_constant)]

use core::f64::consts::{FRAC_2_PI, FRAC_PI_2, FRAC_PI_4, PI};

use super::SPECFUN_INF;
use crate::utils::EGAMMA;
use crate::xsf::c64::{C64, c64};

fn gam0(x: f64) -> f64 {
    // ================================================
    // Purpose: Compute gamma function Г(x)
    // Input :  x  --- Argument of Г(x)  ( |x| ≤ 1 )
    // Output:  GA --- Г(x)
    // ================================================
    const G: [f64; 25] = [
        1.0e0,
        0.5772156649015329e0,
        -0.6558780715202538e0,
        -0.420026350340952e-1,
        0.1665386113822915e0,
        -0.421977345555443e-1,
        -0.96219715278770e-2,
        0.72189432466630e-2,
        -0.11651675918591e-2,
        -0.2152416741149e-3,
        0.1280502823882e-3,
        -0.201348547807e-4,
        -0.12504934821e-5,
        0.11330272320e-5,
        -0.2056338417e-6,
        0.61160950e-8,
        0.50020075e-8,
        -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
        -0.36968e-11,
        0.51e-12,
        -0.206e-13,
        -0.54e-14,
        0.14e-14,
    ];
    let mut gr = G[24];
    for k in (0..=23).rev() {
        gr = gr * x + G[k];
    }
    1.0 / (gr * x)
}

fn bjndd(x: f64, n: usize, bj: &mut [f64], dj: &mut [f64], fj: &mut [f64]) {
    // =====================================================
    // Purpose: Compute Bessel functions Jn(x) and their
    //          first and second derivatives ( 0 <= n <= 100)
    // Input:   x ---  Argument of Jn(x)  ( x ≥ 0 )
    //          n ---  Order of Jn(x)
    // Output:  BJ(n+1) ---  Jn(x)
    //          DJ(n+1) ---  Jn'(x)
    //          FJ(n+1) ---  Jn"(x)
    // =====================================================

    let mut bs = 0.0;
    let mut f = 0.0;
    let mut f0 = 0.0;
    let mut f1 = 1e-35;

    let m = (1..901)
        .find(|&m| {
            let m = m as f64;
            (0.5 * (6.28 * m).log10() - m * (1.36 * x.abs() / m).log10()) as i32 > 20
        })
        .unwrap_or(900);
    for k in (0..=m).rev() {
        f = 2.0 * (k as f64 + 1.0) * f1 / x - f0;
        if k <= n {
            bj[k] = f;
        }
        if k % 2 == 0 {
            bs += 2.0 * f;
        }
        f0 = f1;
        f1 = f;
    }
    for b in &mut bj[..=n] {
        *b /= bs - f;
    }
    dj[0] = -bj[1];
    fj[0] = -bj[0] - dj[0] / x;
    for k in 1..=n {
        dj[k] = bj[k - 1] - k as f64 * bj[k] / x;
        fj[k] = ((k * k) as f64 / (x * x) - 1.0) * bj[k] - dj[k] / x;
    }
}

/// `cy01`, with the outputs `(zf, zd)` returned instead
fn cy01(kf: i32, z: C64) -> (C64, C64) {
    // ===========================================================
    // Purpose: Compute complex Bessel functions Y0(z), Y1(z)
    //          and their derivatives
    // Input :  z  --- Complex argument of Yn(z) ( n=0,1 )
    //          KF --- Function choice code
    //              KF=0 for ZF=Y0(z) and ZD=Y0'(z)
    //              KF=1 for ZF=Y1(z) and ZD=Y1'(z)
    //              KF=2 for ZF=Y1'(z) and ZD=Y1''(z)
    // Output:  ZF --- Y0(z) or Y1(z) or Y1'(z)
    //          ZD --- Y0'(z) or Y1'(z) or Y1''(z)
    // ===========================================================

    // Unlike in xsf, the coefficients have the exponents of the Fortran original (and of `jynbh`),
    // which the C++ lost from the 3rd or 4th coefficient on, e.g. `a[3] = 0.607...` instead of
    // `6.07...`, so that the zeros with `|z| > 12` had relative errors up to 1e-9.
    #[rustfmt::skip]
    const A: [f64; 12] = [
        -0.703125e-01, 0.112152099609375, -0.5725014209747314, 0.6074042001273483e+01,
        -0.1100171402692467e+03, 0.3038090510922384e+04, -0.1188384262567832e+06,
        0.6252951493434797e+07, -0.4259392165047669e+09, 0.3646840080706556e+11,
        -0.3833534661393944e+13, 0.4854014686852901e+15,
    ];
    #[rustfmt::skip]
    const B: [f64; 12] = [
        0.732421875e-01, -0.2271080017089844, 0.1727727502584457e+01, -0.2438052969955606e+02,
        0.5513358961220206e+03, -0.1825775547429318e+05, 0.8328593040162893e+06,
        -0.5006958953198893e+08, 0.3836255180230433e+10, -0.3649010818849833e+12,
        0.4218971570284096e+14, -0.5827244631566907e+16,
    ];
    #[rustfmt::skip]
    const A1: [f64; 12] = [
        0.1171875, -0.144195556640625, 0.6765925884246826, -0.6883914268109947e+01,
        0.1215978918765359e+03, -0.3302272294480852e+04, 0.1276412726461746e+06,
        -0.6656367718817688e+07, 0.4502786003050393e+09, -0.3833857520742790e+11,
        0.4011838599133198e+13, -0.5060568503314727e+15,
    ];
    #[rustfmt::skip]
    const B1: [f64; 12] = [
        -0.1025390625, 0.2775764465332031, -0.1993531733751297e+01, 0.2724882731126854e+02,
        -0.6038440767050702e+03, 0.1971837591223663e+05, -0.8902978767070678e+06,
        0.5310411010968522e+08, -0.4043620325107754e+10, 0.3827011346598605e+12,
        -0.4406481417852278e+14, 0.6065091351222699e+16,
    ];

    let a0 = z.norm();
    let mut z1 = z;
    let z2 = z * z;
    // The C++ returns early for `a0 == 0`, with the same outputs as below
    let (cby0, cby1, cdy0, cdy1) = if a0 == 0.0 {
        (
            c64(-SPECFUN_INF, 0.0),
            c64(-SPECFUN_INF, 0.0),
            c64(SPECFUN_INF, 0.0),
            c64(SPECFUN_INF, 0.0),
        )
    } else {
        if z.re < 0.0 {
            z1 = -z;
        }

        let (mut cbj0, mut cbj1, mut cby0, mut cby1);
        if a0 <= 12.0 {
            cbj0 = C64::ONE;
            let mut cr = C64::ONE;
            for k in 1..=40 {
                cr = -0.25 * cr * z2 / f64::from(k * k);
                cbj0 += cr;
                if cr.norm() < cbj0.norm() * 1.0e-15 {
                    break;
                }
            }

            cbj1 = C64::ONE;
            cr = C64::ONE;
            for k in 1..=40 {
                cr = -0.25 * cr * z2 / (f64::from(k) * (f64::from(k) + 1.0));
                cbj1 += cr;
                if cr.norm() < cbj1.norm() * 1.0e-15 {
                    break;
                }
            }

            cbj1 *= 0.5 * z1;
            let mut w0 = 0.0;
            cr = C64::ONE;
            let mut cs = C64::ZERO;
            for k in 1..=40 {
                w0 += 1.0 / f64::from(k);
                cr = -0.25 * cr / f64::from(k * k) * z2;
                let cp = cr * w0;
                cs += cp;
                if cp.norm() < cs.norm() * 1.0e-15 {
                    break;
                }
            }

            cby0 = FRAC_2_PI * ((z1 / 2.0).ln() + EGAMMA) * cbj0 - FRAC_2_PI * cs;
            let mut w1 = 0.0;
            cr = C64::ONE;
            cs = C64::ONE;
            for k in 1..=40 {
                w1 += 1.0 / f64::from(k);
                cr = -0.25 * cr / f64::from(k * (k + 1)) * z2;
                let cp = cr * (2.0 * w1 + 1.0 / (f64::from(k) + 1.0));
                cs += cp;
                if cp.norm() < cs.norm() * 1.0e-15 {
                    break;
                }
            }
            cby1 = FRAC_2_PI * (((z1 / 2.0).ln() + EGAMMA) * cbj1 - 1.0 / z1 - 0.25 * z1 * cs);
        } else {
            let mut k0 = 12;
            if a0 >= 35.0 {
                k0 = 10;
            }
            if a0 >= 50.0 {
                k0 = 8;
            }

            let ct1 = z1 - FRAC_PI_4;
            let mut cp0 = C64::ONE;
            for k in 1..=k0 {
                cp0 += A[k as usize - 1] * z1.powi(-2 * k);
            }
            let mut cq0 = -0.125 / z1;
            for k in 1..=k0 {
                cq0 += B[k as usize - 1] * z1.powi(-2 * k - 1);
            }

            let cu = (FRAC_2_PI / z1).sqrt();
            cbj0 = cu * (cp0 * ct1.cos() - cq0 * ct1.sin());
            cby0 = cu * (cp0 * ct1.sin() + cq0 * ct1.cos());

            let ct2 = z1 - 0.75 * PI;
            let mut cp1 = C64::ONE;
            for k in 1..=k0 {
                cp1 += A1[k as usize - 1] * z1.powi(-2 * k);
            }

            let mut cq1 = 0.375 / z1;
            for k in 1..=k0 {
                cq1 += B1[k as usize - 1] * z1.powi(-2 * k - 1);
            }
            cbj1 = cu * (cp1 * ct2.cos() - cq1 * ct2.sin());
            cby1 = cu * (cp1 * ct2.sin() + cq1 * ct2.cos());
        }

        if z.re < 0.0 {
            if z.im < 0.0 {
                cby0 -= 2.0 * C64::I * cbj0;
            }
            if z.im > 0.0 {
                cby0 += 2.0 * C64::I * cbj0;
            }
            if z.im < 0.0 {
                cby1 = -(cby1 - 2.0 * C64::I * cbj1);
            }
            if z.im > 0.0 {
                cby1 = -(cby1 + 2.0 * C64::I * cbj1);
            }
            // (the C++ also negates `cbj1` here, which is not used after)
        }

        (cby0, cby1, -cby1, cby0 - 1.0 / z * cby1)
    };

    if kf == 0 {
        (cby0, cdy0)
    } else if kf == 1 {
        (cby1, cdy1)
    } else {
        // kf == 2
        (cdy1, -cdy1 / z - (1.0 - 1.0 / (z * z)) * cby1)
    }
}

/// `cyzo`; `nt` is the length of `zo` and `zv`
pub(crate) fn cyzo(kf: i32, kc: i32, zo: &mut [C64], zv: &mut [C64]) {
    // ===========================================================
    // Purpose : Compute the complex zeros of Y0(z), Y1(z) and
    //           Y1'(z), and their associated values at the zeros
    //           using the modified Newton's iteration method
    // Input:    NT --- Total number of zeros/roots
    //           KF --- Function choice code
    //                  KF=0 for  Y0(z) & Y1(z0)
    //                  KF=1 for  Y1(z) & Y0(z1)
    //                  KF=2 for  Y1'(z) & Y1(z1')
    //           KC --- Choice code
    //                  KC=0 for complex roots
    //                  KC=1 for real roots
    // Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
    //           ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
    //                     at the L-th zero
    // Routine called: CY01 for computing Y0(z) and Y1(z), and
    //                 their derivatives
    // ===========================================================

    let nt = zo.len();
    let mut x = 0.0;
    let mut y = 0.0;
    let mut h = 0.0;

    if kc == 0 {
        x = -2.4;
        y = 0.54;
        h = 3.14;
    } else if kc == 1 {
        x = 0.89;
        y = 0.0;
        h = -3.14;
    }

    if kf == 1 {
        x = -0.503;
    }

    if kf == 2 {
        x = 0.577;
    }
    let mut z = c64(x, y);
    let mut w = 0.0;
    for nr in 1..=nt {
        if nr > 1 {
            z = zo[nr - 2] - h;
        }
        let mut it = 0;
        loop {
            it += 1;
            let (zf, zd) = cy01(kf, z);
            // Unlike in xsf, the Newton step for `f(z) / prod(z - z_i)`, with the zeros `z_i` found
            // so far, is `f / (f' - f sum 1 / (z - z_i))` instead of the quotient of the products,
            // which overflow from the 140th zero on (that zero was inaccurate, the next ones NaN).
            let mut zq = C64::ZERO;
            for &zi in &zo[..nr - 1] {
                zq += 1.0 / (z - zi);
            }
            z -= zf / (zd - zq * zf);
            let w0 = w;
            w = z.norm();
            if !(it <= 50 && ((w - w0) / w).abs() > 1.0e-12) {
                break;
            }
        }

        zo[nr - 1] = z;
    }

    for (&z, v) in zo.iter().zip(zv.iter_mut()) {
        if kf == 0 || kf == 2 {
            *v = cy01(1, z).0;
        } else if kf == 1 {
            *v = cy01(0, z).0;
        }
    }
}

/// `jdzo`; `zo` needs `nt + 25` elements (`n`, `m` and `p` `nt + 24`), and `p` is `true` for 1 (TE)
pub(crate) fn jdzo(nt: i32, zo: &mut [f64], n: &mut [u32], m: &mut [u32], p: &mut [bool]) {
    // ===========================================================
    // Purpose: Compute the zeros of Bessel functions Jn(x) and
    //          Jn'(x), and arrange them in the order of their
    //          magnitudes
    // Input :  NT    --- Number of total zeros ( NT ≤ 1200 )
    // Output:  ZO(L) --- Value of the L-th zero of Jn(x)
    //                    and Jn'(x)
    //          N(L)  --- n, order of Jn(x) or Jn'(x) associated
    //                    with the L-th zero
    //          M(L)  --- m, serial number of the zeros of Jn(x)
    //                    or Jn'(x) associated with the L-th zero
    //                    ( L is the serial number of all the
    //                      zeros of Jn(x) and Jn'(x) )
    //          P(L)  --- 0 (TM) or 1 (TE), a code for designating the
    //                    zeros of Jn(x)  or Jn'(x).
    //                    In the waveguide applications, the zeros
    //                    of Jn(x) correspond to TM modes and
    //                    those of Jn'(x) correspond to TE modes
    //
    // Routine called:    BJNDD for computing Jn(x), Jn'(x) and
    //                    Jn''(x)
    // =============================================================

    // Unlike in xsf, the work arrays are on the stack, as in specfun.f (`m1`, `n1` and `zoc` are
    // `mnzoc`, and `bj`, `dj` and `fj` are `bdfj` in C++), so there's no `Status::NoMemory`.
    // Note: ZO and ZOC arrays are 0-indexed in specfun.f
    let mut p1 = [false; 70];
    let mut m1 = [0; 70];
    let mut n1 = [0; 70];
    let mut zoc = [0.0; 71];
    let mut bj = [0.0; 101];
    let mut dj = [0.0; 101];
    let mut fj = [0.0; 101];

    let mut x = 0.0;
    let ntf = f64::from(nt);

    let (xm, nm, mm) = if nt < 600 {
        (
            -1.0 + 2.248485 * ntf.sqrt() - 0.0159382 * ntf + 3.208775e-4 * ntf.powf(1.5),
            (14.5 + 0.05875 * ntf) as usize,
            (0.02 * ntf) as usize + 6,
        )
    } else {
        (
            5.0 + 1.445389 * ntf.sqrt() + 0.01889876 * ntf - 2.147763e-4 * ntf.powf(1.5),
            (27.8 + 0.0327 * ntf) as usize,
            (0.01088 * ntf) as usize + 10,
        )
    };
    // Unlike in xsf, all orders `n <= xm` are searched (and `j'_{n,1} > n` for the others), since
    // the `nm` above misses some, e.g. `j'_{26,1}` for `nt = 210`.
    let nm = nm.max(xm as usize + 1);

    let mut l0 = 0;
    /* 45 LOOP */
    for i in 1..=nm {
        let i1 = (i - 1) as f64;
        let mut x1 = 0.407658 + 0.4795504 * i1.sqrt() + 0.983618 * i1;
        let mut x2 = 1.99535 + 0.8333883 * i1.sqrt() + 0.984584 * i1;
        let mut l1 = 0;

        /* 30 LOOP */
        'l30: for j in 1..=mm {
            'l20: {
                if i != 1 || j != 1 {
                    x = x1;
                    loop {
                        bjndd(x, i, &mut bj, &mut dj, &mut fj);
                        let x0 = x;
                        x -= dj[i - 1] / fj[i - 1];
                        // Unlike in xsf, `x` is checked instead of the initial guess `x1` (as for
                        // the zeros of Jn below), since `x1 > xm` misses e.g. `j'_{56,2}` for
                        // `nt = 1121`.
                        if x > xm {
                            break 'l20;
                        }
                        if !((x - x0).abs() > 1e-10) {
                            break;
                        }
                    }
                }
                /* 15 */
                l1 += 1;
                n1[l1 - 1] = i as u32 - 1;
                m1[l1 - 1] = j as u32;
                if i == 1 {
                    m1[l1 - 1] = j as u32 - 1;
                }
                p1[l1 - 1] = true;
                zoc[l1] = x;
                if i <= 15 {
                    x1 =
                        x + 3.057 + 0.0122 * i1 + (1.555 + 0.41575 * i1) / ((j + 1) as f64).powi(2);
                } else {
                    x1 =
                        x + 2.918 + 0.01924 * i1 + (6.26 + 0.13205 * i1) / ((j + 1) as f64).powi(2);
                }
            }
            // L20:
            x = x2;
            loop {
                bjndd(x, i, &mut bj, &mut dj, &mut fj);
                let x0 = x;
                x -= bj[i - 1] / dj[i - 1];
                if x > xm {
                    continue 'l30;
                }
                if !((x - x0).abs() > 1e-10) {
                    break;
                }
            }
            l1 += 1;
            n1[l1 - 1] = i as u32 - 1;
            m1[l1 - 1] = j as u32;
            p1[l1 - 1] = false;
            zoc[l1] = x;
            if i <= 15 {
                x2 = x + 3.11 + 0.0138 * i1 + (0.04832 + 0.2804 * i1) / ((j + 1) as f64).powi(2);
            } else {
                x2 = x + 3.001 + 0.0105 * i1 + (11.52 + 0.48525 * i1) / ((j + 3) as f64).powi(2);
            }
        }
        let l = l0 + l1;
        loop {
            if l0 == 0 {
                p[..l].copy_from_slice(&p1[..l]);
                m[..l].copy_from_slice(&m1[..l]);
                n[..l].copy_from_slice(&n1[..l]);
                zo[1..=l].copy_from_slice(&zoc[1..=l]);
                l1 = 0;
            } else if zo[l0] >= zoc[l1] {
                p[l0 + l1 - 1] = p[l0 - 1];
                m[l0 + l1 - 1] = m[l0 - 1];
                n[l0 + l1 - 1] = n[l0 - 1];
                zo[l0 + l1] = zo[l0];
                l0 -= 1;
            } else {
                p[l0 + l1 - 1] = p1[l1 - 1];
                m[l0 + l1 - 1] = m1[l1 - 1];
                n[l0 + l1 - 1] = n1[l1 - 1];
                zo[l0 + l1] = zoc[l1];
                l1 -= 1;
            }
            if l1 == 0 {
                break;
            }
        }
        /* 45 */
        l0 = l;
    }
}

pub(crate) fn jynb(
    n: i32,
    x: f64,
    nm: &mut i32,
    bj: &mut [f64],
    dj: &mut [f64],
    by: &mut [f64],
    dy: &mut [f64],
) {
    // =====================================================
    // Purpose: Compute Bessel functions Jn(x), Yn(x) and
    //          their derivatives
    // Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
    //          n --- Order of Jn(x) and Yn(x)
    // Output:  BJ(n) --- Jn(x)
    //          DJ(n) --- Jn'(x)
    //          BY(n) --- Yn(x)
    //          DY(n) --- Yn'(x)
    //          NM --- Highest order computed
    // Routines called:
    //          JYNBH to calculate the Jn and Yn
    // =====================================================

    jynbh(n, 0, x, nm, bj, by);
    // Compute derivatives by differentiation formulas
    if x < 1.0e-100 {
        for k in 0..=n {
            dj[k as usize] = 0.0;
            dy[k as usize] = SPECFUN_INF;
        }
        dj[1] = 0.5;
    } else {
        dj[0] = -bj[1];
        for k in 1..=*nm {
            dj[k as usize] = bj[(k - 1) as usize] - f64::from(k) / x * bj[k as usize];
        }

        dy[0] = -by[1];
        for k in 1..=*nm {
            dy[k as usize] = by[(k - 1) as usize] - f64::from(k) * by[k as usize] / x;
        }
    }
}

fn jynbh(n: i32, nmin: i32, x: f64, nm: &mut i32, bj: &mut [f64], by: &mut [f64]) {
    // =====================================================
    // Purpose: Compute Bessel functions Jn(x), Yn(x)
    // Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
    //          n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
    //          nmin -- Lowest order computed  ( nmin ≥ 0 )
    // Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
    //          BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 to calculate the starting
    //          point for backward recurrence
    // =====================================================

    let r2p = 0.63661977236758;
    let mut bj0: f64;
    let mut bj1: f64;
    let mut by0: f64;
    let mut by1: f64;

    let a: [f64; 4] = [
        -0.0703125,
        0.112152099609375,
        -0.5725014209747314,
        0.6074042001273483e+01,
    ];
    let b: [f64; 4] = [
        0.0732421875,
        -0.2271080017089844,
        0.1727727502584457e+01,
        -0.2438052969955606e+02,
    ];
    let a1: [f64; 4] = [
        0.1171875,
        -0.144195556640625,
        0.6765925884246826,
        -0.6883914268109947e+01,
    ];
    let b1: [f64; 4] = [
        -0.1025390625,
        0.2775764465332031,
        -0.1993531733751297e+01,
        0.2724882731126854e+02,
    ];

    *nm = n;
    if x < 1.0e-100 {
        for k in nmin..=n {
            bj[(k - nmin) as usize] = 0.0;
            by[(k - nmin) as usize] = -SPECFUN_INF;
        }

        if nmin == 0 {
            bj[0] = 1.0;
        }
        return;
    }

    let ky: i32 = if (x <= 300.0) || (f64::from(n) > 0.9 * x) {
        // Backward recurrence for Jn
        if n == 0 {
            *nm = 1;
        }
        let mut m = msta1(x, 200);
        if m < *nm {
            *nm = m;
        } else {
            m = msta2(x, *nm, 15);
        }
        let mut bs = 0.0;
        let mut su = 0.0;
        let mut sv = 0.0;
        let mut f2 = 0.0;
        let mut f1 = 1.0e-100;
        let mut f = 0.0;

        for k in (0..=m).rev() {
            f = 2.0 * (f64::from(k) + 1.0) / x * f1 - f2;
            if (k <= *nm) && (k >= nmin) {
                bj[(k - nmin) as usize] = f;
            }
            if k == 2 * (k / 2) && k != 0 {
                bs += 2.0 * f;
                su += (-1.0_f64).powf(f64::from(k / 2)) * f / f64::from(k);
            } else if k > 1 {
                sv += (-1.0_f64).powf(f64::from(k / 2)) * f64::from(k)
                    / (f64::from(k) * f64::from(k) - 1.0)
                    * f;
            }
            f2 = f1;
            f1 = f;
        }
        let s0 = bs + f;

        for k in nmin..=*nm {
            bj[(k - nmin) as usize] /= s0;
        }
        // Estimates for Yn at start of recurrence
        bj0 = f1 / s0;
        bj1 = f2 / s0;
        let ec = (x / 2.0).ln() + 0.5772156649015329;
        by0 = r2p * (ec * bj0 - 4.0 * su / s0);
        by1 = r2p * ((ec - 1.0) * bj1 - bj0 / x - 4.0 * sv / s0);

        if 0 >= nmin {
            by[(0 - nmin) as usize] = by0;
        }
        if 1 >= nmin {
            by[(1 - nmin) as usize] = by1;
        }
        2
    } else {
        // Hankel expansion
        let t1 = x - FRAC_PI_4;
        let mut p0 = 1.0;
        let mut q0 = -0.125 / x;

        for k in 1..=4 {
            p0 += a[(k - 1) as usize] * x.powf(f64::from(-2 * k));
            q0 += b[(k - 1) as usize] * x.powf(f64::from(-2 * k - 1));
        }

        let cu = (r2p / x).sqrt();
        bj0 = cu * (p0 * t1.cos() - q0 * t1.sin());
        by0 = cu * (p0 * t1.sin() + q0 * t1.cos());

        if 0 >= nmin {
            bj[(0 - nmin) as usize] = bj0;
            by[(0 - nmin) as usize] = by0;
        }

        let t2 = x - 0.75 * PI;
        let mut p1 = 1.0;
        let mut q1 = 0.375 / x;

        for k in 1..=4 {
            p1 += a1[(k - 1) as usize] * x.powf(f64::from(-2 * k));
            q1 += b1[(k - 1) as usize] * x.powf(f64::from(-2 * k - 1));
        }

        bj1 = cu * (p1 * t2.cos() - q1 * t2.sin());
        by1 = cu * (p1 * t2.sin() + q1 * t2.cos());

        if 1 >= nmin {
            bj[(1 - nmin) as usize] = bj1;
            by[(1 - nmin) as usize] = by1;
        }

        for k in 2..=*nm {
            let bjk = 2.0 * (f64::from(k) - 1.0) / x * bj1 - bj0;
            if k >= nmin {
                bj[(k - nmin) as usize] = bjk;
            }
            bj0 = bj1;
            bj1 = bjk;
        }
        2
    };
    // Forward recurrence for Yn
    for k in ky..=*nm {
        let byk = 2.0 * (f64::from(k) - 1.0) * by1 / x - by0;

        if k >= nmin {
            by[(k - nmin) as usize] = byk;
        }

        by0 = by1;
        by1 = byk;
    }
}

fn jyndd(
    n: i32,
    x: f64,
    bjn: &mut f64,
    djn: &mut f64,
    fjn: &mut f64,
    byn: &mut f64,
    dyn_: &mut f64,
    fyn: &mut f64,
) {
    // ===========================================================
    // purpose: compute bessel functions jn(x) and yn(x), and
    //          their first and second derivatives
    // input:   x   ---  argument of jn(x) and yn(x) ( x > 0 )
    //          n   ---  order of jn(x) and yn(x)
    // output:  bjn ---  jn(x)
    //          djn ---  jn'(x)
    //          fjn ---  jn"(x)
    //          byn ---  yn(x)
    //          dyn ---  yn'(x)
    //          fyn ---  yn"(x)
    // routines called:
    //          jynbh to compute jn and yn
    // ===========================================================

    // NOTE: in C++ `bj` and `by` are uninitialized, and partly read uninitialized if
    // msta1(x, 200) <= n (x << n). jyzo only calls this near the zeros (x > n), so that's fine.
    let mut nm = 0;
    let mut bj = [0.0; 2];
    let mut by = [0.0; 2];

    jynbh(n + 1, n, x, &mut nm, &mut bj, &mut by);
    // compute derivatives by differentiation formulas
    *bjn = bj[0];
    *byn = by[0];
    *djn = -bj[1] + f64::from(n) * bj[0] / x;
    *dyn_ = -by[1] + f64::from(n) * by[0] / x;
    *fjn = (f64::from(n) * f64::from(n) / (x * x) - 1.0) * (*bjn) - (*djn) / x;
    *fyn = (f64::from(n) * f64::from(n) / (x * x) - 1.0) * (*byn) - (*dyn_) / x;
}

#[allow(non_snake_case)]
pub(crate) fn jyzo(
    n: i32,
    nt: i32,
    rj0: &mut [f64],
    rj1: &mut [f64],
    ry0: &mut [f64],
    ry1: &mut [f64],
) {
    // ======================================================
    // Purpose: Compute the zeros of Bessel functions Jn(x),
    //          Yn(x), and their derivatives
    // Input :  n  --- Order of Bessel functions  (n >= 0)
    //          NT --- Number of zeros (roots)
    // Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
    //          RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
    //          RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
    //          RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
    // Routine called: JYNDD for computing Jn(x), Yn(x), and
    //                 their first and second derivatives
    // ======================================================

    /*
     * SciPy Note:
     * See GH-18859 for additional changes done by SciPy for
     * better initial condition selection in Newton iteration
     */

    // NOTE: each of rj0, rj1, ry0, ry1 gets max(nt, 1) elements written (index 0 is written even
    // if nt < 1), and only rj0[0..nt] / ry0[0..nt] are read, so buffers of length max(nt, 1) are
    // in bounds. The caller (bessel_zeros) returns early for N == 0.
    // NOTE: like in C++, the Newton iterations loop forever if they don't converge.
    let mut L: i32;
    let mut x: f64;
    let mut bjn = 0.0;
    let mut djn = 0.0;
    let mut fjn = 0.0;
    let mut byn = 0.0;
    let mut dyn_ = 0.0;
    let mut fyn = 0.0;
    // -- Newton method for j_{N,L}
    // initial guess for j_{N,1}
    if n == 0 {
        x = 2.4;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 1.85576 * f64::from(n).powf(0.33333)
            + 1.03315 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l10: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= bjn / djn;
        if (x - x0).abs() > 1e-11 {
            continue 'l10;
        }

        L += 1;
        rj0[(L - 1) as usize] = x;
        // initial guess for j_{N,L+1}
        if L == 1 {
            if n == 0 {
                x = 5.52;
            } else {
                // Expansion from https://dlmf.nist.gov/10.21#E32 and
                // coefficients from Olver 1951
                x = f64::from(n)
                    + 3.24460 * f64::from(n).powf(0.33333)
                    + 3.15824 / f64::from(n).powf(0.33333);
            }
        } else {
            // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
            x = rj0[(L - 1) as usize] + (rj0[(L - 1) as usize] - rj0[(L - 2) as usize]);
        }
        // NOTE: `n + 10` overflows i32 for n > i32::MAX - 10 (C++ UB; Rust panics in debug builds).
        if L <= (n + 10) {
            jyndd(
                n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
            );
            let h = djn.abs().atan2((fjn * bjn).abs().sqrt());
            // NOTE: `atan(h)` here, but `tan(h)` in the corresponding step for y_{N,L} below;
            // reproduced as is.
            let b = -djn / (bjn * h.atan());
            x -= (h - FRAC_PI_2) / b;
        }

        if L < nt {
            continue 'l10;
        }
        break;
    }

    // -- Newton method for j_{N,L+1}'
    if n == 0 {
        x = 3.8317;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 0.80861 * f64::from(n).powf(0.33333)
            + 0.07249 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l15: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= djn / fjn;
        if (x - x0).abs() > 1e-11 {
            continue 'l15;
        }
        L += 1;
        rj1[(L - 1) as usize] = x;
        if L < nt {
            // https://dlmf.nist.gov/10.21#E20
            x = rj1[(L - 1) as usize] + (rj0[L as usize] - rj0[(L - 1) as usize]);
            continue 'l15;
        }
        break;
    }

    // -- Newton method for y_{N,L}
    // initial guess for y_{N,1}
    if n == 0 {
        x = 0.89357697;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 0.93158 * f64::from(n).powf(0.33333)
            + 0.26035 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l20: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= byn / dyn_;
        if (x - x0).abs() > 1.0e-11 {
            continue 'l20;
        }
        L += 1;
        ry0[(L - 1) as usize] = x;
        // initial guess for y_{N,L+1}
        if L == 1 {
            if n == 0 {
                x = 3.957678419314858;
            } else {
                // Expansion from https://dlmf.nist.gov/10.21#E33 and
                // coefficients from Olver 1951
                x = f64::from(n)
                    + 2.59626 * f64::from(n).powf(0.33333)
                    + 2.022183 / f64::from(n).powf(0.33333);
            }
        } else {
            // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
            x = ry0[(L - 1) as usize] + (ry0[(L - 1) as usize] - ry0[(L - 2) as usize]);
        }
        if L <= n + 10 {
            jyndd(
                n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
            );
            let h = dyn_.abs().atan2((fyn * byn).abs().sqrt());
            let b = -dyn_ / (byn * h.tan());
            x -= (h - FRAC_PI_2) / b;
        }

        if L < nt {
            continue 'l20;
        }
        break;
    }

    // -- Newton method for y_{N,L+1}'
    if n == 0 {
        x = 2.67257;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = f64::from(n)
            + 1.8211 * f64::from(n).powf(0.33333)
            + 0.94001 / f64::from(n).powf(0.33333);
    }
    // iterate
    L = 0;
    'l25: loop {
        let x0 = x;
        jyndd(
            n, x, &mut bjn, &mut djn, &mut fjn, &mut byn, &mut dyn_, &mut fyn,
        );
        x -= dyn_ / fyn;
        if (x - x0).abs() > 1.0e-11 {
            continue 'l25;
        }
        L += 1;
        ry1[(L - 1) as usize] = x;
        if L < nt {
            // https://dlmf.nist.gov/10.21#E20
            x = ry1[(L - 1) as usize] + (ry0[L as usize] - ry0[(L - 1) as usize]);
            continue 'l25;
        }
        break;
    }
}

pub(crate) fn lamn(n: i32, x: f64, nm: &mut i32, bl: &mut [f64], dl: &mut [f64]) {
    // =========================================================
    // Purpose: Compute lambda functions and their derivatives
    // Input:   x --- Argument of lambda function
    //          n --- Order of lambda function
    // Output:  BL(n) --- Lambda function of order n
    //          DL(n) --- Derivative of lambda function
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 for computing the start
    //          point for backward recurrence
    // =========================================================

    // NOTE: bl[0..=n] and dl[0..=n] can be written, and also dl[1] (for |x| < 1e-100) or bl[1] /
    // dl[1] (for x > 12, where *nm is 1 if n == 0), which is out of bounds for n == 0 with buffers
    // of length n + 1. The caller (jahnke_emden_lambda) passes n >= 1, so that's in bounds.
    *nm = n;
    if x.abs() < 1e-100 {
        for k in 0..=n {
            bl[k as usize] = 0.0;
            dl[k as usize] = 0.0;
        }
        bl[0] = 1.0;
        dl[1] = 0.5;
        return;
    }
    if x <= 12.0 {
        let x2 = x * x;
        for k in 0..=n {
            let mut bk = 1.0;
            let mut r = 1.0;
            for i in 1..=50 {
                r = -0.25 * r * x2 / (f64::from(i) * f64::from(i + k));
                bk += r;

                if r.abs() < bk.abs() * 1.0e-15 {
                    break;
                }
            }
            bl[k as usize] = bk;
            if k >= 1 {
                dl[(k - 1) as usize] = -0.5 * x / f64::from(k) * bk;
            }
        }
        let mut uk = 1.0;
        let mut r = 1.0;
        for i in 1..=50 {
            r = -0.25 * r * x2 / (f64::from(i) * (f64::from(i + n) + 1.0));
            uk += r;

            if r.abs() < uk.abs() * 1.0e-15 {
                break;
            }
        }
        dl[n as usize] = -0.5 * x / (f64::from(n) + 1.0) * uk;
        return;
    }
    if n == 0 {
        *nm = 1;
    }
    let mut m = msta1(x, 200);
    if m < *nm {
        *nm = m;
    } else {
        m = msta2(x, *nm, 15);
    }
    let mut bs = 0.0;
    let mut f = 0.0;
    let mut f0 = 0.0;
    let mut f1 = 1e-100;
    for k in (0..=m).rev() {
        f = 2.0 * (f64::from(k) + 1.0) * f1 / x - f0;
        if k <= *nm {
            bl[k as usize] = f;
        }
        if k % 2 == 0 {
            bs += 2.0 * f;
        }
        f0 = f1;
        f1 = f;
    }
    let bg = bs - f;
    for k in 0..=*nm {
        bl[k as usize] /= bg;
    }
    let mut r0 = 1.0;
    for k in 1..=*nm {
        r0 = 2.0 * r0 * f64::from(k) / x;
        bl[k as usize] *= r0;
    }
    dl[0] = -0.5 * x * bl[1];
    for k in 1..=*nm {
        dl[k as usize] = 2.0 * f64::from(k) / x * (bl[(k - 1) as usize] - bl[k as usize]);
    }
}

pub(crate) fn lamv(v: f64, x: f64, vm: &mut f64, vl: &mut [f64], dl: &mut [f64]) {
    // =========================================================
    // Purpose: Compute lambda function with arbitrary order v,
    //          and their derivative
    // Input :  x --- Argument of lambda function
    //          v --- Order of lambda function
    // Output:  VL(n) --- Lambda function of order n+v0
    //          DL(n) --- Derivative of lambda function
    //          VM --- Highest order computed
    // Routines called:
    //      (1) MSTA1 and MSTA2 for computing the starting
    //          point for backward recurrence
    //      (2) GAM0 for computing gamma function (|x| ≤ 1)
    // =========================================================

    // NOTE: vl[0..=n] and dl[0..=n] are written with n = (int)v, and for x > 12 also vl[1] / dl[1]
    // (out of bounds for v < 1 with buffers of length n + 1). The caller
    // (jahnke_emden_lambda) passes 1 <= v < 2^31 with buffers of length n + 1, so that's in bounds,
    // and `(int)v` is not UB.
    let rp2 = 0.63661977236758;

    let x = x.abs();
    let x2 = x * x;
    let mut n = v as i32;
    let v0 = v - f64::from(n);
    *vm = v;

    if x <= 12.0 {
        for k in 0..=n {
            let vk = v0 + f64::from(k);
            let mut bk = 1.0;
            let mut r = 1.0;

            for i in 1..=50 {
                r = -0.25 * r * x2 / (f64::from(i) * (f64::from(i) + vk));
                bk += r;

                if r.abs() < bk.abs() * 1.0e-15 {
                    break;
                }
            }
            vl[k as usize] = bk;

            let mut uk = 1.0;
            r = 1.0;

            for i in 1..=50 {
                r = -0.25 * r * x2 / (f64::from(i) * (f64::from(i) + vk + 1.0));
                uk += r;

                if r.abs() < uk.abs() * 1.0e-15 {
                    break;
                }
            }
            dl[k as usize] = -0.5 * x / (vk + 1.0) * uk;
        }
        return;
    }

    let k0 = if x >= 50.0 {
        8
    } else if x >= 35.0 {
        10
    } else {
        11
    };
    let mut bjv0 = 0.0;
    let mut bjv1 = 0.0;

    for j in 0..=1 {
        let vv = 4.0 * (f64::from(j) + v0) * (f64::from(j) + v0);
        let mut px = 1.0;
        let mut rp = 1.0;
        for k in 1..=k0 {
            // NOTE: `pow(y, 2.0)` is written as `y * y`, which is what GCC compiles it to.
            rp = -0.78125e-2
                * rp
                * (vv - (4.0 * f64::from(k) - 3.0) * (4.0 * f64::from(k) - 3.0))
                * (vv - (4.0 * f64::from(k) - 1.0) * (4.0 * f64::from(k) - 1.0))
                / (f64::from(k) * (2.0 * f64::from(k) - 1.0) * x2);
            px += rp;
        }

        let mut qx = 1.0;
        rp = 1.0;
        for k in 1..=k0 {
            rp = -0.78125e-2
                * rp
                * (vv - (4.0 * f64::from(k) - 1.0) * (4.0 * f64::from(k) - 1.0))
                * (vv - (4.0 * f64::from(k) + 1.0) * (4.0 * f64::from(k) + 1.0))
                / (f64::from(k) * (2.0 * f64::from(k) + 1.0) * x2);
            qx += rp;
        }

        qx = 0.125 * (vv - 1.0) * qx / x;
        let xk = x - (0.5 * (f64::from(j) + v0) + 0.25) * PI;
        let a0 = (rp2 / x).sqrt();
        let ck = xk.cos();
        let sk = xk.sin();

        if j == 0 {
            bjv0 = a0 * (px * ck - qx * sk);
        }
        if j == 1 {
            bjv1 = a0 * (px * ck - qx * sk);
        }
    }

    let mut ga;
    if v0 == 0.0 {
        ga = 1.0;
    } else {
        ga = gam0(v0);
        ga *= v0;
    }

    let fac = (2.0 / x).powf(v0) * ga;
    vl[0] = bjv0;
    dl[0] = -bjv1 + v0 / x * bjv0;
    vl[1] = bjv1;
    dl[1] = bjv0 - (1.0 + v0) / x * bjv1;
    let mut r0 = 2.0 * (1.0 + v0) / x;

    if n <= 1 {
        vl[0] *= fac;
        dl[0] = fac * dl[0] - v0 / x * vl[0];
        vl[1] *= fac * r0;
        dl[1] = fac * r0 * dl[1] - (1.0 + v0) / x * vl[1];
        return;
    }

    if n >= 2 && f64::from(n) <= 0.9 * x {
        let mut f0 = bjv0;
        let mut f1 = bjv1;

        for k in 2..=n {
            let f = 2.0 * (f64::from(k) + v0 - 1.0) / x * f1 - f0;
            f0 = f1;
            f1 = f;
            vl[k as usize] = f;
        }
    } else if n >= 2 {
        let mut m = msta1(x, 200);
        if m < n {
            n = m;
        } else {
            m = msta2(x, n, 15);
        }

        let mut f = 0.0;
        let mut f2 = 0.0;
        let mut f1 = 1.0e-100;

        for k in (0..=m).rev() {
            f = 2.0 * (v0 + f64::from(k) + 1.0) / x * f1 - f2;
            if k <= n {
                vl[k as usize] = f;
            }
            f2 = f1;
            f1 = f;
        }

        let cs = if bjv0.abs() > bjv1.abs() {
            bjv0 / f
        } else {
            bjv1 / f2
        };

        for k in 0..=n {
            vl[k as usize] *= cs;
        }
    }

    vl[0] *= fac;
    for j in 1..=n {
        vl[j as usize] *= fac * r0;
        dl[(j - 1) as usize] = -0.5 * x / (f64::from(j) + v0) * vl[j as usize];
        r0 *= 2.0 * (f64::from(j) + v0 + 1.0) / x;
    }

    dl[n as usize] = 2.0 * (v0 + f64::from(n)) * (vl[(n - 1) as usize] - vl[n as usize]) / x;
    *vm = f64::from(n) + v0;
}

pub(crate) fn msta1(x: f64, mp: i32) -> i32 {
    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that the magnitude of
    //          Jn(x) at that point is about 10^(-MP)
    // Input :  x     --- Argument of Jn(x)
    //          MP    --- Value of magnitude
    // Output:  MSTA1 --- Starting point
    // ===================================================

    // NOTE: the `(int)` casts below are UB in C++ for |x| >= 2^31 / 1.1 (~1.95e9) or NaN, and if
    // the secant step is NaN, +-inf (e.g. for x == 0) or huge. Rust's `as i32` saturates
    // (NaN -> 0), whereas g++ on x86-64 yields INT_MIN. The `+ 1` / `+ 5` / `nn - n1` below can
    // then overflow i32 (panics in debug builds). jynbh only calls this for x >= 1e-100 with
    // n > 0.9 * x or x <= 300.
    let a0 = x.abs();
    let mut n0 = (1.1 * a0) as i32 + 1;
    let mut f0 = 0.5 * (6.28 * f64::from(n0)).log10()
        - f64::from(n0) * (1.36 * a0 / f64::from(n0)).log10()
        - f64::from(mp);
    let mut n1 = n0 + 5;
    let mut f1 = 0.5 * (6.28 * f64::from(n1)).log10()
        - f64::from(n1) * (1.36 * a0 / f64::from(n1)).log10()
        - f64::from(mp);
    let mut nn = 0;
    for _it in 1..=20 {
        nn = (f64::from(n1) - f64::from(n1 - n0) / (1.0 - f0 / f1)) as i32;
        let f = 0.5 * (6.28 * f64::from(nn)).log10()
            - f64::from(nn) * (1.36 * a0 / f64::from(nn)).log10()
            - f64::from(mp);
        if (nn - n1).abs() < 1 {
            break;
        }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    nn
}

pub(crate) fn msta2(x: f64, n: i32, mp: i32) -> i32 {
    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that all Jn(x) has MP
    //          significant digits
    // Input :  x  --- Argument of Jn(x)
    //          n  --- Order of Jn(x)
    //          MP --- Significant digit
    // Output:  MSTA2 --- Starting point
    // ===================================================

    // NOTE: the `(int)` casts below have the same UB as in msta1 (see there), and `nn + 10` can
    // overflow as well.
    let mut n0: i32;
    let obj: f64;

    let a0 = x.abs();
    let hmp = 0.5 * f64::from(mp);
    let ejn =
        0.5 * (6.28 * f64::from(n)).log10() - f64::from(n) * (1.36 * a0 / f64::from(n)).log10();
    if ejn <= hmp {
        obj = f64::from(mp);
        n0 = (1.1 * a0) as i32 + 1;
    } else {
        obj = hmp + ejn;
        n0 = n;
    }
    let mut f0 = 0.5 * (6.28 * f64::from(n0)).log10()
        - f64::from(n0) * (1.36 * a0 / f64::from(n0)).log10()
        - obj;
    let mut n1 = n0 + 5;
    let mut f1 = 0.5 * (6.28 * f64::from(n1)).log10()
        - f64::from(n1) * (1.36 * a0 / f64::from(n1)).log10()
        - obj;
    let mut nn = 0;
    for _it in 1..=20 {
        nn = (f64::from(n1) - f64::from(n1 - n0) / (1.0 - f0 / f1)) as i32;
        let f = 0.5 * (6.28 * f64::from(nn)).log10()
            - f64::from(nn) * (1.36 * a0 / f64::from(nn)).log10()
            - obj;
        if (nn - n1).abs() < 1 {
            break;
        }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    nn + 10
}

pub(crate) fn sphj(x: f64, n: i32, nm: &mut i32, sj: &mut [f64], dj: &mut [f64]) {
    //  MODIFIED to ALLOW N=0 CASE (ALSO IN SPHY)
    //
    //  =======================================================
    //  Purpose: Compute spherical Bessel functions jn(x) and
    //           their derivatives
    //  Input :  x --- Argument of jn(x)
    //           n --- Order of jn(x)  ( n = 0,1,… )
    //  Output:  SJ(n) --- jn(x)
    //           DJ(n) --- jn'(x)
    //           NM --- Highest order computed
    //  Routines called:
    //           MSTA1 and MSTA2 for computing the starting
    //           point for backward recurrence
    //  =======================================================

    // NOTE: `sj[0..=*nm]` / `dj[0..=*nm]` are written, with *nm = min(n, msta1(x, 200)) for
    // |x| >= 1e-100. The only caller, rmn1, uses 252-element buffers, so this overflows if both n
    // and msta1(x, 200) exceed 251 (c * x >~ 30), e.g. prolate_radial1_nocv(2, 200, 1, 300)
    // (confirmed with ASan).
    // NOTE: for x NaN, msta1/msta2 hit their `(int)NaN` UB (see there), so `*nm` (and therefore
    // which elements of `sj` / `dj` are written) may differ from C++; the values are all NaN.
    *nm = n;
    if x.abs() < 1e-100 {
        for k in 0..=n {
            sj[k as usize] = 0.0;
            dj[k as usize] = 0.0;
        }
        sj[0] = 1.0;
        if n > 0 {
            dj[1] = 1.0 / 3.0;
        }
        return;
    }
    sj[0] = x.sin() / x;
    dj[0] = (x.cos() - x.sin() / x) / x;
    if n < 1 {
        return;
    }
    sj[1] = (sj[0] - x.cos()) / x;
    if n >= 2 {
        let sa = sj[0];
        let sb = sj[1];
        let mut m = msta1(x, 200);
        if m < n {
            *nm = m;
        } else {
            m = msta2(x, n, 15);
        }
        let mut f = 0.0;
        let mut f0 = 0.0;
        let mut f1 = 1e-100;
        for k in (0..=m).rev() {
            f = (2.0 * f64::from(k) + 3.0) * f1 / x - f0;
            if k <= *nm {
                sj[k as usize] = f;
            }
            f0 = f1;
            f1 = f;
        }
        let cs = if sa.abs() > sb.abs() { sa / f } else { sb / f0 };
        for k in 0..=*nm {
            sj[k as usize] *= cs;
        }
    }
    for k in 1..=*nm {
        dj[k as usize] = sj[(k - 1) as usize] - (f64::from(k) + 1.0) * sj[k as usize] / x;
    }
}

pub(crate) fn sphy(x: f64, n: i32, nm: &mut i32, sy: &mut [f64], dy: &mut [f64]) {
    // ======================================================
    // Purpose: Compute spherical Bessel functions yn(x) and
    //          their derivatives
    // Input :  x --- Argument of yn(x) ( x ≥ 0 )
    //          n --- Order of yn(x) ( n = 0,1,… )
    // Output:  SY(n) --- yn(x)
    //          DY(n) --- yn'(x)
    //          NM --- Highest order computed
    // ======================================================

    // NOTE: `sy[0..=n]` / `dy[0..=n]` are written, but the loop stops after the first
    // |sy[k]| >= 1e300. The only caller, rmn2l, uses 252-element buffers, so this overflows if
    // n > 251 and sy[251] is finite (c * x >~ 11.75), e.g. prolate_radial2_nocv(2, 200, 1, 300)
    // (confirmed with ASan).
    //
    // NOTE: unlike the Fortran original, the C++ computes no dy[k], k >= 1, if |f| >= 1e300, and
    // otherwise sets *nm = n - 1, so dy[n] is never computed. Reproduced as is.
    if x < 1.0e-60 {
        for k in 0..=n {
            sy[k as usize] = -SPECFUN_INF;
            dy[k as usize] = SPECFUN_INF;
        }
        *nm = n;
        return;
    }
    sy[0] = -x.cos() / x;
    let mut f0 = sy[0];
    dy[0] = (x.sin() + x.cos() / x) / x;

    if n < 1 {
        *nm = n;
        return;
    }

    sy[1] = (sy[0] - x.sin()) / x;
    let mut f1 = sy[1];

    for k in 2..=n {
        let f = ((2.0 * f64::from(k) - 1.0) * f1 / x) - f0;
        sy[k as usize] = f;
        if f.abs() >= SPECFUN_INF {
            *nm = k - 1;
            return;
        }
        f0 = f1;
        f1 = f;
    }
    *nm = n - 1;
    for k in 1..=*nm {
        dy[k as usize] = sy[(k - 1) as usize] - (f64::from(k) + 1.0) * sy[k as usize] / x;
    }
}
