//! Translated into pure Rust from `xsf/cephes/zeta.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the asymptotic expansion for `q > 1e8` includes the `B_2` term, which reduces
//! its relative error from `x (x - 1) / (12 q^2)` (up to ~1e-14) to below machine precision.

// cephes/const.h
const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53

/// Expansion coefficients
/// for Euler-Maclaurin summation formula
/// (2k)! / B2k
/// where B2k are Bernoulli numbers
const ZETA_A: [f64; 12] = [
    12.0,
    -720.0,
    30_240.0,
    -1_209_600.0,
    47_900_160.0,
    -1.892_437_580_318_379_2e9, /*1.307674368e12/691 */
    7.472_424_96e10,
    -2.950_130_727_918_164e12,  /*1.067062284288e16/3617 */
    1.164_678_281_435_006_7e14, /*5.109094217170944e18/43867 */
    -4.597_978_722_407_473e15,  /*8.028576626982912e20/174611 */
    1.815_210_540_194_354_6e17, /*1.5511210043330985984e23/854513 */
    -7.166_165_256_175_667e18,  /*1.6938241367317436694528e27/236364091 */
];

/* 30 Nov 86 -- error in third coefficient fixed */

/// Riemann zeta function of two arguments
///
/// ```text
///                 inf.
///                  -        -x
///   zeta(x,q)  =   >   (k+q)
///                  -
///                 k=0
/// ```
///
/// where x > 1 and q is not a negative integer or zero.
/// The Euler-Maclaurin summation formula is used to obtain
/// the expansion
///
/// ```text
///                n
///                -       -x
/// zeta(x,q)  =   >  (k+q)
///                -
///               k=1
///
///           1-x                 inf.  B   x(x+1)...(x+2j)
///      (n+q)           1         -     2j
///  +  ---------  -  -------  +   >    --------------------
///        x-1              x      -                   x+2j+1
///                   2(n+q)      j=1       (2j)! (n+q)
/// ```
///
/// where the B2j are Bernoulli numbers.  Note that (see zetac.c)
/// zeta(x,1) = zetac(x) + 1.
///
/// REFERENCE:
///
/// Gradshteyn, I. S., and I. M. Ryzhik, Tables of Integrals,
/// Series, and Products, p. 1073; Academic Press, 1980.
#[allow(clippy::float_cmp)]
pub(crate) fn zeta(x: f64, q: f64) -> f64 {
    if x == 1.0 {
        // retinf
        return f64::INFINITY;
    }

    if x < 1.0 {
        // domerr
        return f64::NAN;
    }

    if q <= 0.0 {
        if q == q.floor() {
            // retinf
            return f64::INFINITY;
        }
        if x != x.floor() {
            // domerr
            return f64::NAN; /* because q^-x not defined */
        }
    }

    /* Asymptotic expansion
     * https://dlmf.nist.gov/25.11#E43
     */
    if q > 1e8 {
        // (unlike in xsf) the B_2 term, which is skipped for x = inf to avoid inf * 0
        let b2 = if x.is_finite() {
            x / (12.0 * q * q)
        } else {
            0.0
        };
        return (1.0 / (x - 1.0) + 1.0 / (2.0 * q) + b2) * q.powf(1.0 - x);
    }

    /* Euler-Maclaurin summation formula */

    /* Permit negative q but continue sum until n+q > +9 .
     * This case should be handled by a reflection formula.
     * If q<0 and x is an integer, there is a relation to
     * the polyGamma function.
     */
    let mut s = q.powf(-x);
    let mut a = q;
    let mut i = 0;
    let mut b = 0.0;
    while (i < 9) || (a <= 9.0) {
        i += 1;
        a += 1.0;
        b = a.powf(-x);
        s += b;
        // Unlike in xsf, this also requires the tail sum `~ b * a / (x - 1)` to be negligible,
        // which it isn't for `x` near 1 and tiny `q`, e.g. `zeta(1 + 2^-52, 1e-15)` was 1.0e15
        // instead of 5.5e15.
        if (b / s).abs() < MACHEP && (b * a / ((x - 1.0) * s)).abs() < MACHEP {
            return s;
        }
    }

    let w = a;
    s += b * w / (x - 1.0);
    s -= 0.5 * b;
    a = 1.0;
    let mut k = 0.0;
    for zeta_a in ZETA_A {
        a *= x + k;
        b /= w;
        let mut t = a * b / zeta_a;
        s += t;
        t = (t / s).abs();
        if t < MACHEP {
            return s;
        }
        k += 1.0;
        a *= x + k;
        b /= w;
        k += 1.0;
    }
    s
}
