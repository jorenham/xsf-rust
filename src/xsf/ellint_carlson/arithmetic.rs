//! Translated into pure Rust from `xsf/ellint_carlson/ellint_arithmetic.hh` (xsf main).
//!
//! ```text
//! Algorithms for accurate floating-point summation, dot-product, and
//! polynomial evaluation.
//!
//! Ref:
//!
//! [1] S. Gaillat, V. Ménissier-Morain: Compensated Horner scheme in complex
//!      floating point arithmetic.
//!      2008-07-07, Proc. 8th Conf. Real Number Comput., pp. 133--146
//!      https://www-pequan.lip6.fr/~graillat/papers/rnc08.pdf
//!
//! [2] S. Gaillat, V. Ménissier-Morain: Accurate summation, dot product and
//!      polynomial evaluation in complex floating point arithmetic.
//!      2012-03-30, Inf. Comput., vol. 216 pp. 57--71
//!      https://doi.org/10.1016/j.ic.2011.09.003
//!      https://web.stanford.edu/group/SOL/software/qdotdd/IC2012.pdf
//!
//! [3] S. Graillat, Ph.Langlois, N.Louve: Compensated Horner Scheme.
//!      2005-07-24, Report No. RR2005-04, Université de Perpignan « Via
//!      Domitia »
//!      https://www-pequan.lip6.fr/~jmc/polycopies/Compensation-horner.pdf
//! ```
//!
//! The output arguments that the C++ functions take by (non-const) reference are returned instead,
//! except for the accumulators.

use num_complex::Complex64;

use super::arith_aux as aux;
use super::common::NPOLY;
use super::typing::Num;
use crate::xsf::complex::cmul;

/* Knuth's TwoSum algorithm, EFT on addition, with the sum written to s
 * and the correction term to p. This algorithm is generic enough to
 * write for both real and complex types (because there's no error in the
 * complex "i").
 * Algorithms 2.1 and 3.1 in Ref. [1] */
#[inline]
pub(super) fn eft_sum<T: Num>(x: T, y: T) -> (T, T) {
    let s = x + y;
    let z = s - x;
    let corr = (x - (s - z)) + (y - z);
    (s, corr)
}

/* TwoSum in accumulator style, with sum accumulated to acc and the
 * correction term to corr
 * Algorithm 4.1 in Ref. [2] */
#[inline]
pub(super) fn sum2_acc<T: Num>(summand: T, acc: &mut T, corr: &mut T) {
    let (tmp_s, tmp_e) = eft_sum(summand, *acc);
    *acc = tmp_s;
    *corr = *corr + tmp_e;
}

/* Sum at most the first n elements of iterable thing, including fixed-size
 * array passed by reference (both real and complex). */
#[inline]
pub(super) fn nsum2<T: Num>(x: &[T], n: usize) -> T {
    let mut p = T::from(0.0);
    let mut s = T::from(0.0);

    for &xi in x.iter().take(n) {
        sum2_acc(xi, &mut p, &mut s);
    }

    p + s
}

/* Convenient wrapper for summing fixed-size array. */
#[inline]
pub(super) fn sum2<T: Num>(arr: &[T]) -> T {
    nsum2(arr, arr.len())
}

/* EFT for multiplication, with real-type arguments.
 * NOTE: This implementation uses the FMA function for correct rounding.
 * (mandated by standard since C++11)
 * Algorithm 2.5 in Ref. [1] */
#[inline]
pub(super) fn eft_prod(x: f64, y: f64) -> (f64, f64) {
    let prod = x * y;
    let corr = x.mul_add(y, -prod);
    (prod, corr)
}

/* Buffer holding the correction terms for multiplication.  For real
 * floating-point types, this is just the same type ("just" numeric type).
 * For complex types, there must be two buffers, each with length 4, for
 * holding the EFT terms. */
/// `typing::corrbuf<std::complex<double>>`
#[derive(Default)]
struct CorrBuf {
    value_real: [f64; NPOLY],
    value_imag: [f64; NPOLY],
}

/* EFT for multiplication, with complex-type arguments.
 * Algorithm 3.4 in Ref. [1] */
#[inline]
fn eft_prod_cplx(x: Complex64, y: Complex64, corr: &mut CorrBuf) -> Complex64 {
    /*  p                  e    f    g
     * z5 | value_real : [h1, -h2,  h5, xx]
     * z6 | value_imag : [h3,  h4,  h6, xx]
     */
    let a = x.re;
    let b = x.im;
    let c = y.re;
    let d = y.im;

    let z1;
    let z2;
    let z3;
    let z4;
    let z5;
    let z6;
    (z1, corr.value_real[0]) = eft_prod(a, c);
    (z2, corr.value_real[1]) = eft_prod(b, d);
    (z3, corr.value_imag[0]) = eft_prod(a, d);
    (z4, corr.value_imag[1]) = eft_prod(b, c);
    (z5, corr.value_real[2]) = eft_sum(z1, -z2);
    (z6, corr.value_imag[2]) = eft_sum(z3, z4);

    corr.value_real[1] = -corr.value_real[1];
    Complex64::new(z5, z6)
}

/* Accumulator version of dot-product
 * Algorithm 4.3 in Ref. [2] */
#[inline]
pub(super) fn fdot2_acc(x: f64, y: f64, acc: &mut f64, corr: &mut f64) {
    let (h, r) = eft_prod(x, y);
    let (t, q) = eft_sum(*acc, h);
    *acc = t;
    *corr += q + r;
}

/* Dot-product of at most first n elements, for real floating type only
 * (based on what we see by peeking "begin"). */
pub(super) fn ndot2_real(x: &[f64], y: &[f64], n: usize) -> f64 {
    let mut p = 0.0;
    let mut s = 0.0;

    for (&xi, &yi) in x.iter().zip(y).take(n) {
        fdot2_acc(xi, yi, &mut p, &mut s);
    }

    p + s
}

/* Dot-product, but for complex types.
 * Algorithm 4.4 in Ref. [2] */
pub(super) fn ndot2_cplx(x: &[Complex64], y: &[Complex64], n: usize) -> Complex64 {
    let mut pr = 0.0;
    let mut pi = 0.0;
    let mut cr = 0.0;
    let mut ci = 0.0;

    for (xi, yi) in x.iter().zip(y).take(n) {
        let a = xi.re;
        let b = yi.re;
        let c = xi.im;
        let d = yi.im;

        fdot2_acc(a, b, &mut pr, &mut cr);
        fdot2_acc(c, -d, &mut pr, &mut cr);
        fdot2_acc(a, d, &mut pi, &mut ci);
        fdot2_acc(b, c, &mut pi, &mut ci);
    }

    Complex64::new(pr + cr, pi + ci)
}

/// `ndot2(x, y, n)` for real or complex `T`
#[inline]
pub(super) fn ndot2<T: Num>(x: &[T], y: &[T], n: usize) -> T {
    T::ndot2(x, y, n)
}

/* Convenient wrapper for dot-product */
#[inline]
pub(super) fn dot2<T: Num>(x: &[T], y: &[T]) -> T {
    let len = x.len().min(y.len());
    ndot2(x, y, len)
}

/* Polynomial evaluation with compensated Horner scheme.
 * This instance for the real-only inputs is essentially the Algorithm 9 in
 * Ref. [3] */
pub(super) fn dcomp_horner_real(x: f64, poly: &[f64], degree: usize) -> f64 {
    let mut s = poly[degree];
    let mut r = 0.0;

    for i in (0..degree).rev() {
        let (ptmp, pp) = eft_prod(s, x);
        let ps;
        (s, ps) = eft_sum(ptmp, poly[i]);
        r = r * x + (pp + ps);
    }
    s + r
}

/* For complex-valued input: Algorithm 5.4 in Ref. [2].
 * Notice that this makes use of the aux::acc_sum() algorithm, which is
 * named "Accsum" in the referenced paper. */
pub(super) fn dcomp_horner_cplx(x: Complex64, poly: &[f64], degree: usize) -> Complex64 {
    let mut s = Complex64::from(poly[degree]);
    let mut r = Complex64::from(0.0);

    for i in (0..degree).rev() {
        let mut ws = CorrBuf::default();
        let mut mask_r = [false; NPOLY];
        let mut mask_i = [false; NPOLY];

        let p = eft_prod_cplx(s, x, &mut ws);
        let tmp; /* tmp is for copying into the buffer arrays. */
        (s, tmp) = eft_sum(p, Complex64::from(poly[i]));
        ws.value_real[NPOLY - 1] = tmp.re;
        ws.value_imag[NPOLY - 1] = tmp.im;
        aux::make(&ws.value_real, &mut mask_r);
        aux::make(&ws.value_imag, &mut mask_i);
        r = cmul(r, x)
            + Complex64::new(
                aux::acc_sum(&mut ws.value_real, &mut mask_r),
                aux::acc_sum(&mut ws.value_imag, &mut mask_i),
            );
    }

    s + r
}

/* Convenient version for fixed-size array as polynomial coefficients, with
 * type promotion. */
#[inline]
pub(super) fn comp_horner<T: Num>(x: T, poly: &[f64]) -> T {
    x.dcomp_horner(poly, poly.len() - 1)
}

/* Small integer roots for real types. */
#[inline]
pub(super) fn ocrt(x: f64) -> f64 /* 8th (octic) root */ {
    x.sqrt().sqrt().sqrt()
}
