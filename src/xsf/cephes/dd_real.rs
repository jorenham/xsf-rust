//! The subset of `xsf/cephes/dd_real.h` (xsf v0.2.2) that is needed by the complex `log1p`,
//! `smirnov`, and `iv_ratio_c`.
//!
//! Translated into C++ by SciPy developers in 2024. The parts of the qd double-double floating
//! point package used in SciPy have been reworked in a more modern C++ style using operator
//! overloading. Original header:
//!
//! > This work was supported by the Director, Office of Science, Division of Mathematical,
//! > Information, and Computational Sciences of the U.S. Department of Energy under contract
//! > numbers DE-AC03-76SF00098 and DE-AC02-05CH11231.
//! >
//! > Copyright (c) 2003-2009, The Regents of the University of California, through Lawrence
//! > Berkeley National Laboratory (subject to receipt of any required approvals from U.S. Dept.
//! > of Energy) All rights reserved.
//! >
//! > By downloading or using this software you are agreeing to the modified BSD license
//! > "BSD-LBNL-License.doc" (see LICENSE.txt).
//!
//! Double-double precision (>= 106-bit significand) floating point arithmetic package based on
//! David Bailey's Fortran-90 double-double package, with some changes, by Yozo Hida. This code was
//! taken from v2.3.18 of the qd package.
//!
//! The C++ code uses `volatile` to prevent compilers from reassociating these error-free
//! transformations; Rust never reassociates or contracts floating point operations, so that's not
//! needed here. Note that the error terms of the products are only exact if they don't underflow.

use core::cmp::Ordering;
use core::f64::consts::{E as E_HI, LN_2};
use core::ops::{Add, Div, Mul, Neg, Sub};

/*************************************************************************
 * The basic routines taking double arguments, returning 1 (or 2) doubles
 *************************************************************************/

/// Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|.
#[inline]
fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let c = s - a;
    (s, b - c)
}

/// Computes fl(a+b) and err(a+b).
#[inline]
pub(crate) fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let c = s - a;
    let d = b - c;
    let e = s - c;
    (s, (a - e) + d)
}

/// Computes fl(a*b) and err(a*b).
#[inline]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

/// Computes fl(a*a) and err(a*a).  Faster than the above method.
#[inline]
fn two_sqr(a: f64) -> (f64, f64) {
    let p = a * a;
    (p, a.mul_add(a, -p))
}

/// `std::ldexp`, i.e. `x * 2^n` without intermediate overflow or double rounding (as in musl).
#[allow(clippy::cast_sign_loss)]
pub(crate) fn ldexp_f64(x: f64, mut n: i32) -> f64 {
    const P1023: f64 = f64::from_bits(0x7fe0_0000_0000_0000); // 2^1023
    const M969: f64 = f64::from_bits(0x0360_0000_0000_0000); // 2^-1022 * 2^53

    let mut y = x;
    if n > 1023 {
        y *= P1023;
        n -= 1023;
        if n > 1023 {
            y *= P1023;
            n -= 1023;
            n = n.min(1023);
        }
    } else if n < -1022 {
        // make sure the final n < -53 to avoid double rounding in the subnormal range
        y *= M969;
        n += 1022 - 53;
        if n < -1022 {
            y *= M969;
            n += 1022 - 53;
            n = n.max(-1022);
        }
    }
    y * f64::from_bits(((0x3ff + n) as u64) << 52)
}

/// `std::frexp`, i.e. `(m, e)` such that `x = m * 2^e` with `0.5 <= |m| < 1` (as in glibc).
#[allow(clippy::cast_possible_truncation)]
pub(crate) fn frexp_f64(x: f64) -> (f64, i32) {
    let bits = x.to_bits();
    let ee = ((bits >> 52) & 0x7ff) as i32;
    if ee == 0 {
        if x == 0.0 {
            return (x, 0);
        }
        // subnormal
        let (m, e) = frexp_f64(x * f64::from_bits(0x43f0_0000_0000_0000)); // 2^64
        return (m, e - 64);
    }
    if ee == 0x7ff {
        // inf or nan
        return (x, 0);
    }
    let m = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | 0x3fe0_0000_0000_0000);
    (m, ee - 0x3fe)
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct DoubleDouble {
    pub(crate) hi: f64,
    pub(crate) lo: f64,
}

impl DoubleDouble {
    #[inline]
    pub(crate) const fn new(high: f64) -> Self {
        Self { hi: high, lo: 0.0 }
    }
}

impl From<DoubleDouble> for f64 {
    #[inline]
    fn from(x: DoubleDouble) -> Self {
        x.hi
    }
}

// Arithmetic operations

impl Neg for DoubleDouble {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

impl Add for DoubleDouble {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        /* This one satisfies IEEE style error bound,
        due to K. Briggs and W. Kahan.                   */
        let (mut s1, mut s2) = two_sum(self.hi, rhs.hi);
        let (t1, t2) = two_sum(self.lo, rhs.lo);
        s2 += t1;
        (s1, s2) = quick_two_sum(s1, s2);
        s2 += t2;
        (s1, s2) = quick_two_sum(s1, s2);
        Self { hi: s1, lo: s2 }
    }
}

impl Add<f64> for DoubleDouble {
    type Output = Self;

    #[inline]
    fn add(self, rhs: f64) -> Self {
        let (mut s1, mut s2) = two_sum(self.hi, rhs);
        s2 += self.lo;
        (s1, s2) = quick_two_sum(s1, s2);
        Self { hi: s1, lo: s2 }
    }
}

impl Sub for DoubleDouble {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl Sub<f64> for DoubleDouble {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: f64) -> Self {
        let (mut s1, mut s2) = two_sum(self.hi, -rhs);
        s2 += self.lo;
        (s1, s2) = quick_two_sum(s1, s2);
        Self { hi: s1, lo: s2 }
    }
}

impl Sub<DoubleDouble> for f64 {
    type Output = DoubleDouble;

    #[inline]
    fn sub(self, rhs: DoubleDouble) -> DoubleDouble {
        let (mut s1, mut s2) = two_sum(self, -rhs.hi);
        s2 -= rhs.lo;
        (s1, s2) = quick_two_sum(s1, s2);
        DoubleDouble { hi: s1, lo: s2 }
    }
}

impl Mul for DoubleDouble {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let (mut p1, mut p2) = two_prod(self.hi, rhs.hi);
        p2 += self.hi * rhs.lo + self.lo * rhs.hi;
        (p1, p2) = quick_two_sum(p1, p2);
        Self { hi: p1, lo: p2 }
    }
}

impl Mul<f64> for DoubleDouble {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: f64) -> Self {
        let (p1, e1) = two_prod(self.hi, rhs);
        let (p2, e2) = two_prod(self.lo, rhs);
        let (p1, e1) = quick_two_sum(p1, e2 + p2 + e1);
        Self { hi: p1, lo: e1 }
    }
}

impl Mul<DoubleDouble> for f64 {
    type Output = DoubleDouble;

    #[inline]
    fn mul(self, rhs: DoubleDouble) -> DoubleDouble {
        let (p1, e1) = two_prod(self, rhs.hi);
        let (p2, e2) = two_prod(self, rhs.lo);
        let (p1, e1) = quick_two_sum(p1, e2 + p2 + e1);
        DoubleDouble { hi: p1, lo: e1 }
    }
}

impl Div for DoubleDouble {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self {
        let mut q1 = self.hi / rhs.hi; /* approximate quotient */

        let mut r = self - rhs * q1;

        let mut q2 = r.hi / rhs.hi;
        r = r - rhs * q2;

        let q3 = r.hi / rhs.hi;

        (q1, q2) = quick_two_sum(q1, q2);
        Self { hi: q1, lo: q2 } + q3
    }
}

impl Div<f64> for DoubleDouble {
    type Output = Self;

    #[inline]
    fn div(self, rhs: f64) -> Self {
        self / Self::new(rhs)
    }
}

impl Div<DoubleDouble> for f64 {
    type Output = DoubleDouble;

    #[inline]
    fn div(self, rhs: DoubleDouble) -> DoubleDouble {
        DoubleDouble::new(self) / rhs
    }
}

impl PartialEq for DoubleDouble {
    #[inline]
    fn eq(&self, rhs: &Self) -> bool {
        self.hi == rhs.hi && self.lo == rhs.lo
    }
}

impl PartialEq<f64> for DoubleDouble {
    #[inline]
    fn eq(&self, rhs: &f64) -> bool {
        self.hi == *rhs && self.lo == 0.0
    }
}

// Like the C++ `operator<` and `operator<=`, the low words are compared if the high words are
// neither `<` nor `>`, i.e. also if one of them is NaN.

impl PartialOrd for DoubleDouble {
    #[inline]
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        if self.hi < rhs.hi {
            Some(Ordering::Less)
        } else if self.hi > rhs.hi {
            Some(Ordering::Greater)
        } else {
            self.lo.partial_cmp(&rhs.lo)
        }
    }
}

impl PartialOrd<f64> for DoubleDouble {
    #[inline]
    fn partial_cmp(&self, rhs: &f64) -> Option<Ordering> {
        if self.hi < *rhs {
            Some(Ordering::Less)
        } else if self.hi > *rhs {
            Some(Ordering::Greater)
        } else {
            self.lo.partial_cmp(&0.0)
        }
    }
}

// Math functions

impl DoubleDouble {
    #[inline]
    fn mul_pwr2(self, rhs: f64) -> Self {
        /* double-double * double,  where double is a power of 2. */
        Self {
            hi: self.hi * rhs,
            lo: self.lo * rhs,
        }
    }

    #[inline]
    pub(crate) fn is_finite(self) -> bool {
        self.hi.is_finite()
    }

    #[allow(clippy::float_cmp)]
    pub(crate) fn floor(self) -> Self {
        let mut hi = self.hi.floor();
        let mut lo = 0.0;

        if hi == self.hi {
            /* High word is integer already.  Round the low word. */
            lo = self.lo.floor();
            (hi, lo) = quick_two_sum(hi, lo);
        }

        Self { hi, lo }
    }

    #[inline]
    pub(crate) fn abs(self) -> Self {
        if self.hi < 0.0 { -self } else { self }
    }

    #[inline]
    pub(crate) fn fma(self, b: Self, c: Self) -> Self {
        // TODO: make an accurate fma
        self * b + c
    }

    #[inline]
    fn square(self) -> Self {
        let (p1, mut p2) = two_sqr(self.hi);
        p2 += 2.0 * self.hi * self.lo;
        p2 += self.lo * self.lo;
        let (s1, s2) = quick_two_sum(p1, p2);
        Self { hi: s1, lo: s2 }
    }

    #[inline]
    pub(crate) fn ldexp(self, expt: i32) -> Self {
        // float128 * (2.0 ^ expt)
        Self {
            hi: ldexp_f64(self.hi, expt),
            lo: ldexp_f64(self.lo, expt),
        }
    }

    #[allow(clippy::float_cmp)]
    pub(crate) fn frexp(self) -> (Self, i32) {
        //    r"""return b and l s.t. 0.5<=|b|<1 and 2^l == a
        //    0.5<=|b[0]|<1.0 or |b[0]| == 1.0 and b[0]*b[1]<0
        //    """
        let (mut man, mut exponent) = frexp_f64(self.hi);
        let mut b1 = ldexp_f64(self.lo, -exponent);
        if man.abs() == 0.5 && man * b1 < 0.0 {
            man *= 2.0;
            b1 *= 2.0;
            exponent -= 1;
        }
        (Self { hi: man, lo: b1 }, exponent)
    }

    // Numeric limits

    pub(crate) const NAN: Self = Self {
        hi: f64::NAN,
        lo: f64::NAN,
    };

    pub(crate) const INFINITY: Self = Self {
        hi: f64::INFINITY,
        lo: f64::INFINITY,
    };
}

#[rustfmt::skip]
const INV_FACT: [DoubleDouble; 15] = [
    DoubleDouble { hi: 1.666_666_666_666_666_6e-1, lo: 9.251_858_538_542_97e-18 },
    DoubleDouble { hi: 4.166_666_666_666_666_4e-2, lo: 2.312_964_634_635_742_7e-18 },
    DoubleDouble { hi: 8.333_333_333_333_333e-3, lo: 1.156_482_317_317_871_4e-19 },
    DoubleDouble { hi: 1.388_888_888_888_889e-3, lo: -5.300_543_954_373_577e-20 },
    DoubleDouble { hi: 1.984_126_984_126_984e-4, lo: 1.720_955_829_342_070_5e-22 },
    DoubleDouble { hi: 2.480_158_730_158_73e-5, lo: 2.151_194_786_677_588_2e-23 },
    DoubleDouble { hi: 2.755_731_922_398_589_3e-6, lo: -1.858_393_274_046_472e-22 },
    DoubleDouble { hi: 2.755_731_922_398_589e-7, lo: 2.376_771_462_225_029_7e-23 },
    DoubleDouble { hi: 2.505_210_838_544_172e-8, lo: -1.448_814_070_935_912e-24 },
    DoubleDouble { hi: 2.087_675_698_786_81e-9, lo: -1.207_345_059_113_26e-25 },
    DoubleDouble { hi: 1.605_904_383_682_161_3e-10, lo: 1.258_529_458_875_209_8e-26 },
    DoubleDouble { hi: 1.147_074_559_772_972_5e-11, lo: 2.065_551_275_283_074_5e-28 },
    DoubleDouble { hi: 7.647_163_731_819_816e-13, lo: 7.038_728_777_334_53e-30 },
    DoubleDouble { hi: 4.779_477_332_387_385e-14, lo: 4.399_205_485_834_081e-31 },
    DoubleDouble { hi: 2.811_457_254_345_520_6e-15, lo: 1.650_884_273_086_143_3e-31 },
];

// Math constants
const E: DoubleDouble = DoubleDouble {
    hi: E_HI,
    lo: 1.445_646_891_729_250_2e-16,
};
const LOG2: DoubleDouble = DoubleDouble {
    hi: LN_2,
    lo: 2.319_046_813_846_299_6e-17,
};
const EPS: f64 = 4.930_380_657_631_32e-32; // 2^-104

impl DoubleDouble {
    /// Exponential.  Computes exp(x) in double-double precision.
    #[allow(clippy::float_cmp, clippy::cast_possible_truncation)]
    pub(crate) fn exp(self) -> Self {
        /* Strategy:  We first reduce the size of x by noting that

        exp(kr + m * log(2)) = 2^m * exp(r)^k

        where m and k are integers.  By choosing m appropriately
        we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
        evaluated using the familiar Taylor series.  Reducing the
        argument substantially speeds up the convergence.       */

        const K: f64 = 512.0;
        const INV_K: f64 = 1.0 / K;
        let mut i = 0;

        if self.hi <= -709.0 {
            return Self::new(0.0);
        }

        if self.hi >= 709.0 {
            return Self::INFINITY;
        }

        if self == 0.0 {
            return Self::new(1.0);
        }

        if self == 1.0 {
            return E;
        }

        let m = (self.hi / LOG2.hi + 0.5).floor();
        let r = (self - LOG2 * m).mul_pwr2(INV_K);

        let mut p = r.square();
        let mut s = r + p.mul_pwr2(0.5);
        p = p * r;
        let mut t = p * INV_FACT[0];
        loop {
            s = s + t;
            p = p * r;
            i += 1;
            t = p * INV_FACT[i];
            if !(f64::from(t).abs() > INV_K * EPS && i < 5) {
                break;
            }
        }

        s = s + t;

        for _ in 0..9 {
            s = s.mul_pwr2(2.0) + s.square();
        }
        s = s + 1.0;

        s.ldexp(m as i32)
    }

    /// Logarithm.  Computes log(x) in double-double precision.
    /// This is a natural logarithm (i.e., base e).
    #[allow(clippy::float_cmp)]
    pub(crate) fn log(self) -> Self {
        /* Strategy.  The Taylor series for log converges much more
        slowly than that of exp, due to the lack of the factorial
        term in the denominator.  Hence this routine instead tries
        to determine the root of the function

        f(x) = exp(x) - a

        using Newton iteration.  The iteration is given by

        x' = x - f(x)/f'(x)
        = x - (1 - a * exp(-x))
        = x + a * exp(-x) - 1.

        Only one iteration is needed, since Newton's iteration
        approximately doubles the number of digits per iteration. */

        if self == 1.0 {
            return Self::new(0.0);
        }

        if self.hi <= 0.0 {
            return Self::NAN;
        }

        let mut x = Self::new(self.hi.ln()); /* Initial approximation */

        /* x = x + a * exp(-x) - 1.0; */
        x = x + self * (-x).exp() - 1.0;
        x
    }

    pub(crate) fn log1p(self) -> Self {
        if self.hi <= -1.0 {
            return -Self::INFINITY;
        }
        let la = self.hi.ln_1p();
        let elam1 = la.exp_m1();
        let mut ll = (self.lo / (1.0 + self.hi)).ln_1p();
        if self.hi > 0.0 {
            ll -= (elam1 - self.hi) / (elam1 + 1.0);
        }
        Self::new(la) + ll
    }
}

#[cfg(test)]
mod tests {
    use super::DoubleDouble;

    #[test]
    fn test_cmp_nan_hi() {
        // like in C++, the low words are compared if the high words are unordered
        let nan = DoubleDouble {
            hi: f64::NAN,
            lo: 0.0,
        };
        let zero = DoubleDouble::new(0.0);
        assert!(nan <= zero);
        assert!(nan <= 0.0);
    }
}
