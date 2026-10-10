//! The subset of `xsf/cephes/dd_real.h` (xsf v0.2.2) that is needed by the complex `log1p`.
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

use core::ops::{Add, Mul};

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
fn two_sum(a: f64, b: f64) -> (f64, f64) {
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
