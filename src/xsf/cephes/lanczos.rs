//! The subset of `xsf/cephes/lanczos.h` (xsf v0.2.2) that is used in Rust, translated into pure
//! Rust.
//!
//! (C) Copyright John Maddock 2006.
//! Use, modification and distribution are subject to the
//! Boost Software License, Version 1.0. (See accompanying file
//! LICENSE_1_0.txt or copy at <https://www.boost.org/LICENSE_1_0.txt>)
//!
//! Both lanczos.h and lanczos.c were formed from Boost's lanczos.hpp
//!
//! Scipy changes:
//! - 06-22-2016: Removed all code not related to double precision and
//!   ported to c for use in Cephes. Note that the order of the
//!   coefficients is reversed to match the behavior of polevl.
//!
//! Optimal values for G for each N are taken from
//! <https://web.viu.ca/pughg/phdThesis/phdThesis.pdf>,
//! as are the theoretical error bounds.
//!
//! Constants calculated using the method described by Godfrey
//! <https://my.fit.edu/~gabdo/gamma.txt> and elaborated by Toth at
//! <https://www.rskey.org/gamma.htm> using NTL::RR at 1000 bit precision.
//!
//! Lanczos Coefficients for N=13 G=6.024680040776729583740234375
//! Max experimental error (with arbitrary precision arithmetic) 2.852e-17 for 0 < x <= 172
//! Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//!
//! Use for double precision.

use crate::xsf::cephes::polevl::ratevl;

const LANCZOS_SUM_EXPG_SCALED_NUM: [f64; 13] = [
    6.061_842_346_248_907e-3,
    5.098_416_655_656_676e-1,
    19.519_927_882_476_175,
    449.944_556_906_316_8,
    6_955.999_602_515_376,
    75_999.293_040_145_42,
    601_859.617_168_109_9,
    3_481_712.154_980_646,
    14_605_578.087_685_067,
    43_338_889.324_676_14,
    86_363_131.288_138_6,
    103_794_043.116_344_54,
    56_906_521.913_471_565,
];

const LANCZOS_SUM_EXPG_SCALED_DENOM: [f64; 13] = [
    1.0,
    66.0,
    1_925.0,
    32_670.0,
    357_423.0,
    2_637_558.0,
    13_339_535.0,
    45_995_730.0,
    105_258_076.0,
    150_917_976.0,
    120_543_840.0,
    39_916_800.0,
    0.0,
];

#[doc(hidden)]
#[must_use]
#[inline]
pub fn lanczos_sum_expg_scaled(x: f64) -> f64 {
    ratevl(
        x,
        &LANCZOS_SUM_EXPG_SCALED_NUM,
        &LANCZOS_SUM_EXPG_SCALED_DENOM,
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_lanczos_sum_expg_scaled() {
        xsref::test("lanczos_sum_expg_scaled", "d-d", |x| {
            crate::lanczos_sum_expg_scaled(x[0])
        });
    }
}
