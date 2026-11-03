//! Translated into pure Rust from `xsf/trigamma.h` (xsf main), which is based on
//! [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl):
//!
//! > The MIT License (MIT)
//! >
//! > Copyright (c) 2017 Jeff Bezanson, Stefan Karpinski, Viral B. Shah, and others:
//! >
//! > <https://github.com/JuliaMath/SpecialFunctions.jl/graphs/contributors>
//! >
//! > Permission is hereby granted, free of charge, to any person obtaining a copy
//! > of this software and associated documentation files (the "Software"), to deal
//! > in the Software without restriction, including without limitation the rights
//! > to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//! > copies of the Software, and to permit persons to whom the Software is
//! > furnished to do so, subject to the following conditions:
//! >
//! > The above copyright notice and this permission notice shall be included in all
//! > copies or substantial portions of the Software.
//! >
//! > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//! > IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//! > FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//! > AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//! > LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//! > OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//! > SOFTWARE.

use core::f64::consts::PI;

use num_complex::{Complex, c64};

use crate::xsf::cephes::polevl::polevl;
use crate::xsf::complex::{cdiv, cmul};
use crate::xsf::evalpoly::cevalpoly;
use crate::xsf::trig::cephes_sinpi;

const N: i32 = 10;

// order of the Bernoulli numbers is reversed due to `evalpoly` convention
const COEFFS: [f64; 8] = [
    -7.092_156_862_745_098,
    1.166_666_666_666_666_7,
    -0.253_113_553_113_553_1,
    0.075_757_575_757_575_76,
    -0.033_333_333_333_333_33,
    0.023_809_523_809_523_808,
    -0.033_333_333_333_333_33,
    0.166_666_666_666_666_66,
];

// https://github.com/JuliaMath/SpecialFunctions.jl/blob/c9ff36085c9008b9bdaa457e745f459a1827e52e/src/gamma.jl#L65
// Implementation based on differentiation of equation 12 of
// K.S. Kolbig, Programs for computing the logarithm of the gamma function, and the
// digamma function, for complex argument, Computer Physics Communications, 1972
// https://doi.org/10.1016/0010-4655(72)90012-4.

/// `xsf::trigamma` for real `z`
#[allow(clippy::cast_possible_truncation)]
fn trigamma_f64(mut z: f64) -> f64 {
    let x = z;
    if x <= 0.0 {
        // reflection formula http://dlmf.nist.gov/5.15.E6
        let s = cephes_sinpi(z);
        let c = PI / s;
        return c * c - trigamma_f64(1.0 - z);
    }

    let mut psi = 0.0;
    if x < f64::from(N) {
        // recursion based on https://dlmf.nist.gov/5.15.E1
        let n = (f64::from(N) - x.floor()) as i32;
        psi += 1.0 / (z * z);
        for v in 1..n {
            let w = z + f64::from(v);
            psi += 1.0 / (w * w);
        }
        z += f64::from(n);
    }

    // https://dlmf.nist.gov/5.15.E8
    let t = 1.0 / z;
    let w = t * t;
    psi += t + 0.5 * w;
    psi += t * w * polevl(w, &COEFFS);
    psi
}

/// `xsf::trigamma` for complex `z`
#[allow(clippy::cast_possible_truncation)]
fn trigamma_c64(mut z: Complex<f64>) -> Complex<f64> {
    let x = z.re;
    if x <= 0.0 {
        // reflection formula http://dlmf.nist.gov/5.15.E6
        let s = crate::sinpi(z);
        let c = cdiv(PI.into(), s);
        return cmul(c, c) - trigamma_c64(c64(1.0, 0.0) - z);
    }

    let mut psi = Complex::ZERO;
    if x < f64::from(N) {
        // recursion based on https://dlmf.nist.gov/5.15.E1
        let n = (f64::from(N) - x.floor()) as i32;
        psi += cdiv(Complex::ONE, z * z);
        for v in 1..n {
            let w = z + f64::from(v);
            psi += cdiv(Complex::ONE, w * w);
        }
        z += f64::from(n);
    }

    // https://dlmf.nist.gov/5.15.E8
    let t = cdiv(Complex::ONE, z);
    let w = t * t;
    psi += t + 0.5 * w;
    psi += t * w * cevalpoly(&COEFFS, w);
    psi
}

pub trait TrigammaArg: crate::sealed::Sealed {
    fn trigamma(self) -> Self;
}

impl TrigammaArg for f64 {
    #[inline]
    fn trigamma(self) -> Self {
        trigamma_f64(self)
    }
}

impl TrigammaArg for Complex<f64> {
    #[inline]
    fn trigamma(self) -> Self {
        trigamma_c64(self)
    }
}

/// Trigamma function $\psi\'(z)$ for real or complex input
///
/// This is the derivative of the [digamma function](crate::digamma) [^DLMF]:
///
/// $$
/// \psi\'(z) = \sum_{k=0}^\infty {1 \over (z+k)^2} .
/// $$
///
/// It equals [`polygamma(1, x)`](crate::polygamma) for real $x$.
///
/// [^DLMF]: NIST Digital Library of Mathematical Functions, <https://dlmf.nist.gov/5.15.E1>.
///
/// # See also
/// - [`digamma`](crate::digamma): Digamma function
/// - [`polygamma`](crate::polygamma): Polygamma function (n-th derivative)
#[must_use]
#[inline]
pub fn trigamma<T: TrigammaArg>(z: T) -> T {
    z.trigamma()
}
