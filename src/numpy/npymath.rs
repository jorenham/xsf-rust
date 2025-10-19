//! <https://github.com/numpy/numpy/blob/4992eaf/numpy/_core/src/npymath/npy_math_internal.h.src>

use core::f64::consts::{LN_2, LOG2_E};

mod sealed {
    pub trait Sealed {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

pub trait LogAddExpArg: sealed::Sealed {
    fn npy_logaddexp(self, other: Self) -> Self;
    fn npy_logaddexp2(self, other: Self) -> Self;
}

impl LogAddExpArg for f32 {
    #[inline(always)]
    fn npy_logaddexp(self, other: Self) -> Self {
        (self as f64).npy_logaddexp(other as f64) as f32
    }

    #[inline(always)]
    fn npy_logaddexp2(self, other: Self) -> Self {
        (self as f64).npy_logaddexp2(other as f64) as f32
    }
}

impl LogAddExpArg for f64 {
    #[inline(always)]
    fn npy_logaddexp(self, other: Self) -> Self {
        if self == other {
            // Handles infinities of the same sign without warnings
            self + LN_2
        } else {
            let tmp = self - other;
            if tmp > 0.0 {
                self + (-tmp).exp().ln_1p()
            } else if tmp <= 0.0 {
                other + tmp.exp().ln_1p()
            } else {
                // NaNs
                tmp
            }
        }
    }

    #[inline(always)]
    fn npy_logaddexp2(self, other: Self) -> Self {
        if self == other {
            // Handles infinities of the same sign without warnings
            self + 1.0
        } else {
            let tmp = self - other;
            if tmp > 0.0 {
                self + (-tmp).exp2().ln_1p() * LOG2_E
            } else if tmp <= 0.0 {
                other + tmp.exp2().ln_1p() * LOG2_E
            } else {
                // NaNs
                tmp
            }
        }
    }
}

/// ln(e<sup>x</sup> + e<sup>y</sup>)
///
/// This function is useful in statistics where the calculated probabilities of events may be so
/// small as to exceed the range of normal floating point numbers. In such cases the logarithm of
/// the calculated probability is stored. This function allows adding probabilities stored in such
/// a fashion.
///
/// # See also
/// - [`logaddexp2`]
/// - [NumPy documentation](https://numpy.org/doc/stable/reference/generated/numpy.logaddexp.html)
///
#[inline]
pub fn logaddexp<T: LogAddExpArg>(x: T, y: T) -> T {
    x.npy_logaddexp(y)
}

/// log<sub>2</sub>(2<sup>x</sup> + 2<sup>y</sup>)
///
/// This function is useful in machine learning when the calculated probabilities of events may be
/// so small as to exceed the range of normal floating point numbers. In such cases the base-2
/// logarithm of the calculated probability can be used instead. This function allows adding
/// probabilities stored in such a fashion.
///
/// # See also
/// - [`logaddexp`]
/// - [NumPy documentation](https://numpy.org/doc/stable/reference/generated/numpy.logaddexp2.html)
///
#[inline]
pub fn logaddexp2<T: LogAddExpArg>(x: T, y: T) -> T {
    x.npy_logaddexp2(y)
}

#[cfg(test)]
mod tests {
    use crate::np_assert_allclose;

    #[test]
    fn test_logaddexp_f32() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [5.0, 4.0, 3.0, 2.0, 1.0];
        let z = [6.0, 6.0, 6.0, 6.0, 6.0];

        let xf = x.map(f64::log2);
        let yf = y.map(f64::log2);
        let zf = z.map(f64::log2);

        let zr: [f32; 5] = std::array::from_fn(|i| crate::logaddexp2(xf[i] as f32, yf[i] as f32));
        np_assert_allclose!(zr.map(|z| z as f64), zf, atol = 1.5e-7);
    }

    #[test]
    fn test_logaddexp_f64() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [5.0, 4.0, 3.0, 2.0, 1.0];
        let z = [6.0, 6.0, 6.0, 6.0, 6.0];

        let xf = x.map(f64::log2);
        let yf = y.map(f64::log2);
        let zf = z.map(f64::log2);

        let zr: [f64; 5] = std::array::from_fn(|i| crate::logaddexp2(xf[i], yf[i]));
        np_assert_allclose!(zr, zf, atol = 1e-15);
    }
}
