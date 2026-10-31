//! Carlson symmetric elliptic integrals, translated into pure Rust from
//! `xsf/cpu/ellint_carlson.h` and `xsf/ellint_carlson/*.hh` (xsf main).
//!
//! The C++ functions return an `ExitStatus` and write the result to an output argument; here they
//! return both as a tuple. `set_error` is a no-op in our build of xsf, so the status is discarded.

mod argcheck;
mod arith_aux;
mod arithmetic;
mod common;
mod rc;
mod rd;
mod rf;
mod rg;
mod rj;
mod typing;

use num_complex::Complex64;

use typing::Num;

const ELLINT_CARLSON_RERR: f64 = 5e-16;

fn elliprc_impl<T: Num>(x: T, y: T) -> T {
    rc::rc(x, y, ELLINT_CARLSON_RERR).1
}

fn elliprd_impl<T: Num>(x: T, y: T, z: T) -> T {
    rd::rd(x, y, z, ELLINT_CARLSON_RERR).1
}

fn elliprf_impl<T: Num>(x: T, y: T, z: T) -> T {
    rf::rf(x, y, z, ELLINT_CARLSON_RERR).1
}

fn elliprg_impl<T: Num>(x: T, y: T, z: T) -> T {
    rg::rg(x, y, z, ELLINT_CARLSON_RERR).1
}

fn elliprj_impl<T: Num>(x: T, y: T, z: T, p: T) -> T {
    rj::rj(x, y, z, p, ELLINT_CARLSON_RERR, false).1
}

pub trait CarlsonArg: crate::sealed::Sealed + Sized {
    fn elliprc(self, y: Self) -> Self;
    fn elliprd(self, y: Self, z: Self) -> Self;
    fn elliprf(self, y: Self, z: Self) -> Self;
    fn elliprg(self, y: Self, z: Self) -> Self;
    fn elliprj(self, y: Self, z: Self, p: Self) -> Self;
}

impl CarlsonArg for f64 {
    #[inline]
    fn elliprc(self, y: Self) -> Self {
        elliprc_impl(self, y)
    }
    #[inline]
    fn elliprd(self, y: Self, z: Self) -> Self {
        elliprd_impl(self, y, z)
    }
    #[inline]
    fn elliprf(self, y: Self, z: Self) -> Self {
        elliprf_impl(self, y, z)
    }
    #[inline]
    fn elliprg(self, y: Self, z: Self) -> Self {
        elliprg_impl(self, y, z)
    }
    #[inline]
    fn elliprj(self, y: Self, z: Self, p: Self) -> Self {
        elliprj_impl(self, y, z, p)
    }
}

impl CarlsonArg for Complex64 {
    #[inline]
    fn elliprc(self, y: Self) -> Self {
        elliprc_impl(self, y)
    }
    #[inline]
    fn elliprd(self, y: Self, z: Self) -> Self {
        elliprd_impl(self, y, z)
    }
    #[inline]
    fn elliprf(self, y: Self, z: Self) -> Self {
        elliprf_impl(self, y, z)
    }
    #[inline]
    fn elliprg(self, y: Self, z: Self) -> Self {
        elliprg_impl(self, y, z)
    }
    #[inline]
    fn elliprj(self, y: Self, z: Self, p: Self) -> Self {
        elliprj_impl(self, y, z, p)
    }
}

/// Degenerate symmetric elliptic integral $R_C(x, y)$ for real or complex input
///
/// $$
/// R_C(x, y) = {1 \over 2} \int_0^\infty (t + x)^{-1/2} (t + y)^{-1} \dd t
/// = R_F(x, y, y)
/// $$
///
/// For $y < 0$, the Cauchy principal value is returned.
///
/// Corresponds to [`scipy.special.elliprc`][elliprc] in SciPy
///
/// [elliprc]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprc.html
///
/// # See also
/// - [`elliprf`]: Completely-symmetric elliptic integral of the first kind
#[must_use]
#[inline]
pub fn elliprc<T: CarlsonArg>(x: T, y: T) -> T {
    x.elliprc(y)
}

/// Symmetric elliptic integral of the second kind $R_D(x, y, z)$ for real or complex input
///
/// $$
/// R_D(x, y, z) = {3 \over 2} \int_0^\infty
/// \left[(t + x) (t + y)\right]^{-1/2} (t + z)^{-3/2} \dd t
/// = R_J(x, y, z, z)
/// $$
///
/// Corresponds to [`scipy.special.elliprd`][elliprd] in SciPy
///
/// [elliprd]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprd.html
///
/// # See also
/// - [`elliprj`]: Symmetric elliptic integral of the third kind
#[must_use]
#[inline]
pub fn elliprd<T: CarlsonArg>(x: T, y: T, z: T) -> T {
    x.elliprd(y, z)
}

/// Completely-symmetric elliptic integral of the first kind $R_F(x, y, z)$ for real or complex
/// input
///
/// $$
/// R_F(x, y, z) = {1 \over 2} \int_0^\infty \left[(t + x) (t + y) (t + z)\right]^{-1/2} \dd t
/// $$
///
/// Corresponds to [`scipy.special.elliprf`][elliprf] in SciPy
///
/// [elliprf]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprf.html
///
/// # See also
/// - [`elliprg`]: Completely-symmetric elliptic integral of the second kind
/// - [`ellipk`](crate::ellipk): Complete elliptic integral of the first kind, $K(m) = R_F(0, 1 - m,
///   1)$
#[must_use]
#[inline]
pub fn elliprf<T: CarlsonArg>(x: T, y: T, z: T) -> T {
    x.elliprf(y, z)
}

/// Completely-symmetric elliptic integral of the second kind $R_G(x, y, z)$ for real or complex
/// input
///
/// $$
/// R_G(x, y, z) = {1 \over 4} \int_0^\infty \left[(t + x) (t + y) (t + z)\right]^{-1/2}
/// \left({x \over t + x} + {y \over t + y} + {z \over t + z}\right) t \dd t
/// $$
///
/// Corresponds to [`scipy.special.elliprg`][elliprg] in SciPy
///
/// [elliprg]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprg.html
///
/// # See also
/// - [`elliprf`]: Completely-symmetric elliptic integral of the first kind
/// - [`ellipe`](crate::ellipe): Complete elliptic integral of the second kind,
///   $E(m) = 2 R_G(0, 1 - m, 1)$
#[must_use]
#[inline]
pub fn elliprg<T: CarlsonArg>(x: T, y: T, z: T) -> T {
    x.elliprg(y, z)
}

/// Symmetric elliptic integral of the third kind $R_J(x, y, z, p)$ for real or complex input
///
/// $$
/// R_J(x, y, z, p) = {3 \over 2} \int_0^\infty
/// \left[(t + x) (t + y) (t + z)\right]^{-1/2} (t + p)^{-1} \dd t
/// $$
///
/// For real $x, y, z \ge 0$ and $p < 0$, the Cauchy principal value is returned.
///
/// Corresponds to [`scipy.special.elliprj`][elliprj] in SciPy
///
/// [elliprj]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprj.html
///
/// # See also
/// - [`elliprd`]: Symmetric elliptic integral of the second kind
#[must_use]
#[inline]
pub fn elliprj<T: CarlsonArg>(x: T, y: T, z: T, p: T) -> T {
    x.elliprj(y, z, p)
}

#[cfg(test)]
mod tests {
    use core::f64::consts::{LN_2, PI};

    use num_complex::{Complex64, c64};

    const INF: f64 = f64::INFINITY;
    const ZERO: Complex64 = Complex64::new(0.0, 0.0);
    const ONE: Complex64 = Complex64::new(1.0, 0.0);
    const I: Complex64 = Complex64::new(0.0, 1.0);

    // The cases are from `scipy.special.tests.test_basic.TestCarlson` (also used in the xsf tests),
    // with the reference values recomputed with mpmath.

    #[test]
    fn test_elliprc() {
        let actual = [
            crate::elliprc(0.0, 0.25),
            crate::elliprc(2.25, 2.0),
            crate::elliprc(0.25, -2.0), // Cauchy principal value
        ];
        let expected = [PI, LN_2, 0.231_049_060_186_648_45];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        let actual = [
            crate::elliprc(ZERO, I),
            crate::elliprc(-I, I),
            crate::elliprc(I, -ONE),
        ];
        let expected = [
            c64(1.110_720_734_539_591_5, -1.110_720_734_539_591_5),
            c64(1.226_084_956_907_219_9, -0.344_711_369_887_676_8),
            c64(0.777_785_969_204_473_8, 0.198_324_849_934_287_73),
        ];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        assert_eq!(crate::elliprc(1.0, INF), 0.0);
        assert!(crate::elliprc(1.0, 0.0).is_nan());
        assert_eq!(crate::elliprc(ONE, c64(1.0, INF)), ZERO);
    }

    #[test]
    fn test_elliprd() {
        let actual = [crate::elliprd(0.0, 2.0, 1.0), crate::elliprd(2.0, 3.0, 4.0)];
        let expected = [1.797_210_352_103_388_4, 0.165_105_272_942_610_54];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        let actual = [
            crate::elliprd(I, -I, c64(2.0, 0.0)),
            crate::elliprd(ZERO, I, -I),
            crate::elliprd(ZERO, c64(-1.0, 1.0), I),
            crate::elliprd(c64(-2.0, -1.0), -I, c64(-1.0, 1.0)),
        ];
        let expected = [
            c64(0.659_338_541_542_197_7, 0.0),
            c64(1.270_819_627_190_968_6, 2.781_112_015_952_057_7),
            c64(-1.857_723_543_923_906, -0.961_934_508_888_385_6),
            c64(1.824_902_739_370_380_5, -1.221_847_578_482_703_6),
        ];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        assert_eq!(crate::elliprd(1.0, 1.0, INF), 0.0);
        assert_eq!(crate::elliprd(1.0, 1.0, 0.0), INF);
        assert!(crate::elliprd(1.0, 1.0, -1.0).is_nan());
        assert!(crate::elliprd(1.0, 1.0, -f64::MIN_POSITIVE / 2.0).is_nan());
        assert_eq!(crate::elliprd(ZERO, ONE, ZERO).re, INF);
        assert!(crate::elliprd(ONE, ONE, -ONE).re.is_nan());
    }

    #[test]
    fn test_elliprf() {
        let actual = [crate::elliprf(0.0, 1.0, 2.0), crate::elliprf(2.0, 3.0, 4.0)];
        let expected = [1.311_028_777_146_059_8, 0.584_082_841_677_151_7];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        let actual = [
            crate::elliprf(I, -I, ZERO),
            crate::elliprf(c64(-1.0, 1.0), I, ZERO),
            crate::elliprf(I, -I, c64(2.0, 0.0)),
            crate::elliprf(c64(-1.0, 1.0), I, c64(1.0, -1.0)),
        ];
        let expected = [
            c64(1.854_074_677_301_371_9, 0.0),
            c64(0.796_125_865_842_339_1, -1.213_856_669_836_496),
            c64(1.044_144_565_406_436, 0.0),
            c64(0.939_120_502_186_193_7, -0.532_962_520_186_352_7),
        ];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        assert_eq!(crate::elliprf(1.0, INF, 1.0), 0.0);
        assert_eq!(crate::elliprf(0.0, 1.0, 0.0), INF);
        assert!(crate::elliprf(1.0, 1.0, -1.0).is_nan());
        assert_eq!(crate::elliprf(c64(INF, 0.0), ZERO, ONE), ZERO);
        assert!(crate::elliprf(ONE, ONE, c64(-INF, 1.0)).re.is_nan());

        // K(m) = R_F(0, 1 - m, 1)
        for m in [-1e300, -1e10, -1.0, 0.0, 0.5, 0.99, 1.0 - 1e-12] {
            crate::np_assert_allclose!(
                [crate::elliprf(0.0, 1.0 - m, 1.0)],
                [crate::ellipk(m)],
                rtol = 1e-15
            );
        }
    }

    #[test]
    fn test_elliprg() {
        let actual = [
            crate::elliprg(0.0, 16.0, 16.0),
            crate::elliprg(2.0, 3.0, 4.0),
            crate::elliprg(0.0, 0.0796, 4.0),
        ];
        let expected = [PI, 1.725_503_028_069_227_8, 1.028_475_809_028_804];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        let actual = [
            crate::elliprg(ZERO, I, -I),
            crate::elliprg(c64(-1.0, 1.0), I, ZERO),
            crate::elliprg(-I, c64(-1.0, 1.0), I),
        ];
        let expected = [
            c64(0.423_606_542_396_989_5, 0.0),
            c64(0.446_605_916_770_183_74, 0.707_683_523_575_153_8),
            c64(0.360_233_921_844_733_1, 0.403_486_234_017_221_1),
        ];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        assert_eq!(crate::elliprg(0.0, 0.0, 0.0), 0.0);
        assert_eq!(crate::elliprg(1.0, INF, 1.0), INF);
        assert_eq!(crate::elliprg(c64(INF, 0.0), ONE, ONE).re, INF);

        // E(m) = 2 R_G(0, 1 - m, 1)
        for m in [-1e300, -1e10, -1.0, 0.0, 0.5, 0.99, 1.0] {
            crate::np_assert_allclose!(
                [2.0 * crate::elliprg(0.0, 1.0 - m, 1.0)],
                [crate::ellipe(m)],
                rtol = 1e-15
            );
        }
    }

    #[test]
    fn test_elliprj() {
        let actual = [
            crate::elliprj(0.0, 1.0, 2.0, 3.0),
            crate::elliprj(2.0, 3.0, 4.0, 5.0),
            // Cauchy principal values
            crate::elliprj(2.0, 3.0, 4.0, -0.5),
            crate::elliprj(2.0, 3.0, 4.0, -5.0),
            // "difficult arguments"
            crate::elliprj(
                6.483_625_725_195_452e-8,
                1.164_913_652_819_688_6e-27,
                3.676_734_016_716_8e13,
                0.493_704_617_023_468,
            ),
            crate::elliprj(
                14.375_105_857_849_121,
                9.993_988_969_725_365e-11,
                1.728_442_622_699_44e-26,
                5.898_871_222_598_245e-6,
            ),
        ];
        let expected = [
            0.776_886_237_785_823_3,
            0.142_975_796_671_567_54,
            0.247_238_197_030_515_64,
            -0.127_112_300_429_639_1,
            8.634_269_206_442_419e-6,
            829_774.142_480_162_7,
        ];
        crate::np_assert_allclose!(&actual, &expected, rtol = 5e-15);

        let actual = [
            crate::elliprj(c64(2.0, 0.0), c64(3.0, 0.0), c64(4.0, 0.0), c64(-1.0, 1.0)),
            crate::elliprj(I, -I, ZERO, c64(2.0, 0.0)),
            crate::elliprj(c64(-1.0, 1.0), c64(-1.0, -1.0), ONE, c64(2.0, 0.0)),
            crate::elliprj(I, -I, ZERO, c64(1.0, -1.0)),
            crate::elliprj(c64(-1.0, 1.0), c64(-1.0, -1.0), ONE, c64(-3.0, 1.0)),
        ];
        let expected = [
            c64(0.136_139_458_277_705_35, -0.382_075_616_244_271_67),
            c64(1.649_001_166_271_088_4, 0.0),
            c64(0.941_483_588_412_202_4, 0.0),
            c64(1.826_011_522_900_931_7, 1.229_066_190_864_347_2),
            c64(-0.611_279_708_120_281_8, -1.068_403_839_000_680_7),
        ];
        crate::np_assert_allclose!(&actual, &expected, rtol = 1e-15);

        assert_eq!(crate::elliprj(1.0, 1.0, 1.0, INF), 0.0);
        assert_eq!(crate::elliprj(1.0, 1.0, INF, 1.0), 0.0);
        assert!(crate::elliprj(-1.0, 1.0, 1.0, 1.0).is_nan());
        assert!(crate::elliprj(1.0, 0.0, 0.0, 0.0).is_nan());
    }
}
