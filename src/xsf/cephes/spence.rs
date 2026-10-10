//! Translated into pure Rust from `xsf/cephes/spence.h` (xsf v0.2.2).
//!
//! Translated into C++ by SciPy developers in 2024 from the Cephes Math Library Release 2.1
//! (January, 1989), Copyright 1985, 1987, 1989 by Stephen L. Moshier.
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.

use crate::xsf::cephes::polevl::polevl;

use core::f64::consts::PI;

const SPENCE_A: [f64; 8] = [
    4.651_285_860_739_900_3e-5,
    7.315_890_452_380_947e-3,
    1.338_476_395_783_090_3e-1,
    8.796_913_117_545_303e-1,
    2.711_498_511_965_534_6,
    4.256_971_560_081_218,
    3.297_713_409_852_251,
    1.0,
];

const SPENCE_B: [f64; 8] = [
    6.909_904_889_125_533e-4,
    2.540_437_639_325_444e-2,
    2.829_748_606_025_681e-1,
    1.411_725_977_518_310_6,
    3.638_005_333_451_370_7,
    5.032_788_801_433_17,
    3.547_713_409_852_251,
    1.0,
];

/// Dilogarithm
///
/// Computes the integral
///
/// ```text
///                    x
///                    -
///                   | | log t
/// spence(x)  =  -   |   ----- dt
///                 | |   t - 1
///                  -
///                  1
/// ```
///
/// for x >= 0.  A rational approximation gives the integral in
/// the interval (0.5, 1.5).  Transformation formulas for 1/x
/// and 1-x are employed outside the basic expansion range.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE      0,4         30000       3.8e-15     5.5e-16
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn spence(mut x: f64) -> f64 {
    let w;

    if x < 0.0 {
        return f64::NAN;
    }

    if x == 1.0 {
        return 0.0;
    }

    if x == 0.0 {
        return PI * PI / 6.0;
    }

    let mut flag = 0;

    if x > 2.0 {
        x = 1.0 / x;
        flag |= 2;
    }

    if x > 1.5 {
        w = (1.0 / x) - 1.0;
        flag |= 2;
    } else if x < 0.5 {
        w = -x;
        flag |= 1;
    } else {
        w = x - 1.0;
    }

    let mut y = -w * polevl(w, &SPENCE_A) / polevl(w, &SPENCE_B);

    if flag & 1 != 0 {
        y = (PI * PI) / 6.0 - x.ln() * (1.0 - x).ln() - y;
    }

    if flag & 2 != 0 {
        let z = x.ln();
        y = -0.5 * z * z - y;
    }

    y
}
