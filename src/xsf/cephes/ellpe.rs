//! Translated into pure Rust from `xsf/cephes/ellpe.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library, Release 2.1:  February, 1989
//! Copyright 1984, 1987, 1989 by Stephen L. Moshier
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//!
//! Feb, 2002:  altered by Travis Oliphant
//! so that it is called with argument m
//! (which gets immediately converted to m1 = 1-m)
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the domain error condition and return value are corrected in the comments.

use crate::xsf::cephes::polevl::polevl;

const ELLPE_P: [f64; 11] = [
    1.535_525_773_010_133e-4,
    2.508_884_921_636_020_4e-3,
    8.687_868_165_658_896e-3,
    1.073_509_490_560_761_9e-2,
    7.773_954_925_167_871e-3,
    7.583_952_894_135_147e-3,
    1.156_884_368_105_741_2e-2,
    2.183_179_960_155_572_4e-2,
    5.680_519_456_178_606e-2,
    4.431_471_805_609_908_4e-1,
    1.0,
];

const ELLPE_Q: [f64; 10] = [
    3.279_548_985_764_858_5e-5,
    1.009_627_926_793_567_2e-3,
    6.506_094_899_769_275e-3,
    1.688_621_639_933_113_3e-2,
    2.617_697_424_544_936_4e-2,
    3.348_339_048_882_249e-2,
    4.271_809_265_189_315e-2,
    5.859_366_344_711_01e-2,
    9.374_999_971_976_443e-2,
    2.499_999_999_998_883e-1,
];

/// Complete elliptic integral of the second kind
///
/// Approximates the integral
///
/// ```text
///            pi/2
///             -
///            | |                 2
/// E(m)  =    |    sqrt( 1 - m sin t ) dt
///          | |
///           -
///            0
/// ```
///
/// Where m = 1 - m1, using the approximation
///
/// ```text
///      P(x)  -  x log x Q(x).
/// ```
///
/// Though there are no singularities, the argument m1 is used
/// internally rather than m for compatibility with ellpk().
///
/// E(1) = 1; E(0) = pi/2.
///
/// ACCURACY:
///
/// ```text
///                      Relative error:
/// arithmetic   domain     # trials      peak         rms
///    IEEE       0, 1       20000       2.2e-16     7.2e-17
/// ```
///
/// ERROR MESSAGES:
///
/// ```text
///   message         condition      value returned
/// ellpe domain         m>1               NaN
/// ```
#[allow(clippy::float_cmp)]
pub(crate) fn ellpe(mut x: f64) -> f64 {
    x = 1.0 - x;
    if x <= 0.0 {
        if x == 0.0 {
            return 1.0;
        }
        return f64::NAN;
    }
    if x > 1.0 {
        return ellpe(1.0 - 1.0 / x) * x.sqrt();
    }
    polevl(x, &ELLPE_P) - x.ln() * (x * polevl(x, &ELLPE_Q))
}
