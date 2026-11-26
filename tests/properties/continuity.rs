//! No jumps between neighbouring inputs, e.g. at a branch switch or an exact-parameter guard

use num_complex::{Complex64, c64};
use proptest::prelude::*;
use proptest::sample::Index;

use super::strategies::{moderate, small, unit};
use super::symmetry::is_pole;

/// Within `sqrt(eps)` of an integer, e.g. 3 + 1e-12, where terms with opposite poles at the
/// integers cancel more than half the digits
fn near_int(x: f64) -> bool {
    (x - x.round()).abs() <= f64::EPSILON.sqrt()
}

/// `near_int` for a pole of Gamma, e.g. -3 + 1e-12 or -1e-15
fn near_pole(x: f64) -> bool {
    x.round() <= 0.0 && near_int(x)
}

/// `f(x)` with `x[i]` moved by `-h`, `0` and `h`, for `h = 2^-k |x[i]|`, but at most half the
/// distance to the nearest integer if `x[i]` is not one.
///
/// Poles in the parameters are at integers, and a pole closer than `h` with a small residue would
/// look like a jump.
fn neighbours<const N: usize, T: Into<Complex64>>(
    x: [f64; N],
    i: Index,
    k: i32,
    f: impl Fn([f64; N]) -> T,
) -> [Complex64; 3] {
    let i = i.index(N);
    let offset = (x[i] - x[i].round()).abs();
    let h = x[i].abs() * f64::from(-k).exp2();
    let h = if offset > 0.0 { h.min(offset / 2.0) } else { h };
    [-h, 0.0, h].map(|t| {
        let mut y = x;
        y[i] += t;
        f(y).into()
    })
}

/// A finite middle value within a quarter of the digits of both neighbours, if these are normal
/// and agree to half the digits.
///
/// A smooth function lies between its neighbours, and steep ones, poles, and zeros make the
/// neighbours disagree. Exact poles, e.g. of `hyp2f1` at `c = -2`, can be infinite or NaN.
fn continuous([lo, mid, hi]: [Complex64; 3]) -> Result<(), TestCaseError> {
    let normal = |w: Complex64| w.is_finite() && w.norm() >= f64::MIN_POSITIVE;
    let scale = lo.norm().max(hi.norm());
    let agree = f64::EPSILON.sqrt();
    if mid.is_finite() && normal(lo) && normal(hi) && (hi - lo).norm() <= agree * scale {
        let jump = (mid - lo).norm().max((mid - hi).norm());
        prop_assert!(jump <= agree.sqrt() * scale, "{lo:e}, {mid:e}, {hi:e}");
    }
    Ok(())
}

/// `property!(continuous, ...)` that perturbs one random argument as in `neighbours`, and skips
/// (rather than rejects) the inputs where `$skip` holds
macro_rules! continuous {
    ($name:ident($($arg:ident $(in $strategy:expr)?),+) $(unless $skip:expr)? => $call:expr) => {
        property!(
            continuous,
            $name($($arg $(in $strategy)?),+, i in any::<Index>(), k in 26..=52_i32) => {
                $(if $skip {
                    return Ok(());
                })?
                neighbours([$($arg),+], i, k, |[$($arg),+]| $call)
            }
        );
    };
    ($name:ident($($arg:ident $(in $strategy:expr)?),+) $(unless $skip:expr)?) => {
        continuous!($name($($arg $(in $strategy)?),+) $(unless $skip)? => xsf::$name($($arg),+));
    };
}

// known: polynomials at x = 1 (#2198, TODO for c > 0), for c < 0 (TODO), and near x = 1 for a tiny
// c (TODO); near the poles of c, e.g. a tiny c for integers a and b at x << 0 (TODO)
continuous!(
    hyp2f1(a in moderate(), b in moderate(), c in moderate(), x)
        unless is_pole(a) || is_pole(b) || near_pole(c)
);
// known: polynomials (#2130), Euler polynomials (#2107), near-integer a - b at |z| > 1 (#2104,
// TODO for an integer a - b), near the poles of c (TODO), and cancellation for larger parameters
// (#2068, #2196)
continuous!(
    hyp2f1_complex(a in small(), b in small(), c in small(), re, im)
        unless is_pole(a) || is_pole(b) || is_pole(c - a) || is_pole(c - b)
            || near_int(a - b) || near_pole(c)
        => xsf::hyp2f1(a, b, c, c64(re, im))
);
// known: finite at the poles b = 0, -1, ... (TODO), unlike the limit for integers b <= a <= 0, and
// wrong for x << 0 near them (TODO); the underflow of e^x after Kummer's transformation (#1824)
continuous!(
    hyp1f1(a in moderate(), b in moderate(), x)
        unless near_pole(b) || (is_pole(a) && x.exp() == 0.0)
);
continuous!(hypu(a in moderate(), b in moderate(), x));

// known: AMOS loses all digits beyond 2^30 (accepted, scipy/xsf#216)
const MAX_BESSEL: f64 = (1 << 30) as f64;

continuous!(hyp0f1(b in moderate(), x) unless x.abs().sqrt() > MAX_BESSEL);

continuous!(gamma(x));
continuous!(loggamma(x));
continuous!(rgamma(x));

continuous!(bessel_j(v in moderate(), x) unless x.abs() > MAX_BESSEL);
continuous!(bessel_y(v in moderate(), x) unless x.abs() > MAX_BESSEL);
continuous!(bessel_i(v in moderate(), x) unless x.abs() > MAX_BESSEL);
continuous!(bessel_k(v in moderate(), x) unless x.abs() > MAX_BESSEL);

continuous!(betainc(a in moderate(), b in moderate(), x in unit()));
continuous!(betaincinv(a in moderate(), b in moderate(), y in unit()));
continuous!(gammainc(a in moderate(), x));
continuous!(gammaincinv(a in moderate(), p in unit()));
continuous!(gdtrib(a in moderate(), p in unit(), x));
