//! No jumps between neighbouring inputs, e.g. at a branch switch or an exact-parameter guard

// exact inputs, e.g. x = 1
#![allow(clippy::float_cmp)]

use num_complex::{Complex64, c64};
use proptest::prelude::*;
use proptest::sample::Index;

use super::strategies::{moderate, small, unit};
use super::symmetry::{int_distance, is_pole};

/// Relative steps `2^-k` from half the digits, `2^-MIN_K`, down to an ulp, `2^-MAX_K`
const MIN_K: u32 = f64::MANTISSA_DIGITS / 2;
const MAX_K: u32 = f64::MANTISSA_DIGITS - 1;

/// A jump must persist at a `2^-16` times smaller step, where the curvature of a smooth function
/// changes the middle value `2^-32` times less
const CONFIRM_K: u32 = 16;

/// Within `sqrt(eps)` of an integer, e.g. 3 + 1e-12, where terms with opposite poles at the
/// integers cancel more than half the digits
fn near_int(x: f64) -> bool {
    int_distance(x) <= f64::EPSILON.sqrt()
}

/// `near_int` for a pole of Gamma, e.g. -3 + 1e-12 or -1e-15
fn near_pole(x: f64) -> bool {
    x.round() <= 0.0 && near_int(x)
}

/// `f(x)` with `x[i]` moved by `-h`, `0` and `h`, for `h = 2^-k |x[i]|`, but at most half the
/// distance to the nearest integer if `x[i]` is not one, or `None` if these inputs are not
/// distinct or round to that integer.
///
/// Poles in the parameters are at integers, and a pole closer than `h` with a small residue would
/// look like a jump.
fn bracket<const N: usize, T: Into<Complex64>>(
    x: [f64; N],
    i: usize,
    k: u32,
    f: impl Fn([f64; N]) -> T,
) -> Option<[Complex64; 3]> {
    let offset = int_distance(x[i]);
    let h = x[i].abs() * (-f64::from(k)).exp2();
    let h = if offset > 0.0 { h.min(offset / 2.0) } else { h };
    let ys = [x[i] - h, x[i], x[i] + h];
    let off_int = offset == 0.0 || ys.iter().all(|&y| int_distance(y) > 0.0);
    (ys[0] < ys[1] && ys[1] < ys[2] && off_int).then(|| {
        ys.map(|t| {
            let mut y = x;
            y[i] = t;
            f(y).into()
        })
    })
}

/// `bracket` at step `k`, and at the `2^CONFIRM_K` times smaller step that confirms a jump
fn neighbours<const N: usize, T: Into<Complex64>>(
    x: [f64; N],
    i: usize,
    k: u32,
    f: impl Fn([f64; N]) -> T + Copy,
) -> [Option<[Complex64; 3]>; 2] {
    [k, (k + CONFIRM_K).min(MAX_K)].map(|k| bracket(x, i, k, f))
}

/// The middle value is not within a quarter of the digits of both neighbours, while these are
/// normal and agree to half the digits.
///
/// A smooth function lies between its neighbours, and steep ones, poles, and zeros make the
/// neighbours disagree.
fn jumps([lo, mid, hi]: [Complex64; 3]) -> bool {
    let normal = |w: Complex64| (f64::MIN_POSITIVE..f64::INFINITY).contains(&w.norm());
    let scale = lo.norm().max(hi.norm());
    let agree = f64::EPSILON.sqrt();
    let close = (mid - lo).norm().max((mid - hi).norm()) <= agree.sqrt() * scale;
    normal(lo) && normal(hi) && (hi - lo).norm() <= agree * scale && !close
}

/// No jump at both steps of `neighbours`, where a NaN or infinite middle value is a jump
fn continuous(brackets: [Option<[Complex64; 3]>; 2]) -> Result<(), TestCaseError> {
    match brackets.map(|w| w.filter(|&w| jumps(w))) {
        [Some([lo, mid, hi]), Some(_)] => {
            Err(TestCaseError::fail(format!("{lo:e}, {mid:e}, {hi:e}")))
        },
        _ => Ok(()),
    }
}

/// `property!(continuous, ...)` that perturbs one random argument as in `neighbours`, and skips
/// (rather than rejects) the inputs where `$skip` holds
macro_rules! continuous {
    ($name:ident($($arg:ident $(in $strategy:expr)?),+) $(unless $skip:expr)? => $call:expr) => {
        property!(
            continuous,
            $name($($arg $(in $strategy)?),+, i in any::<Index>(), k in MIN_K..=MAX_K) => {
                $(if $skip {
                    return Ok(());
                })?
                let x = [$($arg),+];
                neighbours(x, i.index(x.len()), k, |[$($arg),+]| $call)
            }
        );
    };
    ($name:ident($($arg:ident $(in $strategy:expr)?),+) $(unless $skip:expr)?) => {
        continuous!($name($($arg $(in $strategy)?),+) $(unless $skip)? => xsf::$name($($arg),+));
    };
}

/// A smooth extremum, where the curvature over a step of 0.16 looks like a jump
#[test]
fn bessel_j_extremum() {
    let x = [0.0, 52_707_179.318_687_3];
    let brackets = neighbours(x, 1, MIN_K, |[v, x]| xsf::bessel_j(v, x));
    assert!(brackets[0].is_some_and(jumps));
    continuous(brackets).unwrap();
}

// singular at the poles of c, and at x = 1 for c <= a + b; known: polynomials at x = 1 (#2198,
// fixed by #2214), next to the poles of c (#2241), a tiny c near x = 1 (#2242) and at x << 0
// (#2243), a tiny a or b at x << 0, and c = a or b at x < -1
continuous!(
    hyp2f1(a in moderate(), b in moderate(), c in moderate(), x)
        unless near_pole(c)
            || (x == 1.0 && (c - a - b <= 0.0 || is_pole(a) || is_pole(b)))
            || (a.abs().min(b.abs()) < f64::EPSILON && x < -1.0)
            || ((a == c || b == c) && x < -1.0)
);
// singular at the poles of c, and at z = 1 for c <= a + b; known: polynomials (#2130), Euler
// polynomials (#2107), near-integer a - b at |z| > 1 (#2104, #2244, #2245), next to the poles of
// c (#2246), and cancellation for larger parameters (#2196)
continuous!(
    hyp2f1_complex(a in small(), b in small(), c in small(), re, im)
        unless is_pole(a) || is_pole(b) || is_pole(c - a) || is_pole(c - b) || near_pole(c)
            || (re == 1.0 && c - a - b <= 0.0)
            || (near_int(a - b) && re.hypot(im) > 1.0)
        => xsf::hyp2f1(a, b, c, c64(re, im))
);
// known: finite at the poles b = 0, -1, ... (#2239), unlike the limit for integers b <= a <= 0,
// and wrong for x << 0 next to them (#2240); polynomials a = 0, -1, ... where e^x underflows
continuous!(
    hyp1f1(a in moderate(), b in moderate(), x)
        unless near_pole(b) || (is_pole(a) && x.exp() == 0.0)
);
// known: NaN for an integer a at a tiny x
continuous!(hypu(a in moderate(), b in moderate(), x) unless a.fract() == 0.0 && x < f64::EPSILON);

// known: AMOS loses all digits beyond 2^30
const MAX_BESSEL: f64 = (1 << 30) as f64;

// singular at the poles b = 0, -1, ...; known: beyond the Bessel cutoff at 2 sqrt(|x|) (#1472)
continuous!(hyp0f1(b in moderate(), x) unless is_pole(b) || 2.0 * x.abs().sqrt() > MAX_BESSEL);

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
