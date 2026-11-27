//! No jumps between neighbouring inputs, e.g. at a branch switch or an exact-parameter guard

// exact inputs, e.g. x = 1
#![allow(clippy::float_cmp)]

use num_complex::{Complex64, c64};
use proptest::prelude::*;
use proptest::sample::Index;

use super::strategies::{arg, balanced, imag, moderate, near_threshold_z, small, unit};
use super::{NEAR_POLE, almost_int, config, int_distance, is_pole, pole_distance};

/// Relative steps `2^-k` from half the digits, `2^-MIN_K`, down to an ulp, `2^-MAX_K`
const MIN_K: u32 = f64::MANTISSA_DIGITS / 2;
const MAX_K: u32 = f64::MANTISSA_DIGITS - 1;

/// A jump must persist at a `2^-16` times smaller step, where the curvature of a smooth function
/// changes the middle value `2^-32` times less
const CONFIRM_K: u32 = 16;

/// `f` at `x - h`, `x` and `x + h`
type Bracket = [Complex64; 3];

/// `f(x)` with `x[i]` moved by `-h`, `0` and `h`, for `h = 2^-k scale`, but at most half the
/// distance to the nearest integer if `x[i]` is not one, or `None` if these inputs are not finite
/// and distinct, or round to that integer.
///
/// Poles in the parameters are at integers, and a pole closer than `h` with a small residue would
/// look like a jump.
fn bracket<const N: usize, T: Into<Complex64>>(
    x: [f64; N],
    i: usize,
    scale: f64,
    k: u32,
    f: impl Fn([f64; N]) -> T,
) -> Option<Bracket> {
    let offset = int_distance(x[i]);
    let h = scale * (-f64::from(k)).exp2();
    let h = if offset > 0.0 { h.min(offset / 2.0) } else { h };
    let ys = [x[i] - h, x[i], x[i] + h];
    let valid = |y: f64| y.is_finite() && (offset == 0.0 || int_distance(y) > 0.0);
    (ys[0] < ys[1] && ys[1] < ys[2] && ys.into_iter().all(valid)).then(|| {
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
    scale: f64,
    k: u32,
    f: impl Fn([f64; N]) -> T + Copy,
) -> [Option<Bracket>; 2] {
    [k, k + CONFIRM_K].map(|k| bracket(x, i, scale, k, f))
}

/// The middle value is not within a quarter of the digits of both neighbours, while these are
/// normal and agree to half the digits, relative to the larger neighbour in `wide`.
///
/// A smooth function lies between its neighbours, and steep ones, poles, and zeros make the
/// neighbours disagree.
fn jumps([lo, mid, hi]: Bracket, [wide_lo, _, wide_hi]: Bracket) -> bool {
    let scale = wide_lo.norm().max(wide_hi.norm());
    let normal = |w: Complex64| (f64::MIN_POSITIVE..f64::INFINITY).contains(&w.norm());
    let agree = f64::EPSILON.sqrt();
    let close = (mid - lo).norm().max((mid - hi).norm()) <= agree.sqrt() * scale;
    normal(lo) && normal(hi) && (hi - lo).norm() <= agree * scale && !close
}

/// No jump at both steps of `neighbours`, relative to the first, so that a smooth zero, where the
/// middle value approaches its neighbours with the step, is not one. A NaN or infinite middle
/// value is a jump.
fn continuous([wide, narrow]: [Option<Bracket>; 2]) -> Result<(), TestCaseError> {
    match (wide, narrow) {
        (Some(w @ [lo, mid, hi]), Some(n)) if jumps(w, w) && jumps(n, w) => {
            Err(TestCaseError::fail(format!("{lo:e}, {mid:e}, {hi:e}")))
        },
        _ => Ok(()),
    }
}

/// `property!(continuous, ...)` that perturbs one random argument `x` as in `neighbours`, with
/// `scale = |x|`, and skips (rather than rejects) the inputs where `$skip` holds
macro_rules! continuous {
    ($name:ident($($arg:ident $(in $strategy:expr)?),+) $(unless $skip:expr)? => $call:expr) => {
        property!(
            continuous,
            $name(
                $($arg $(in $strategy)?),+, i in any::<Index>(), k in MIN_K..=MAX_K - CONFIRM_K
            ) => {
                $(if $skip {
                    return Ok(());
                })?
                let x = [$($arg),+];
                let i = i.index(x.len());
                neighbours(x, i, x[i].abs(), k, |[$($arg),+]| $call)
            }
        );
    };
    ($name:ident($($arg:ident $(in $strategy:expr)?),+) $(unless $skip:expr)?) => {
        continuous!($name($($arg $(in $strategy)?),+) $(unless $skip)? => xsf::$name($($arg),+));
    };
}

proptest! {
    #![proptest_config(config())]

    // a nearly integer c - a - b
    #[test]
    fn hyp2f1_balanced(
        [a, b, c] in balanced(), x in arg(), i in any::<Index>(), k in MIN_K..=MAX_K - CONFIRM_K
    ) {
        if known_hyp2f1_balanced_x(a, b, c, x) {
            return Ok(());
        }
        let p = [a, b, c, x];
        let i = i.index(p.len());
        continuous(neighbours(p, i, p[i].abs(), k, |[a, b, c, x]| xsf::hyp2f1(a, b, c, x)))?;
    }

    // a nearly integer c - a - b, and z next to the switches between the transformations
    #[test]
    fn hyp2f1_complex_balanced(
        [a, b, c] in balanced(),
        (re, im) in near_threshold_z(),
        i in any::<Index>(),
        k in MIN_K..=MAX_K - CONFIRM_K
    ) {
        if known_hyp2f1_complex_balanced(a, b, c, re, im) {
            return Ok(());
        }
        let p = [a, b, c, re, im];
        let i = i.index(p.len());
        // relative to |z| at Re z = 0, so that it crosses that switch
        let scale = if i == 3 && re == 0.0 { im.abs() } else { p[i].abs() };
        continuous(neighbours(p, i, scale, k, |[a, b, c, re, im]| {
            xsf::hyp2f1(a, b, c, c64(re, im))
        }))?;
    }
}

/// A jump at the first step of `neighbours` that the second does not confirm
fn assert_unconfirmed(brackets: [Option<Bracket>; 2]) {
    assert!(brackets[0].is_some_and(|w| jumps(w, w)));
    continuous(brackets).unwrap();
}

/// A smooth extremum, where the curvature over a step of 0.16 looks like a jump
#[test]
fn bessel_j_extremum() {
    let x = [0.0, 52_707_179.318_687_3];
    assert_unconfirmed(neighbours(x, 1, x[1], MIN_K, |[v, x]| xsf::bessel_j(v, x)));
}

/// A smooth double zero in c, from the two reciprocal gammas `1 / (Gamma(c - a) Gamma(c - b))`
/// of Gauss's sum
#[test]
fn hyp2f1_double_zero() {
    let x = [-1.5, -1.5, -2.5, 1.0];
    assert_unconfirmed(neighbours(x, 2, -x[2], 30, |[a, b, c, x]| xsf::hyp2f1(a, b, c, x)));
}

/// The known failures of `hyp2f1` for real arguments
fn known_hyp2f1(a: f64, b: f64, c: f64, x: f64) -> bool {
    pole_distance(c) <= NEAR_POLE
        || (x == 1.0 && (c - a - b <= 0.0 || is_pole(a) || is_pole(b)))
        || (a.abs().min(b.abs()) < f64::EPSILON && x < -1.0)
        || ((a == c || b == c) && x < -1.0)
}

/// The known failures of `hyp2f1` for real arguments and a nearly integer `c - a - b`
fn known_hyp2f1_balanced_x(a: f64, b: f64, c: f64, x: f64) -> bool {
    // a terminating `p` with `c - q` next to an integer in `(p, 0]` near x = 1 (#2330)
    let terminating = |p: f64, q: f64| {
        is_pole(p) && int_distance(c - q) <= NEAR_POLE && (p + 1.0..=0.0).contains(&(c - q).round())
    };
    known_hyp2f1(a, b, c, x)
        // for 0 < c - a - b <= sqrt(eps), (1 - x)^(c - a - b) is 0 at x = 1 but about 1 at its
        // neighbours
        || (x == 1.0 && c - a - b <= f64::EPSILON.sqrt())
        || (x > 0.9 && (terminating(a, b) || terminating(b, a)))
        // c < 0 at a subnormal x (#2255, #2332)
        || (x.is_subnormal() && c < 0.0)
}

// singular at the poles of c, and at x = 1 for c <= a + b; known: polynomials at x = 1 (#2198,
// fixed by #2214), next to the poles of c (#2241), a tiny c near x = 1 (#2242) and at x << 0
// (#2243), a tiny a or b at x << 0 (#2253), and c = a or b at x < -1 (#2254)
continuous!(
    hyp2f1(a in moderate(), b in moderate(), c in moderate(), x) unless known_hyp2f1(a, b, c, x)
);

/// The known failures of `hyp2f1` for complex `z`
fn known_hyp2f1_complex(a: f64, b: f64, c: f64, re: f64, im: f64) -> bool {
    is_pole(a)
        || is_pole(b)
        || is_pole(c - a)
        || is_pole(c - b)
        || pole_distance(c) <= NEAR_POLE
        || (pole_distance(a).min(pole_distance(b)) <= NEAR_POLE && almost_int(c - a - b))
        || (re == 1.0 && im.abs() < f64::EPSILON.sqrt() && c - a - b <= 0.0)
        || (almost_int(a - b) && re.hypot(im) >= 1.1)
}

/// The largest `pole_distance` of a or b in the known failures #2334 of the complex `hyp2f1`,
/// e.g. 1e-5
const NEAR_POLE_COMPLEX: f64 = 1e-4;

/// The known failures of the complex `hyp2f1` for a nearly integer `c - a - b`
fn known_hyp2f1_complex_balanced(a: f64, b: f64, c: f64, re: f64, im: f64) -> bool {
    let near_pole = |t: f64| pole_distance(t) <= NEAR_POLE_COMPLEX;
    let r = re.hypot(im);
    let gap = c - a - b;
    known_hyp2f1_complex(a, b, c, re, im)
        // wrong or NaN for a or b next to a pole at 0.9 <= |z| <= 1.1 (#2334)
        || ((near_pole(a) || near_pole(b)) && almost_int(gap) && (0.9..=1.1).contains(&r))
        // wrong for negative a and b with c = a + b (#2335)
        || (gap == 0.0 && a.max(b) < 0.0 && c < -32.0)
        // wrong for a or b above 64 and c - a - b next to a negative integer, near |z| = 1 and
        // |1 - z| < 0.9 (#2339)
        || (a.min(b) > 0.0
            && a.max(b) > 64.0
            && int_distance(gap) <= NEAR_POLE
            && gap.round() < 0.0
            && (0.9..1.1).contains(&r)
            && (1.0 - re).hypot(im) < 0.9)
}

// singular at the poles of c, and at z = 1 for c <= a + b; known: polynomials (#2130), Euler
// polynomials (#2107, #2247), near-integer a - b at |z| >= 1.1 (#2104, #2244, #2245), next to the
// poles of c (#2246), a or b next to a pole for a near-integer c - a - b (#2256), and
// cancellation for larger parameters (#2196)
continuous!(
    hyp2f1_complex(a in small(), b in small(), c in small(), re, im in imag())
        unless known_hyp2f1_complex(a, b, c, re, im)
        => xsf::hyp2f1(a, b, c, c64(re, im))
);

// known: finite at the poles b = 0, -1, ... (#2239), unlike the limit for integers b <= a <= 0,
// and wrong for x << 0 next to them (#2240); polynomials a = 0, -1, ... at x << 0 (#2251)
continuous!(
    hyp1f1(a in moderate(), b in moderate(), x)
        unless is_pole(b)
            || (pole_distance(b) <= f64::EPSILON.sqrt() && x < -1e4)
            || (is_pole(a) && x.exp() == 0.0)
);
// known: NaN for an integer a at a tiny x (#2252)
continuous!(
    hypu(a in moderate(), b in moderate(), x) unless a.fract() == 0.0 && 0.0 < x && x < 1e-100
);

// known: AMOS is inaccurate for large x, e.g. bessel_j(-7, 2.71e79)
const MAX_BESSEL: f64 = (1 << 30) as f64;

// singular at the poles b = 0, -1, ...; known: beyond the Bessel cutoff at 2 sqrt(|x|) (#1472)
continuous!(hyp0f1(b in moderate(), x) unless is_pole(b) || 2.0 * x.abs().sqrt() > MAX_BESSEL);

continuous!(gamma(x));
continuous!(loggamma(x));
continuous!(rgamma(x));

continuous!(bessel_j(v in moderate(), x) unless x.abs() > MAX_BESSEL);
continuous!(bessel_y(v in moderate(), x) unless x.abs() > MAX_BESSEL);
continuous!(bessel_i(v in moderate(), x));
continuous!(bessel_k(v in moderate(), x));

continuous!(betainc(a in moderate(), b in moderate(), x in unit()));
continuous!(betaincinv(a in moderate(), b in moderate(), y in unit()));
continuous!(gammainc(a in moderate(), x));
continuous!(gammaincinv(a in moderate(), p in unit()));
continuous!(gdtrib(a in moderate(), p in unit(), x));
