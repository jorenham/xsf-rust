//! Adaptive Gauss-Kronrod quadrature, after QUADPACK's `dqag` with the 21-point rule `dqk21`

/// Maximum number of subintervals, like `scipy.integrate.quad`
const LIMIT: usize = 50;

/// Abscissae of the 21-point Kronrod rule; the odd ones (from 0) are those of the 10-point Gauss
/// rule
#[rustfmt::skip]
const XGK: [f64; 11] = [
    0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
    0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
    0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
    0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
    0.294392862701460198131126603103866, 0.148874338981631210884826001129720,
    0.000000000000000000000000000000000,
];

/// Weights of the 21-point Kronrod rule
#[rustfmt::skip]
const WGK: [f64; 11] = [
    0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821,
];

/// Weights of the 10-point Gauss rule
#[rustfmt::skip]
const WG: [f64; 5] = [
    0.066671344308688137593568809893332, 0.149451349150580593145776339657697,
    0.219086362515982043995534934228163, 0.269266719309996355091226921569469,
    0.295524224714752870173892994651338,
];

/// `dqk21`: the 21-point Kronrod estimate of the integral of `f` over `[a, b]`, and its error
fn qk21(f: &impl Fn(f64) -> f64, a: f64, b: f64) -> (f64, f64) {
    let centr = 0.5 * (a + b);
    let hlgth = 0.5 * (b - a);

    let fc = f(centr);
    let mut resg = 0.0;
    let mut resk = WGK[10] * fc;
    let mut resabs = resk.abs();
    let mut fv = [(0.0, 0.0); 10];
    for (j, fvj) in fv.iter_mut().enumerate() {
        let absc = hlgth * XGK[j];
        let (fval1, fval2) = (f(centr - absc), f(centr + absc));
        *fvj = (fval1, fval2);
        let fsum = fval1 + fval2;
        if j % 2 == 1 {
            resg += WG[j / 2] * fsum;
        }
        resk += WGK[j] * fsum;
        resabs += WGK[j] * (fval1.abs() + fval2.abs());
    }

    let reskh = resk * 0.5;
    let resasc = fv.into_iter().enumerate().fold(
        WGK[10] * (fc - reskh).abs(),
        |resasc, (j, (fval1, fval2))| {
            resasc + WGK[j] * ((fval1 - reskh).abs() + (fval2 - reskh).abs())
        },
    );

    let result = resk * hlgth;
    let resabs = resabs * hlgth.abs();
    let resasc = resasc * hlgth.abs();
    let mut abserr = ((resk - resg) * hlgth).abs();
    if resasc != 0.0 && abserr != 0.0 {
        abserr = resasc * (200.0 * abserr / resasc).powf(1.5).min(1.0);
    }
    if resabs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) {
        abserr = abserr.max(50.0 * f64::EPSILON * resabs);
    }
    (result, abserr)
}

/// The integral of `f` over `[a, b]`, and an estimate of its absolute error
///
/// The subinterval with the largest error is bisected, until the error is at most
/// `max(epsabs, epsrel * |result|)`, or there are `LIMIT` subintervals. Unlike `dqag`, this doesn't
/// stop early on roundoff, and it sums the subintervals afresh in each step.
///
/// The error estimate is at least `50 * f64::EPSILON * integral(|f|)`, so a smaller `epsrel` only
/// uses up all subintervals.
pub(crate) fn quad(f: impl Fn(f64) -> f64, a: f64, b: f64, epsabs: f64, epsrel: f64) -> (f64, f64) {
    // `(a, b, result, abserr)` of each subinterval
    let mut intervals = vec![(a, b, qk21(&f, a, b))];
    loop {
        let result: f64 = intervals.iter().map(|&(_, _, (r, _))| r).sum();
        let abserr: f64 = intervals.iter().map(|&(_, _, (_, e))| e).sum();
        if !(abserr > epsabs.max(epsrel * result.abs())) || intervals.len() == LIMIT {
            return (result, abserr);
        }

        let (i, _) =
            intervals.iter().enumerate().max_by(|(_, x), (_, y)| x.2.1.total_cmp(&y.2.1)).unwrap();
        let (a1, b2, _) = intervals.swap_remove(i);
        let b1 = 0.5 * (a1 + b2);
        intervals.push((a1, b1, qk21(&f, a1, b1)));
        intervals.push((b1, b2, qk21(&f, b1, b2)));
    }
}

/// The integral of a unimodal `f >= 0` over `[a, b]`, with its maximum near `mode`
///
/// It sums [`quad`] with `epsrel` over pieces of width `width`, `2 width`, `4 width`, ... on either
/// side of `mode`, until a piece is no larger than the one before it and at most `f64::EPSILON`
/// times the sum, so `mode` and `width` need not be exact. That assumes the rest of that side is no
/// larger than the last piece, e.g. if `f` decays at least as fast as `|x - mode|^-2`. `a` and `b`
/// can be infinite.
pub(crate) fn quad_unimodal(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    mode: f64,
    width: f64,
    epsrel: f64,
) -> f64 {
    // at least 8 ulps of `mode`, so that the rounding of `x +- w` can't undo its doubling
    let width = width.max(8.0 * f64::EPSILON * mode.abs());
    let mut sum = 0.0;
    for end in [b, a] {
        let (mut x, mut w, mut prev) = (mode, width, 0.0);
        while x != end {
            let y = if end > x { (x + w).min(end) } else { (x - w).max(end) };
            let (piece, _) = quad(&f, x.min(y), x.max(y), epsrel * sum, epsrel);
            sum += piece;
            if !(piece > prev || piece > f64::EPSILON * sum) {
                break;
            }
            (x, w, prev) = (y, 2.0 * w, piece);
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::SQRT_TAU;
    use core::f64::consts::PI;

    #[test]
    fn test_quad_smooth() {
        let (result, abserr) = quad(f64::sin, 0.0, PI, 0.0, 1e-13);
        assert!((result - 2.0).abs() <= 4e-16 && abserr <= 2e-13, "{result} {abserr}");

        // reversed, and with a near-singularity at x = -0.001
        let (result, abserr) = quad(|x| 1.0 / (x + 0.001), 1.0, 0.0, 0.0, 1e-13);
        let expected = -1001.0_f64.ln();
        assert!((result - expected).abs() <= 1e-15 && abserr <= 7e-13, "{result} {abserr}");
    }

    #[test]
    fn test_quad_limit() {
        // the endpoint singularity converges too slowly for 50 subintervals
        let (result, abserr) = quad(|x| 1.0 / x.sqrt(), 0.0, 1.0, 0.0, 1e-14);
        assert!((result - 2.0).abs() <= abserr && abserr > 1e-14 * 2.0, "{result} {abserr}");
    }

    #[test]
    fn test_quad_unimodal() {
        // a Gaussian with sigma = 0.01, from a mode that is off
        let f = |x: f64| (-0.5 * x * x / 1e-4).exp();
        let expected = 0.01 * SQRT_TAU;
        let result = quad_unimodal(f, f64::NEG_INFINITY, f64::INFINITY, 0.03, 0.01, 2e-14);
        assert!((result - expected).abs() <= 1e-15 * expected, "{result} {expected}");

        // the exponential decay on [0, inf), which has its maximum at the end
        let result = quad_unimodal(|x: f64| (-x).exp(), 0.0, f64::INFINITY, 0.0, 1e-3, 2e-14);
        assert!((result - 1.0).abs() <= 1e-15, "{result}");

        // a width far below the scale, so the pieces grow on both sides
        let f = |x: f64| (-0.5 * x * x).exp();
        let result = quad_unimodal(f, f64::NEG_INFINITY, f64::INFINITY, 0.0, 1e-18, 1e-13);
        assert!((result - SQRT_TAU).abs() <= 1e-15, "{result}");

        // a width below the ulp of the mode
        let f = |x: f64| (-0.5 * (x - 1.0) * (x - 1.0)).exp();
        let result = quad_unimodal(f, f64::NEG_INFINITY, f64::INFINITY, 1.0, 1e-20, 1e-13);
        assert!((result - SQRT_TAU).abs() <= 1e-15, "{result}");
    }

    #[test]
    fn test_quad_nan() {
        let (result, abserr) = quad(|x| (x - 0.5).sqrt(), 0.0, 1.0, 0.0, 1e-13);
        assert!(result.is_nan() && abserr.is_nan());
    }
}
