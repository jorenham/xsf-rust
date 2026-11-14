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
    let mut resasc = WGK[10] * (fc - reskh).abs();
    for (j, (fval1, fval2)) in fv.into_iter().enumerate() {
        resasc += WGK[j] * ((fval1 - reskh).abs() + (fval2 - reskh).abs());
    }

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

        let (i, _) = intervals
            .iter()
            .enumerate()
            .max_by(|(_, x), (_, y)| x.2.1.total_cmp(&y.2.1))
            .unwrap();
        let (a1, b2, _) = intervals.swap_remove(i);
        let b1 = 0.5 * (a1 + b2);
        intervals.push((a1, b1, qk21(&f, a1, b1)));
        intervals.push((b1, b2, qk21(&f, b1, b2)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    #[test]
    fn test_quad_smooth() {
        let (result, abserr) = quad(f64::sin, 0.0, PI, 0.0, 1e-13);
        assert!(
            (result - 2.0).abs() <= 4e-16 && abserr <= 2e-13,
            "{result} {abserr}"
        );

        // reversed, and with a near-singularity at x = -0.001
        let (result, abserr) = quad(|x| 1.0 / (x + 0.001), 1.0, 0.0, 0.0, 1e-13);
        let expected = -1001.0_f64.ln();
        assert!(
            (result - expected).abs() <= 1e-15 && abserr <= 7e-13,
            "{result} {abserr}"
        );
    }

    #[test]
    fn test_quad_limit() {
        // the endpoint singularity converges too slowly for 50 subintervals
        let (result, abserr) = quad(|x| 1.0 / x.sqrt(), 0.0, 1.0, 0.0, 1e-14);
        assert!(
            (result - 2.0).abs() <= abserr && abserr > 1e-14 * 2.0,
            "{result} {abserr}"
        );
    }

    #[test]
    fn test_quad_nan() {
        let (result, abserr) = quad(|x| (x - 0.5).sqrt(), 0.0, 1.0, 0.0, 1e-13);
        assert!(result.is_nan() && abserr.is_nan());
    }
}
