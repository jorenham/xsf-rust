//! Translated into pure Rust from `xsf/cdflib.h` (xsf v0.2.2).
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, the comment on the `increasing` argument is corrected ("not increasing" instead of
//! "not decreasing").
//!
//! Note that our `igam` and `igamc` can differ from those in xsf in the last bits (e.g. because
//! `f64::exp_m1` is used instead of the Cephes `expm1`), so the roots found by `gdtrib` can as well.

use crate::xsf::cephes::{igam, igamc};
use crate::xsf::tools::{bracket_root_for_cdf_inversion, find_root_chandrupatla};

/// Inverse of [`gdtr`](crate::gdtr) with respect to `b`
///
/// Corresponds to [`scipy.special.gdtrib`][gdtrib].
///
/// [gdtrib]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gdtrib.html
///
/// # See also
/// - [`gdtr`](crate::gdtr): Gamma CDF
/// - [`gdtrc`](crate::gdtrc): Gamma survival function
#[must_use]
#[allow(
    clippy::float_cmp,
    clippy::manual_range_contains,
    clippy::neg_cmp_op_on_partial_ord
)]
pub fn gdtrib(a: f64, p: f64, x: f64) -> f64 {
    if p.is_nan() || a.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if !((0.0 <= p) && (p <= 1.0)) {
        return f64::NAN;
    }
    if !(a > 0.0) || a.is_infinite() {
        return f64::NAN;
    }
    if !(x >= 0.0) || x.is_infinite() {
        return f64::NAN;
    }
    if x == 0.0 {
        if p == 0.0 {
            return f64::NAN;
        }
        /* gdtrib(a, p, x) tends to 0 as x -> 0 when p > 0 */
        return 0.0;
    }
    if p == 0.0 {
        /* gdtrib(a, p, x) tends to infinity as p -> 0 from the right when x > 0. */
        return f64::INFINITY;
    }
    if p == 1.0 {
        /* gdtrib(a, p, x) tends to 0 as p -> 1.0 from the left when x > 0. */
        return 0.0;
    }
    let q = 1.0 - p;
    let func = |b: f64| {
        if p <= q {
            return igam(b, a * x) - p;
        }
        q - igamc(b, a * x)
    };
    let lower_bound = f64::MIN_POSITIVE;
    let upper_bound = f64::MAX;
    /* To explain the magic constants used below:
     * 1.0 is the initial guess for the root. -0.875 is the initial step size
     * for the leading bracket endpoint if the bracket search will proceed to the
     * left, likewise 7.0 is the initial step size when the bracket search will
     * proceed to the right. 0.125 is the scale factor for a left moving bracket
     * search and 8.0 the scale factor for a right moving bracket search. These
     * constants are chosen so that:
     *
     * 1. The scale factor and bracket endpoints remain powers of 2, allowing for
     *    exact arithmetic, preventing roundoff error from causing numerical catastrophe
     *    which could lead to unexpected results.
     * 2. The bracket sizes remain constant in a relative sense. Each candidate bracket
     *    will contain roughly the same number of floating point values. This means that
     *    the number of necessary function evaluations in the worst case scenario for
     *    Chandrupatla's algorithm will remain constant.
     *
     * false specifies that the function is not increasing. 342 is equal to
     * max(ceil(log_8(DBL_MAX)), ceil(log_(1/8)(DBL_MIN))). An upper bound for the
     * number of iterations needed in this bracket search to check all normalized
     * floating point values.
     */
    let (xl, xr, f_xl, f_xr, bracket_status) = bracket_root_for_cdf_inversion(
        func,
        1.0,
        lower_bound,
        upper_bound,
        -0.875,
        7.0,
        0.125,
        8.0,
        false,
        342,
    );
    if bracket_status == 1 {
        return 0.0;
    }
    if bracket_status == 2 {
        return f64::INFINITY;
    }
    if bracket_status >= 3 {
        return f64::NAN;
    }
    let (result, root_status) =
        find_root_chandrupatla(func, xl, xr, f_xl, f_xr, f64::EPSILON, 1e-100, 100);
    if root_status != 0 {
        /* The root finding return should only fail if there's a bug in our code. */
        return f64::NAN;
    }
    result
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_gdtrib() {
        xsref::test("gdtrib", "d_d_d-d", |x| crate::gdtrib(x[0], x[1], x[2]));
    }
}
