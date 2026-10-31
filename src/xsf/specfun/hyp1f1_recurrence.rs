//! Unlike in xsf: the backward recurrences of Boost.Math for `M(a, b, z)` with non-integer
//! `a < -1` and `b > 0`, translated into pure Rust from
//! `boost/math/special_functions/detail/hypergeometric_1F1_recurrence.hpp` and
//! `boost/math/tools/recurrence.hpp` (Boost.Math develop @ cf949b2).
//!
//! Original header comments appear below.
//!
//! Copyright 2014 Anton Bikineev
//! Copyright 2014 Christopher Kormanyos
//! Copyright 2014 John Maddock
//! Copyright 2014 Paul Bristow
//! Distributed under the Boost
//! Software License, Version 1.0. (See accompanying file
//! LICENSE_1_0.txt or copy at <https://www.boost.org/LICENSE_1_0.txt>)
//!
//! (C) Copyright Anton Bikineev 2014
//! Use, modification and distribution are subject to the
//! Boost Software License, Version 1.0. (See accompanying file
//! LICENSE_1_0.txt or copy at <https://www.boost.org/LICENSE_1_0.txt>)

/// `policies::get_max_series_iterations<Policy>()` (the Boost default, also used by SciPy)
const MAX_SERIES_ITERATIONS: f64 = 1e6;

/// `max_scaling` in `hypergeometric_1F1_imp`, i.e. `lltrunc(log_max_value<T>()) - 2`
const MAX_SCALING: i64 = 707;

/// `hypergeometric_1F1_recurrence_a_coefficients`, i.e. DLMF 13.3.1:
/// `(b - a) M(a - 1, b, z) + (2a - b + z) M(a, b, z) - a M(a + 1, b, z) = 0`
fn recurrence_a_coefficients(a: f64, b: f64, z: f64) -> impl Fn(i32) -> (f64, f64, f64) {
    move |i| {
        let ai = a + f64::from(i);

        let an = b - ai;
        let bn = 2.0 * ai - b + z;
        let cn = -ai;

        (an, bn, cn)
    }
}

/// `hypergeometric_1F1_recurrence_small_b_coefficients`, i.e. DLMF 13.3.2:
/// `b (b - 1) M(a, b - 1, z) + b (1 - b - z) M(a, b, z) + z (b - a) M(a, b + 1, z) = 0`
//
// for use when we're recursing to a small b:
fn recurrence_small_b_coefficients(
    a: f64,
    b: f64,
    z: f64,
    n: i32,
) -> impl Fn(i32) -> (f64, f64, f64) {
    move |i| {
        let bi = b + f64::from(i + n);
        let bi_minus_1 = b + f64::from(i + n - 1);

        let an = bi * bi_minus_1;
        let bn = bi * (-bi_minus_1 - z);
        let cn = z * (bi - a);

        (an, bn, cn)
    }
}

/// `hypergeometric_1F1_recurrence_a_and_b_coefficients`, i.e.
/// `b (b - 1) M(a - 1, b - 1, z) + b (z - b + 1) M(a, b, z) - a z M(a + 1, b + 1, z) = 0`
fn recurrence_a_and_b_coefficients(
    a: f64,
    b: f64,
    z: f64,
    offset: i32,
) -> impl Fn(i32) -> (f64, f64, f64) {
    move |i| {
        let ai = a + f64::from(offset + i);
        let bi = b + f64::from(offset + i);

        let an = bi * (b + f64::from(offset + i - 1));
        let bn = bi * (z - (b + f64::from(offset + i - 1)));
        let cn = -ai * z;

        (an, bn, cn)
    }
}

/// `tools::apply_recurrence_relation_backward`, which solves the recurrence relation
/// `a(n) w(n-1) + b(n) w(n) + c(n) w(n+1) = 0` in the backward direction, starting from
/// `first = w(1)` and `second = w(0)`, and returns `(w(-number_of_steps), w(1-number_of_steps))`.
#[allow(clippy::cast_possible_truncation)]
fn apply_recurrence_relation_backward(
    get_coefs: impl Fn(i32) -> (f64, f64, f64),
    number_of_steps: i32,
    mut first: f64,
    mut second: f64,
    log_scaling: &mut i64,
) -> (f64, f64) {
    for k in 0..number_of_steps {
        let (a, b, c) = get_coefs(-k);

        if (second != 0.0)
            && (((f64::MAX * (a / b) / 2048.0).abs() < second.abs())
                || ((f64::MAX * (a / c) / 2048.0).abs() < first.abs())
                || ((f64::MIN_POSITIVE * (a / b) * 2048.0).abs() > second.abs())
                || ((f64::MIN_POSITIVE * (a / c) * 2048.0).abs() > first.abs()))
        {
            // Rescale everything:
            let log_scale = second.abs().ln() as i32;
            let scale = f64::from(-log_scale).exp();
            second *= scale;
            first *= scale;
            *log_scaling += i64::from(log_scale);
        }
        // scale each part separately to avoid spurious overflow:
        let next = (b / -a) * second + (c / -a) * first;

        first = second;
        second = next;
    }

    (second, first)
}

/// `hypergeometric_1F1_backward_recurrence_for_negative_a`: backward recurrence on `a` from
/// `M(frac(a) + 2, b, z)` and `M(frac(a) + 1, b, z)`, for non-integer `a < 0`
#[allow(clippy::cast_possible_truncation, clippy::float_cmp)]
fn backward_recurrence_for_negative_a(
    a: f64,
    b: f64,
    z: f64,
    m: impl Fn(f64, f64) -> Option<f64>,
    log_scaling: &mut i64,
) -> Option<f64> {
    let mut integer_part = a.trunc();
    let mut ak = a - integer_part;
    //
    // We need ak-1 positive to avoid infinite recursion below:
    //
    // NOTE: the `ak == 0` (integer `a`) branches are omitted, since `a` is never an integer here.
    ak += 2.0;
    integer_part -= 2.0;
    if ak - 1.0 == b {
        // When ak - 1 == b are recursion coefficients disappear to zero and
        // we end up with a NaN result.  Reduce the recursion steps by 1 to
        // avoid this.  We rely on |b| small and therefore no infinite recursion.
        ak -= 1.0;
        integer_part += 1.0;
    }

    // NOTE: Boost raises an evaluation error here.
    if -integer_part > MAX_SERIES_ITERATIONS {
        return None;
    }

    let first = m(ak, b)?;
    ak -= 1.0;
    let second = m(ak, b)?;
    integer_part += 1.0;

    let s = recurrence_a_coefficients(ak, b, z);

    let (second, _) =
        apply_recurrence_relation_backward(s, -integer_part as i32, first, second, log_scaling);
    Some(second)
}

/// `hypergeometric_1F1_backwards_recursion_on_b_for_negative_a`
#[allow(clippy::cast_possible_truncation)]
fn backwards_recursion_on_b_for_negative_a(
    a: f64,
    b: f64,
    z: f64,
    m: impl Fn(f64, f64) -> Option<f64>,
    log_scaling: &mut i64,
) -> Option<f64> {
    //
    // We compute
    //
    // M[a + a_shift, b + b_shift; z]
    //
    // and recurse backwards on a and b down to
    //
    // M[a, b, z]
    //
    // With a + a_shift > 1 and b + b_shift > z
    //
    // There are 3 distinct regions to ensure stability during the recursions:
    //
    // a > 0         :  stable for backwards on a
    // a < 0, b > 0  :  stable for backwards on a and b
    // a < 0, b < 0  :  stable for backwards on b (as long as |b| is small).
    //
    // We could simplify things by ignoring the middle region, but it's more efficient
    // to recurse on a and b together when we can.
    //

    debug_assert!(a < -1.0); // Not tested nor taken for -1 < a < 0

    let b_shift = (z - b).trunc() + 2.0;

    let mut a_shift = (-a).trunc();
    if a + a_shift != 0.0 {
        a_shift += 2.0;
    }
    //
    // If the shifts are so large that we would throw an evaluation_error, try the series instead,
    // even though this will almost certainly throw as well:
    //
    // NOTE: here `None` is returned instead.
    if b_shift > MAX_SERIES_ITERATIONS || a_shift > MAX_SERIES_ITERATIONS {
        return None;
    }
    let (b_shift, a_shift) = (b_shift as i32, a_shift as i32);

    // NOTE: for `b < 0` Boost uses `itrunc(b + b_shift)` instead of `b_shift`.
    let mut a_b_shift = b_shift; // The max we can shift on a and b together
    let mut leading_a_shift = 3.min(a_shift); // Just enough to make a negative
    if a_b_shift > a_shift - 3 {
        a_b_shift = if a_shift < 3 { 0 } else { a_shift - 3 };
    } else {
        // Need to ensure that leading_a_shift is large enough that a will reach it's target
        // after the first 2 phases (-,0) and (-,-) are over:
        leading_a_shift = a_shift - a_b_shift;
    }
    let mut trailing_b_shift = b_shift - a_b_shift;
    if a_b_shift < 5 {
        // Might as well do things in two steps rather than 3:
        if a_b_shift > 0 {
            leading_a_shift += a_b_shift;
            trailing_b_shift += a_b_shift;
        }
        a_b_shift = 0;
        leading_a_shift -= 1;
    }

    debug_assert!(leading_a_shift > 1);
    debug_assert_eq!(
        a_b_shift + leading_a_shift + i32::from(a_b_shift == 0),
        a_shift
    );
    debug_assert_eq!(a_b_shift + trailing_b_shift, b_shift);

    if (trailing_b_shift == 0) && (b.abs() < 0.5) && a_b_shift != 0 {
        // Better to have the final recursion on b alone, otherwise we lose precision when b is very small:
        let diff = a_b_shift.min(3);
        a_b_shift -= diff;
        leading_a_shift += diff;
        trailing_b_shift += diff;
    }

    let mut first = m(a + f64::from(a_shift), b + f64::from(b_shift))?;
    //
    // It would be good to compute "second" from first and the ratio - unfortunately we are right on the cusp
    // recursion on a switching from stable backwards to stable forwards behaviour and so this is not possible here.
    //
    let mut second = m(a + f64::from(a_shift - 1), b + f64::from(b_shift))?;

    //
    // Now we have [a + a_shift, b + b_shift, z] and [a + a_shift - 1, b + b_shift, z]
    // and want to recurse until [a + a_shift - leading_a_shift, b + b_shift, z] and [a + a_shift - leadng_a_shift - 1, b + b_shift, z]
    // which is leading_a_shift -1 steps.
    //
    // NOTE: it is actually `leading_a_shift` steps.
    (second, first) = apply_recurrence_relation_backward(
        recurrence_a_coefficients(a + f64::from(a_shift - 1), b + f64::from(b_shift), z),
        leading_a_shift,
        first,
        second,
        log_scaling,
    );

    if a_b_shift != 0 {
        //
        // Now we need to switch to an a+b shift so that we have:
        // [a + a_shift - leading_a_shift, b + b_shift, z] and [a + a_shift - leadng_a_shift - 1, b + b_shift - 1, z]
        // A&S 13.4.3 gives us what we need:
        //
        {
            // local a's and b's:
            let la = a + f64::from(a_shift - leading_a_shift - 1);
            let lb = b + f64::from(b_shift);
            second = ((1.0 + la - lb) * second - la * first) / (1.0 - lb);
        }
        //
        // Now apply a_b_shift - 1 recursions to get down to
        // [a + 1, b + trailing_b_shift + 1, z] and [a, b + trailing_b_shift, z]
        //
        (second, first) = apply_recurrence_relation_backward(
            recurrence_a_and_b_coefficients(
                a,
                b + f64::from(b_shift - a_b_shift),
                z,
                a_b_shift - 1,
            ),
            a_b_shift - 1,
            first,
            second,
            log_scaling,
        );
        //
        // Now we need to switch to a b shift, a different application of A&S 13.4.3
        // will get us there, we leave "second" where it is, and move "first" sideways:
        //
        {
            let lb = b + f64::from(trailing_b_shift + 1);
            first = (second * (lb - 1.0) - a * first) / -(1.0 + a - lb);
        }
    } else {
        //
        // We have M[a+1, b+b_shift, z] and M[a, b+b_shift, z] and need M[a, b+b_shift-1, z] for
        // recursion on b: A&S 13.4.3 gives us what we need.
        //
        let third = -(second * (1.0 + a - b - f64::from(b_shift)) - first * a)
            / (b + f64::from(b_shift - 1));
        first = second;
        second = third;
        trailing_b_shift -= 1;
    }
    //
    // Finish off by applying trailing_b_shift recursions:
    //
    if trailing_b_shift != 0 {
        (second, _) = apply_recurrence_relation_backward(
            recurrence_small_b_coefficients(a, b, z, trailing_b_shift),
            trailing_b_shift,
            first,
            second,
            log_scaling,
        );
    }
    Some(second)
}

/// `e^s M(a, b, z)`, from the `b > 0` branch of `hypergeometric_1F1_divergent_fallback`.
///
/// `m(a', b')` computes the initial values `M(a', b', z)`, with `-1 < a' < 2` and `b' > z`.
/// Returns `None` if `m` does, or if more than 10^6 steps are needed.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
pub(super) fn hyp1f1_neg_a(
    a: f64,
    b: f64,
    z: f64,
    s: f64,
    m: impl Fn(f64, f64) -> Option<f64>,
) -> Option<f64> {
    let mut log_scaling = 0;
    let mut result = if z < b {
        backward_recurrence_for_negative_a(a, b, z, m, &mut log_scaling)
    } else {
        backwards_recursion_on_b_for_negative_a(a, b, z, m, &mut log_scaling)
    }?;

    // `e^s` is applied like `e^z` for `z < 0` in `hypergeometric_1F1_imp`:
    let scaling = s.trunc();
    result *= (s - scaling).exp();
    let mut log_scaling = log_scaling.saturating_add(scaling as i64);

    // The rest is from `hypergeometric_1F1_imp`:
    //
    // Actual result will be result * e^log_scaling.
    //
    let max_scale_factor = (MAX_SCALING as f64).exp();

    // NOTE: clamped to bound the loops; beyond that the result over- or underflows anyway,
    // since `ln|result|` is in [-745, 710).
    log_scaling = log_scaling.clamp(-3 * MAX_SCALING, 3 * MAX_SCALING);
    while !result.is_infinite() && (log_scaling > MAX_SCALING) {
        result *= max_scale_factor;
        log_scaling -= MAX_SCALING;
    }
    while log_scaling < -MAX_SCALING {
        result /= max_scale_factor;
        log_scaling += MAX_SCALING;
    }
    if log_scaling != 0 {
        result *= (log_scaling as f64).exp();
    }
    Some(result)
}
