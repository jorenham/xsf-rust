//! Translated into pure Rust from `xsf/ellint_carlson/ellint_common.hh` (xsf main).

/* Number of correction terms in polynomial evaluation. */
pub(super) const NPOLY: usize = 4;

/* Configuration constants */
pub(super) mod config {
    pub(crate) const MAX_ITER: u32 = 1000;
    pub(crate) const ASYM_ZERO_UL: f64 = 5.0e-14;
    pub(crate) const ASYM_CLOSE_UL: f64 = 1.0e-9;
}

pub(super) mod constants {
    pub(crate) use core::f64::consts::PI;

    pub(crate) const RC_C: [f64; 8] = [
        80_080.0, //
        0.0,      //
        24_024.0, /* 3 / 10 */
        11_440.0, /* 1 / 7 */
        30_030.0, /* 3 / 8 */
        32_760.0, /* 9 / 22 */
        61_215.0, /* 159 / 208 */
        90_090.0, /* 9 / 8 */
    ];

    pub(crate) const RF_C1: [f64; 4] = [0.0, -24_024.0, 10_010.0, -5775.0];

    pub(crate) const RF_C2: [f64; 3] = [17_160.0, -16_380.0, 15_015.0];

    pub(crate) const RF_C33: f64 = 6930.0; /*  3 / 104 */

    pub(crate) const RF_DENOM: f64 = 240_240.0;

    pub(crate) const RDJ_C1: [f64; 4] = [0.0, -875_160.0, 417_690.0, -255_255.0];

    pub(crate) const RDJ_C2: [f64; 3] = [0.0, 680_680.0, 306_306.0];

    pub(crate) const RDJ_C3: [f64; 3] = [0.0, -706_860.0, 675_675.0];

    pub(crate) const RDJ_C4: [f64; 2] = [-556_920.0, 612_612.0];

    pub(crate) const RDJ_C5: [f64; 2] = [471_240.0, -540_540.0];

    pub(crate) const RDJ_DENOM: f64 = 4_084_080.0;
}

/* Exit statuses for the library functions. */
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)] // the statuses are discarded, since `set_error` is a no-op
pub(super) enum ExitStatus {
    Success = 0,
    Singular,
    Underflow,
    Overflow,
    NIter,
    PrecLoss,
    NoResult,
    BadArgs,
    BadRerr,
    Other,
    Memory,
    Unused,
}

pub(super) const fn is_horrible(n: ExitStatus) -> bool {
    matches!(
        n,
        ExitStatus::NoResult | ExitStatus::BadArgs | ExitStatus::BadRerr | ExitStatus::Other
    )
}

pub(super) const fn is_troublesome(n: ExitStatus) -> bool {
    !(matches!(n, ExitStatus::Success) || is_horrible(n))
}

/* Comparison function based on absolute value. */
#[inline]
pub(super) fn abscmp<T: super::typing::Num>(a: &T, b: &T) -> bool {
    a.abs() < b.abs()
}

/// `std::sort(first, last, comp)` of libstdc++ for at most 16 elements, i.e. an insertion sort
///
/// Unlike Rust's `sort_by`, this only requires a "less than" predicate, so that the order of
/// incomparable (e.g. NaN) elements is the same as in C++.
pub(super) fn sort_by<T: Copy>(arr: &mut [T], comp: impl Fn(&T, &T) -> bool) {
    debug_assert!(arr.len() <= 16);
    for i in 1..arr.len() {
        let val = arr[i];
        if comp(&val, &arr[0]) {
            arr.copy_within(0..i, 1);
            arr[0] = val;
        } else {
            let mut j = i;
            while comp(&val, &arr[j - 1]) {
                arr[j] = arr[j - 1];
                j -= 1;
            }
            arr[j] = val;
        }
    }
}

/// Unlike in xsf, `rd`, `rg` and `rj` rescale their arguments by `f^2` if the largest part is
/// outside [2^-256, 2^256], since e.g. `elliprg(5.17e139, 2.97e11, 1.55e225)` overflowed (NaN,
/// not 1.97e112). `f` is a power of 2, and doesn't downscale the smallest nonzero part below `min`.
#[allow(clippy::float_cmp)]
pub(super) fn scale_factor<T: super::typing::Num>(args: &[T], min: f64) -> Option<f64> {
    let parts = args.iter().flat_map(|a| [a.real().abs(), a.imag().abs()]);
    let amax = parts.clone().fold(0.0, f64::max);
    let amin = parts.filter(|&a| a > 0.0).fold(f64::INFINITY, f64::min);
    let e = amax.log2();
    if !(e.is_finite() && e.abs() > 256.0) {
        return None;
    }
    let mut f = ((254.0 - e) / 2.0).floor().exp2();
    while f < 1.0 && amin * f * f < min {
        f *= 2.0;
    }
    (f != 1.0).then_some(f)
}

/// `std::max({a, ...})` (with `std::max(a, b)` = `a < b ? b : a`), i.e. the first of the largest
/// values, ignoring NaN unless it is the first one
pub(super) fn max(vals: &[f64]) -> f64 {
    let mut largest = vals[0];
    for &v in &vals[1..] {
        if largest < v {
            largest = v;
        }
    }
    largest
}
