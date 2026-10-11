//! Translated into pure Rust from `xsf/cephes/owens_t.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024. Original header comment appears below.
//!
//! Copyright Benjamin Sobotta 2012
//!
//! Use, modification and distribution are subject to the
//! Boost Software License, Version 1.0. (See accompanying file
//! LICENSE_1_0.txt or copy at <https://www.boost.org/LICENSE_1_0.txt>)
//!
//! Reference:
//! Mike Patefield, David Tandy
//! FAST AND ACCURATE CALCULATION OF OWEN'S T-FUNCTION
//! Journal of Statistical Software, 5 (5), 1-25
//!
//! The Cephes `expm1` from `cephes/unity.h` is replaced by `f64::exp_m1` from the Rust standard
//! library, like in `xsf/exp.rs`.

use core::f64::consts::{SQRT_2, TAU};

use crate::xsf::cephes::{erf, erfc, ndtr};

const OWENS_T_SELECT_METHOD: [usize; 120] = [
    0, 0, 1, 12, 12, 12, 12, 12, 12, 12, 12, 15, 15, 15, 8, 0, 1, 1, 2, 2, 4, 4, 13, 13, 14, 14,
    15, 15, 15, 8, 1, 1, 2, 2, 2, 4, 4, 14, 14, 14, 14, 15, 15, 15, 9, 1, 1, 2, 4, 4, 4, 4, 6, 6,
    15, 15, 15, 15, 15, 9, 1, 2, 2, 4, 4, 5, 5, 7, 7, 16, 16, 16, 11, 11, 10, 1, 2, 4, 4, 4, 5, 5,
    7, 7, 16, 16, 16, 11, 11, 11, 1, 2, 3, 3, 5, 5, 7, 7, 16, 16, 16, 16, 16, 11, 11, 1, 2, 3, 3,
    5, 5, 17, 17, 17, 17, 16, 16, 16, 11, 11,
];

const OWENS_T_METHODS: [i32; 18] = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6];

const OWENS_T_HRANGE: [f64; 14] = [
    0.02, 0.06, 0.09, 0.125, 0.26, 0.4, 0.6, 1.6, 1.7, 2.33, 2.4, 3.36, 3.4, 4.8,
];

const OWENS_T_ARANGE: [f64; 7] = [0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.999_99];

const OWENS_T_ORD: [f64; 18] = [
    2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 12.0, 18.0, 10.0, 20.0, 30.0, 0.0, 4.0, 7.0, 8.0, 20.0, 0.0, 0.0,
];

const OWENS_T_C: [f64; 31] = [
    1.0,
    -1.0,
    1.0,
    -0.999_999_999_999_999_8,
    0.999_999_999_999_983_9,
    -0.999_999_999_999_306_3,
    0.999_999_999_979_733_7,
    -0.999_999_999_574_958_4,
    0.999_999_993_322_623_5,
    -0.999_999_918_892_324_2,
    0.999_999_219_514_348_3,
    -0.999_993_935_137_206_7,
    0.999_961_355_976_905_5,
    -0.999_795_563_665_139_4,
    0.999_092_789_629_617_1,
    -0.996_593_837_411_918_2,
    0.989_100_171_383_861_3,
    -0.970_078_558_040_693_3,
    0.929_114_386_832_631_9,
    -0.854_205_869_595_615_6,
    0.737_965_260_330_301,
    -0.585_234_698_828_374,
    0.415_997_776_145_676_3,
    -0.258_821_087_524_194_36,
    0.137_553_582_516_389_27,
    -0.060_795_276_632_595_575,
    0.021_633_768_329_987_153,
    -0.005_934_056_934_551_867,
    0.001_174_341_481_833_294_6,
    -0.000_148_915_561_335_036_9,
    9.072_354_320_794_358e-6,
];

const OWENS_T_PTS: [f64; 13] = [
    0.003_508_203_967_645_171_6,
    0.031_279_042_338_030_756,
    0.085_266_826_283_219_45,
    0.162_450_717_308_122_77,
    0.258_511_960_491_254_36,
    0.368_075_538_406_975_3,
    0.485_010_929_056_047,
    0.602_775_141_526_185_7,
    0.714_778_842_177_532_3,
    0.814_755_109_887_601,
    0.897_110_297_559_489_7,
    0.957_238_080_859_442_6,
    0.991_788_329_746_297,
];

const OWENS_T_WTS: [f64; 13] = [
    0.018_831_438_115_323_503,
    0.018_567_086_243_977_65,
    0.018_042_093_461_223_385,
    0.017_263_829_606_398_752,
    0.016_243_219_975_989_858,
    0.014_994_592_034_116_705,
    0.013_535_474_469_662_09,
    0.011_886_351_605_820_165,
    0.010_070_377_242_777_432,
    0.008_113_054_574_229_958,
    0.006_041_900_952_847_024,
    0.003_886_221_701_074_205_7,
    0.001_679_303_108_454_609,
];

fn get_method(h: f64, a: f64) -> usize {
    let mut ihint = 14;
    let mut iaint = 7;

    for (i, &hr) in OWENS_T_HRANGE.iter().enumerate() {
        if h <= hr {
            ihint = i;
            break;
        }
    }

    for (i, &ar) in OWENS_T_ARANGE.iter().enumerate() {
        if a <= ar {
            iaint = i;
            break;
        }
    }
    OWENS_T_SELECT_METHOD[iaint * 15 + ihint]
}

fn owens_t_norm1(x: f64) -> f64 {
    erf(x / SQRT_2) / 2.0
}

fn owens_t_norm2(x: f64) -> f64 {
    erfc(x / SQRT_2) / 2.0
}

fn owens_t1(h: f64, a: f64, m: f64) -> f64 {
    let mut j = 1;
    let mut jj = 1;

    let hs = -0.5 * h * h;
    let dhs = hs.exp();
    let as_ = a * a;
    let mut aj = a / TAU;
    let mut dj = hs.exp_m1();
    let mut gj = hs * dhs;

    let mut val = a.atan() / TAU;

    loop {
        val += dj * aj / f64::from(jj);

        if m <= f64::from(j) {
            break;
        }
        j += 1;
        jj += 2;
        aj *= as_;
        dj = gj - dj;
        gj *= hs / f64::from(j);
    }

    val
}

#[allow(clippy::cast_possible_truncation)]
fn owens_t2(h: f64, a: f64, ah: f64, m: f64) -> f64 {
    let mut i = 1;
    let maxi = (2.0 * m + 1.0) as i32;
    let hs = h * h;
    let as_ = -a * a;
    let y = 1.0 / hs;
    let mut val = 0.0;
    let mut vi = a * (-0.5 * ah * ah).exp() / TAU.sqrt();
    let mut z = (ndtr(ah) - 0.5) / h;

    loop {
        val += z;
        if maxi <= i {
            break;
        }
        z = y * (vi - f64::from(i) * z);
        vi *= as_;
        i += 2;
    }
    val *= (-0.5 * hs).exp() / TAU.sqrt();

    val
}

fn owens_t3(h: f64, a: f64, ah: f64) -> f64 {
    let aa = a * a;
    let hh = h * h;
    let y = 1.0 / hh;

    let mut vi = a * (-ah * ah / 2.0).exp() / TAU.sqrt();
    let mut zi = owens_t_norm1(ah) / h;
    let mut result = 0.0;

    for (i, &c) in (0..=30).zip(&OWENS_T_C) {
        result += zi * c;
        zi = y * (f64::from(2 * i + 1) * zi - vi);
        vi *= aa;
    }

    result *= (-hh / 2.0).exp() / TAU.sqrt();

    result
}

fn owens_t4(h: f64, a: f64, m: f64) -> f64 {
    let maxi = 2.0 * m + 1.0;
    let hh = h * h;
    let naa = -a * a;

    let mut i = 1;
    let mut ai = a * (-hh * (1.0 - naa) / 2.0).exp() / TAU;
    let mut yi = 1.0;
    let mut result = 0.0;

    loop {
        result += ai * yi;

        if maxi <= f64::from(i) {
            break;
        }

        i += 2;
        yi = (1.0 - hh * yi) / f64::from(i);
        ai *= naa;
    }

    result
}

fn owens_t5(h: f64, a: f64) -> f64 {
    let mut result = 0.0;
    let aa = a * a;
    let nhh = -0.5 * h * h;

    for (pt, wt) in OWENS_T_PTS.iter().zip(&OWENS_T_WTS) {
        let r = 1.0 + aa * pt;
        result += wt * (nhh * r).exp() / r;
    }

    result *= a;

    result
}

#[allow(clippy::float_cmp)]
fn owens_t6(h: f64, a: f64) -> f64 {
    let normh = owens_t_norm2(h);
    let y = 1.0 - a;
    let r = y.atan2(1.0 + a);
    let mut result = normh * (1.0 - normh) / 2.0;

    if r != 0.0 {
        result -= r * (-y * h * h / (2.0 * r)).exp() / TAU;
    }

    result
}

#[allow(clippy::float_cmp)]
fn owens_t_dispatch(h: f64, a: f64, ah: f64) -> f64 {
    if h == 0.0 {
        return a.atan() / TAU;
    }
    if a == 0.0 {
        return 0.0;
    }
    if a == 1.0 {
        return owens_t_norm2(-h) * owens_t_norm2(h) / 2.0;
    }

    let index = get_method(h, a);
    let m = OWENS_T_ORD[index];
    let meth_code = OWENS_T_METHODS[index];

    match meth_code {
        1 => owens_t1(h, a, m),
        2 => owens_t2(h, a, ah, m),
        3 => owens_t3(h, a, ah),
        4 => owens_t4(h, a, m),
        5 => owens_t5(h, a),
        6 => owens_t6(h, a),
        _ => f64::NAN,
    }
}

pub(crate) fn owens_t(mut h: f64, a: f64) -> f64 {
    if h.is_nan() || a.is_nan() {
        return f64::NAN;
    }

    /* exploit that T(-h,a) == T(h,a) */
    h = h.abs();

    /*
     * Use equation (2) in the paper to remap the arguments such that
     * h >= 0 and 0 <= a <= 1 for the call of the actual computation
     * routine.
     */
    let fabs_a = a.abs();
    let fabs_ah = fabs_a * h;

    let result = if fabs_a == f64::INFINITY {
        /* See page 13 in the paper */
        0.5 * owens_t_norm2(h)
    } else if h == f64::INFINITY {
        0.0
    } else if fabs_a <= 1.0 {
        owens_t_dispatch(h, fabs_a, fabs_ah)
    } else if fabs_ah <= 0.67 {
        let normh = owens_t_norm1(h);
        let normah = owens_t_norm1(fabs_ah);
        0.25 - normh * normah - owens_t_dispatch(fabs_ah, 1.0 / fabs_a, h)
    } else {
        let normh = owens_t_norm2(h);
        let normah = owens_t_norm2(fabs_ah);
        f64::midpoint(normh, normah) - normh * normah - owens_t_dispatch(fabs_ah, 1.0 / fabs_a, h)
    };

    if a < 0.0 {
        /* exploit that T(h,-a) == -T(h,a) */
        return -result;
    }

    result
}
