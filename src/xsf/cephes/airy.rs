//! Translated into pure Rust from `xsf/cephes/airy.h` (xsf v0.2.2), which was translated into
//! C++ by SciPy developers in 2024.
//!
//! Cephes Math Library Release 2.8:  June, 2000
//! Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
//!
//! `set_error` is a no-op in our build of xsf, so the calls to it are omitted.
//!
//! Unlike in xsf, `(status, ai, aip, bi, bip)` is returned instead of written to output pointers.

#![allow(clippy::unreadable_literal, clippy::excessive_precision)]

use crate::xsf::cephes::polevl::{p1evl, polevl};
use core::f64::consts::PI;

const MACHEP: f64 = 1.110_223_024_625_156_5e-16; // 2**-53
const SQRT1OPI: f64 = 5.641895835477562869480794515607725858441E-1; // sqrt(1/pi)
const SQRT3: f64 = 1.732050807568877293527446341505872366943; // sqrt(3)

const C1: f64 = 0.35502805388781723926;
const C2: f64 = 0.258819403792806798405;
const MAXAIRY: f64 = 103.892;

const AN: [f64; 8] = [
    3.46538101525629032477E-1,
    1.20075952739645805542E1,
    7.62796053615234516538E1,
    1.68089224934630576269E2,
    1.59756391350164413639E2,
    7.05360906840444183113E1,
    1.40264691163389668864E1,
    9.99999999999999995305E-1,
];

const AD: [f64; 8] = [
    5.67594532638770212846E-1,
    1.47562562584847203173E1,
    8.45138970141474626562E1,
    1.77318088145400459522E2,
    1.64234692871529701831E2,
    7.14778400825575695274E1,
    1.40959135607834029598E1,
    1.00000000000000000470E0,
];

const APN: [f64; 8] = [
    6.13759184814035759225E-1,
    1.47454670787755323881E1,
    8.20584123476060982430E1,
    1.71184781360976385540E2,
    1.59317847137141783523E2,
    6.99778599330103016170E1,
    1.39470856980481566958E1,
    1.00000000000000000550E0,
];

const APD: [f64; 8] = [
    3.34203677749736953049E-1,
    1.11810297306158156705E1,
    7.11727352147859965283E1,
    1.58778084372838313640E2,
    1.53206427475809220834E2,
    6.86752304592780337944E1,
    1.38498634758259442477E1,
    9.99999999999999994502E-1,
];

const BN16: [f64; 5] = [
    -2.53240795869364152689E-1,
    5.75285167332467384228E-1,
    -3.29907036873225371650E-1,
    6.44404068948199951727E-2,
    -3.82519546641336734394E-3,
];

const BD16: [f64; 5] = [
    /* 1.00000000000000000000E0, */
    -7.15685095054035237902E0,
    1.06039580715664694291E1,
    -5.23246636471251500874E0,
    9.57395864378383833152E-1,
    -5.50828147163549611107E-2,
];

const BPPN: [f64; 5] = [
    4.65461162774651610328E-1,
    -1.08992173800493920734E0,
    6.38800117371827987759E-1,
    -1.26844349553102907034E-1,
    7.62487844342109852105E-3,
];

const BPPD: [f64; 5] = [
    /* 1.00000000000000000000E0, */
    -8.70622787633159124240E0,
    1.38993162704553213172E1,
    -7.14116144616431159572E0,
    1.34008595960680518666E0,
    -7.84273211323341930448E-2,
];

const AFN: [f64; 9] = [
    -1.31696323418331795333E-1,
    -6.26456544431912369773E-1,
    -6.93158036036933542233E-1,
    -2.79779981545119124951E-1,
    -4.91900132609500318020E-2,
    -4.06265923594885404393E-3,
    -1.59276496239262096340E-4,
    -2.77649108155232920844E-6,
    -1.67787698489114633780E-8,
];

const AFD: [f64; 9] = [
    /* 1.00000000000000000000E0, */
    1.33560420706553243746E1,
    3.26825032795224613948E1,
    2.67367040941499554804E1,
    9.18707402907259625840E0,
    1.47529146771666414581E0,
    1.15687173795188044134E-1,
    4.40291641615211203805E-3,
    7.54720348287414296618E-5,
    4.51850092970580378464E-7,
];

const AGN: [f64; 11] = [
    1.97339932091685679179E-2,
    3.91103029615688277255E-1,
    1.06579897599595591108E0,
    9.39169229816650230044E-1,
    3.51465656105547619242E-1,
    6.33888919628925490927E-2,
    5.85804113048388458567E-3,
    2.82851600836737019778E-4,
    6.98793669997260967291E-6,
    8.11789239554389293311E-8,
    3.41551784765923618484E-10,
];

const AGD: [f64; 10] = [
    /*  1.00000000000000000000E0, */
    9.30892908077441974853E0,
    1.98352928718312140417E1,
    1.55646628932864612953E1,
    5.47686069422975497931E0,
    9.54293611618961883998E-1,
    8.64580826352392193095E-2,
    4.12656523824222607191E-3,
    1.01259085116509135510E-4,
    1.17166733214413521882E-6,
    4.91834570062930015649E-9,
];

const APFN: [f64; 9] = [
    1.85365624022535566142E-1,
    8.86712188052584095637E-1,
    9.87391981747398547272E-1,
    4.01241082318003734092E-1,
    7.10304926289631174579E-2,
    5.90618657995661810071E-3,
    2.33051409401776799569E-4,
    4.08718778289035454598E-6,
    2.48379932900442457853E-8,
];

const APFD: [f64; 9] = [
    /*  1.00000000000000000000E0, */
    1.47345854687502542552E1,
    3.75423933435489594466E1,
    3.14657751203046424330E1,
    1.09969125207298778536E1,
    1.78885054766999417817E0,
    1.41733275753662636873E-1,
    5.44066067017226003627E-3,
    9.39421290654511171663E-5,
    5.65978713036027009243E-7,
];

const APGN: [f64; 11] = [
    -3.55615429033082288335E-2,
    -6.37311518129435504426E-1,
    -1.70856738884312371053E0,
    -1.50221872117316635393E0,
    -5.63606665822102676611E-1,
    -1.02101031120216891789E-1,
    -9.48396695961445269093E-3,
    -4.60325307486780994357E-4,
    -1.14300836484517375919E-5,
    -1.33415518685547420648E-7,
    -5.63803833958893494476E-10,
];

// NOTE: In C++ this is declared as `double[11]` with only 10 initializers (so the last element is
// an implicit 0.0), but `p1evl(zz, airy_APGD, 10)` only reads the first 10.
const APGD: [f64; 10] = [
    /*  1.00000000000000000000E0, */
    9.85865801696130355144E0,
    2.16401867356585941885E1,
    1.73130776389749389525E1,
    6.17872175280828766327E0,
    1.08848694396321495475E0,
    9.95005543440888479402E-2,
    4.78468199683886610842E-3,
    1.18159633322838625562E-4,
    1.37480673554219441465E-6,
    5.79912514929147598821E-9,
];

/// Airy function
///
/// Solution of the differential equation
///
/// ```text
///     y"(x) = xy.
/// ```
///
/// The function returns the two independent solutions Ai, Bi
/// and their first derivatives Ai'(x), Bi'(x).
///
/// Evaluation is by power series summation for small x,
/// by rational minimax approximations for large x.
///
/// ACCURACY:
///
/// Error criterion is absolute when function <= 1, relative
/// when function > 1, except * denotes relative error criterion.
/// For large negative x, the absolute error increases as x^1.5.
/// For large positive x, the relative error increases as x^1.5.
///
/// ```text
/// Arithmetic  domain   function  # trials      peak         rms
/// IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
/// IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
/// IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
/// IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
/// IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
/// IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16
/// ```
///
/// Returns `(status, ai, aip, bi, bip)`, where `status` is `-1` if `x > 103.892` (in which case
/// `(ai, aip, bi, bip) = (0, 0, inf, inf)`), and `0` otherwise. All four outputs are always
/// written (also for NaN input, which yields NaN outputs through the power series).
#[allow(clippy::too_many_lines)]
pub(crate) fn airy(x: f64) -> (i32, f64, f64, f64, f64) {
    // The C++ writes these through output pointers. Each one is written on every path before
    // returning, so the initial NaN values are never observed.
    let mut ai = f64::NAN;
    let mut aip = f64::NAN;
    let mut bi = f64::NAN;
    let mut bip = f64::NAN;

    let mut domflg: i32 = 0;
    if x > MAXAIRY {
        return (-1, 0.0, 0.0, f64::INFINITY, f64::INFINITY);
    }

    if x < -2.09 {
        // domflg = 15; (unused, since this branch always returns)
        let t = (-x).sqrt();
        let zeta = -2.0 * x * t / 3.0;
        let t = t.sqrt();
        let k = SQRT1OPI / t;
        let z = 1.0 / zeta;
        let zz = z * z;
        let uf = 1.0 + zz * polevl(zz, &AFN) / p1evl(zz, &AFD);
        let ug = z * polevl(zz, &AGN) / p1evl(zz, &AGD);
        let theta = zeta + 0.25 * PI;
        let f = theta.sin();
        let g = theta.cos();
        ai = k * (f * uf - g * ug);
        bi = k * (g * uf + f * ug);
        let uf = 1.0 + zz * polevl(zz, &APFN) / p1evl(zz, &APFD);
        let ug = z * polevl(zz, &APGN) / p1evl(zz, &APGD);
        let k = SQRT1OPI * t;
        aip = -k * (g * uf + f * ug);
        bip = k * (f * uf - g * ug);
        return (0, ai, aip, bi, bip);
    }

    if x >= 2.09 {
        /* cbrt(9) */
        domflg = 5;
        let t = x.sqrt();
        let zeta = 2.0 * x * t / 3.0;
        let g = zeta.exp();
        let t = t.sqrt();
        let k = 2.0 * t * g;
        let z = 1.0 / zeta;
        let f = polevl(z, &AN) / polevl(z, &AD);
        ai = SQRT1OPI * f / k;
        let k = -0.5 * SQRT1OPI * t / g;
        let f = polevl(z, &APN) / polevl(z, &APD);
        aip = f * k;

        if x > 8.3203353 {
            /* zeta > 16 */
            let f = z * polevl(z, &BN16) / p1evl(z, &BD16);
            let k = SQRT1OPI * g;
            bi = k * (1.0 + f) / t;
            let f = z * polevl(z, &BPPN) / p1evl(z, &BPPD);
            bip = k * t * (1.0 + f);
            return (0, ai, aip, bi, bip);
        }
    }

    let mut f = 1.0;
    let mut g = x;
    let mut t = 1.0;
    let mut uf = 1.0;
    let mut ug = x;
    let mut k = 1.0;
    let z = x * x * x;
    while t > MACHEP {
        uf *= z;
        k += 1.0;
        uf /= k;
        ug *= z;
        k += 1.0;
        ug /= k;
        uf /= k;
        f += uf;
        k += 1.0;
        ug /= k;
        g += ug;
        t = (uf / f).abs();
    }
    uf = C1 * f;
    ug = C2 * g;
    if (domflg & 1) == 0 {
        ai = uf - ug;
    }
    if (domflg & 2) == 0 {
        bi = SQRT3 * (uf + ug);
    }

    /* the deriviative of ai */
    k = 4.0;
    uf = x * x / 2.0;
    ug = z / 3.0;
    f = uf;
    g = 1.0 + ug;
    uf /= 3.0;
    t = 1.0;

    while t > MACHEP {
        uf *= z;
        ug /= k;
        k += 1.0;
        ug *= z;
        uf /= k;
        f += uf;
        k += 1.0;
        ug /= k;
        uf /= k;
        g += ug;
        k += 1.0;
        t = (ug / g).abs();
    }

    uf = C1 * f;
    ug = C2 * g;
    if (domflg & 4) == 0 {
        aip = uf - ug;
    }
    if (domflg & 8) == 0 {
        bip = SQRT3 * (uf + ug);
    }
    (0, ai, aip, bi, bip)
}
