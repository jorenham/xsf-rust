//! Translated into pure Rust from `xsf/multivariate_normal.h` (xsf main).
//!
//! Unlike in xsf (scipy/xsf#360), `bivariate_normal_cdf` returns NaN if any argument is NaN, and
//! treats `|h| > 40` and `|k| > 40` as infinite.
//!
//! `bivariate_normal_cdf` is adapted from the original MATLAB survival function implementation
//! by Dr. Alan Genz, and `trivariate_normal_cdf` is a port of `tvnl` in Alan Genz's MATLAB
//! `tvn.m`. Original source:
//! <https://web.archive.org/web/20200205123040/http://www.math.wsu.edu/faculty/genz/software/matlab/tvn.m>
//!
//! > Copyright (C) 2011, Alan Genz, All rights reserved.
//! >
//! > Redistribution and use in source and binary forms, with or without
//! > modification, are permitted provided the following conditions are met:
//! >   1. Redistributions of source code must retain the above copyright
//! >      notice, this list of conditions and the following disclaimer.
//! >   2. Redistributions in binary form must reproduce the above copyright
//! >      notice, this list of conditions and the following disclaimer in the
//! >      documentation and/or other materials provided with the distribution.
//! >   3. The contributor name(s) may not be used to endorse or promote
//! >      products derived from this software without specific prior written
//! >      permission.
//! >
//! > THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//! > "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//! > LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//! > FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//! > COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//! > INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//! > BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
//! > OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//! > ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
//! > TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
//! > USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use core::f64::consts::{FRAC_PI_2, PI};

use crate::xsf::cephes::ndtr;

/// Handles degenerate cases for `bivariate_normal_sf` (infinite arguments or zero correlation).
/// Returns `Some(p)` if a boundary case applies, `None` otherwise.
fn bivariate_normal_sf_boundary(dh: f64, dk: f64, r: f64) -> Option<f64> {
    if dh == f64::INFINITY || dk == f64::INFINITY {
        // if dh ==  inf | dk ==  inf:p = 0;
        return Some(0.0);
    }
    if dh == f64::NEG_INFINITY {
        // elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
        return Some(if dk == f64::NEG_INFINITY {
            1.0
        } else {
            ndtr(-dk)
        });
    }
    if dk == f64::NEG_INFINITY {
        // elseif dk == -inf, p = phid(-dh);
        return Some(ndtr(-dh));
    }
    if r == 0.0 {
        // elseif r == 0, p = phid(-dh)*phid(-dk);
        return Some(ndtr(-dh) * ndtr(-dk));
    }
    None
}

// Positive nodes and weights of symmetric Gauss-Legendre rules on [-1, 1].
const BVN_W6: [f64; 3] = [
    0.171_324_492_379_170_5,
    0.360_761_573_048_138_4,
    0.467_913_934_572_690_4,
];
const BVN_X6: [f64; 3] = [
    0.932_469_514_203_152_2,
    0.661_209_386_466_264_7,
    0.238_619_186_083_197,
];

const BVN_W12: [f64; 6] = [
    0.047_175_336_386_511_77,
    0.106_939_325_995_318_3,
    0.160_078_328_543_346_4,
    0.203_167_426_723_065_9,
    0.233_492_536_538_354_7,
    0.249_147_045_813_402_9,
];
const BVN_X12: [f64; 6] = [
    0.981_560_634_246_719_1,
    0.904_117_256_370_475,
    0.769_902_674_194_305,
    0.587_317_954_286_617_1,
    0.367_831_498_998_180_2,
    0.125_233_408_511_469_2,
];

const BVN_W20: [f64; 10] = [
    0.017_614_007_139_152_12,
    0.040_601_429_800_386_94,
    0.062_672_048_334_109_05,
    0.083_276_741_576_704_75,
    0.101_930_119_817_240_4,
    0.118_194_531_961_518_4,
    0.131_688_638_449_176_6,
    0.142_096_109_318_382_1,
    0.149_172_986_472_603_7,
    0.152_753_387_130_725_9,
];
const BVN_X20: [f64; 10] = [
    0.993_128_599_185_094_9,
    0.963_971_927_277_913_8,
    0.912_234_428_251_326,
    0.839_116_971_822_218_8,
    0.746_331_906_460_150_8,
    0.636_053_680_726_515,
    0.510_867_001_950_827_1,
    0.373_706_088_715_419_6,
    0.227_785_851_141_645_1,
    0.076_526_521_133_497_32,
];

// Nodes and weights of the 20-point Gauss-Laguerre rule on [0, inf), with weight e^-t.
const BVN_LAGUERRE_X20: [f64; 20] = [
    0.070_539_889_691_988_75,
    0.372_126_818_001_611_46,
    0.916_582_102_483_273_5,
    1.707_306_531_028_344,
    2.749_199_255_309_432,
    4.048_925_313_850_887,
    5.615_174_970_861_616_5,
    7.459_017_453_671_064,
    9.594_392_869_581_096,
    12.038_802_546_964_316,
    14.814_293_442_630_74,
    17.948_895_520_519_375,
    21.478_788_240_285_01,
    25.451_702_793_186_904,
    29.932_554_631_700_61,
    35.013_434_240_479,
    40.833_057_056_728_57,
    47.619_994_047_346_5,
    55.810_795_750_063_896,
    66.524_416_525_615_75,
];
const BVN_LAGUERRE_W20: [f64; 20] = [
    0.168_746_801_851_113_88,
    0.291_254_362_006_068_3,
    0.266_686_102_867_001_3,
    0.166_002_453_269_506_83,
    0.074_826_064_668_792_37,
    0.024_964_417_309_283_22,
    0.006_202_550_844_572_237,
    0.001_144_962_386_476_908_2,
    0.000_155_741_773_027_811_97,
    1.540_144_086_522_491_5e-5,
    1.086_486_366_517_982_4e-6,
    5.330_120_909_556_714_6e-8,
    1.757_981_179_050_582e-9,
    3.725_502_402_512_320_6e-11,
    4.767_529_251_578_191e-13,
    3.372_844_243_362_438_6e-15,
    1.155_014_339_500_398_9e-17,
    1.539_522_140_582_343_5e-20,
    5.286_442_725_569_158e-24,
    1.656_456_612_499_023_3e-28,
];

/// `P(X > h, Y > k)` by Gauss-Laguerre quadrature, if `-1 < r < 1`, `h > r k` and `k > r h`
/// (i.e. the density on the region peaks at its corner `(h, k)`), and `d > 5` (see below).
/// Returns `None` otherwise.
///
/// - Since `h + k > 0`, `P = 0` at `r = -1`, so Plackett's identity gives
///   `P = int_-1^r phi2(h, k; s) ds`.
/// - The exponent of `phi2` is `a / (1 + s) + b / (1 - s)`, with `a = (h + k)^2 / 4` and
///   `b = (h - k)^2 / 4`.
/// - Substituting `1 / (1 + s) = 1 / (1 + r) + t / a` gives `int_0^inf e^-t f(t) dt`, with no
///   cancellation, and `f` smooth if `d = a (1 - r) / (2 (1 + r))` is large: its nearest
///   singularity is at `t = -d`. For `d > 5`, the quadrature error is ~1e-14 (plus the
///   rounding error of the exponent, ~1e-13 for tiny `P`).
fn bivariate_normal_sf_tail(h: f64, k: f64, r: f64) -> Option<f64> {
    if !(r > -1.0 && h > r * k && k > r * h && (h + k) * (h + k) * (1.0 - r) > 40.0 * (1.0 + r)) {
        return None;
    }
    let a = (h + k) * (h + k) / 4.0;
    let b = (h - k) * (h - k) / 4.0;
    let w0 = 1.0 / (1.0 + r);
    let mut sum = 0.0;
    for (&t, &wt) in BVN_LAGUERRE_X20.iter().zip(&BVN_LAGUERRE_W20) {
        let w = w0 + t / a;
        sum += wt * (-a * w0 - b * w / (2.0 * w - 1.0)).exp() / (a * w * (2.0 * w - 1.0).sqrt());
    }
    Some(sum / (2.0 * PI))
}

/// Cumulative distribution function of the standard bivariate normal distribution
///
/// Returns $P(X \le h, Y \le k)$, with $h$ = `dh` and $k$ = `dk`, for a standard bivariate
/// normal vector $(X, Y)$ with correlation $-1 \le r \le 1$.
///
/// Corresponds to `xsf::bivariate_normal_cdf` in xsf main, which is used by
/// [`scipy.stats.multivariate_normal`][multivariate_normal].
///
/// [multivariate_normal]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multivariate_normal.html
#[must_use]
#[allow(clippy::many_single_char_names)]
pub fn bivariate_normal_cdf(dh: f64, dk: f64, r: f64) -> f64 {
    // Cumulative distribution function of the standard bivariate normal distribution.
    //
    // Return P(X <= dh, Y <= dk) for a standard bivariate normal vector
    // (X, Y) with correlation r.
    //
    // dh, dk are the upper limits of the lower tail, and r must satisfy
    // -1 <= r <= 1.
    //
    // Adapted from the original MATLAB survival function implementation by Dr. Alan Genz;
    // see license information in LICENSES_bundled.txt.
    // In the comments, phid is the CDF of the standard normal distribution.

    // Unlike in xsf (scipy/xsf#360), NaN input returns NaN, e.g. `(1, -1, NaN)` was 0.68
    if dh.is_nan() || dk.is_nan() || r.is_nan() {
        return f64::NAN;
    }
    // Unlike in xsf (scipy/xsf#360), `|dh| > 40` is treated as infinite, so that `hk` etc. cannot
    // overflow, e.g. `(1e70, 38, 0.97)` was NaN instead of 1.
    let dh = if dh.abs() > 40.0 {
        dh * f64::INFINITY
    } else {
        dh
    };
    let dk = if dk.abs() > 40.0 {
        dk * f64::INFINITY
    } else {
        dk
    };

    // By symmetry, the CDF at (dh, dk) is the upper-tail probability at (-dh, -dk).
    let dh = -dh;
    let dk = -dk;
    if let Some(p) = bivariate_normal_sf_boundary(dh, dk, r) {
        return p;
    }
    // Unlike in xsf (scipy/xsf#360), the tail where the density peaks at the corner (h, k) is
    // computed without cancellation, see `bivariate_normal_sf_tail`
    if let Some(p) = bivariate_normal_sf_tail(dh, dk, r) {
        return p;
    }
    // else, tp = 2*pi; h = dh; k = dk; hk = h*k; bvn = 0;
    let tp = 2.0 * PI;
    let h = dh;
    let mut k = dk;
    let mut hk = h * k;
    let mut bvn = 0.0;
    // Unlike in xsf (scipy/xsf#360), the 20-point rule is also used in the tail
    let tail = r.abs() * (h * h + k * k) > 2.0;
    let (w, x): (&[f64], &[f64]) = if r.abs() < 0.3 && !tail {
        // Gauss Legendre points and weights, n = 6
        (&BVN_W6, &BVN_X6)
    } else if r.abs() < 0.75 && !tail {
        // Gauss Legendre points and weights, n = 12
        (&BVN_W12, &BVN_X12)
    } else {
        // Gauss Legendre points and weights, n = 20
        (&BVN_W20, &BVN_X20)
    };

    // if abs(r) < 0.925, hs = ( h*h + k*k )/2; asr = asin(r)/2;
    if r.abs() < 0.925 {
        let hs = f64::midpoint(h * h, k * k);
        let asr = r.asin() / 2.0;
        // sn = sin(asr*x); bvn = exp((sn*hk-hs)./(1-sn.^2))*w';
        for (&w_i, &x_i) in w.iter().zip(x) {
            for x_ij in [1.0 - x_i, 1.0 + x_i] {
                let sn = (asr * x_ij).sin();
                bvn += ((sn * hk - hs) / (1.0 - sn * sn)).exp() * w_i;
            }
        }
        // bvn = bvn*asr/tp + phid(-h)*phid(-k);
        bvn = bvn * asr / tp + ndtr(-h) * ndtr(-k);
    } else {
        // else, if r < 0, k = -k; hk = -hk; end
        if r < 0.0 {
            k = -k;
            hk = -hk;
        }
        if r.abs() < 1.0 {
            // if abs(r) < 1, as = 1-r^2; a = sqrt(as); bs = (h-k)^2;
            let as_ = 1.0 - r * r;
            let mut a = as_.sqrt();
            let bs = (h - k) * (h - k);
            // asr = -( bs/as + hk )/2; c = (4-hk)/8 ; d = (12-hk)/80;
            let asr = -(bs / as_ + hk) / 2.0;
            let c = (4.0 - hk) / 8.0;
            let d = (12.0 - hk) / 80.0;
            if asr > -100.0 {
                // if asr > -100, bvn = a*exp(asr)*(1-c*(bs-as)*(1-d*bs)/3+c*d*as^2); end
                bvn = a
                    * asr.exp()
                    * (1.0 - c * (bs - as_) * (1.0 - d * bs) / 3.0 + c * d * as_ * as_);
            }
            if hk > -100.0 {
                // if hk  > -100, b = sqrt(bs); sp = sqrt(tp)*phid(-b/a);
                let b = bs.sqrt();
                let sp = tp.sqrt() * ndtr(-b / a);
                // bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1-d*bs)/3 );
                bvn -= (-hk / 2.0).exp() * sp * b * (1.0 - c * bs * (1.0 - d * bs) / 3.0);
            }
            // end, a = a/2; xs = (a*x).^2; asr = -( bs./xs + hk )/2;
            a /= 2.0;

            // ix = find( asr > -100 ); xs = xs(ix); sp = ( 1 + c*xs.*(1+5*d*xs) );
            let mut tmp = 0.0;
            for (&w_i, &x_i) in w.iter().zip(x) {
                for x_ij in [1.0 - x_i, 1.0 + x_i] {
                    let xs_i = (a * x_ij) * (a * x_ij);
                    let asr_i = -(bs / xs_i + hk) / 2.0;
                    if asr_i <= -100.0 {
                        continue;
                    }
                    let sp = 1.0 + c * xs_i * (1.0 + 5.0 * d * xs_i);
                    // rs = sqrt(1-xs); ep = exp( -(hk/2)*xs./(1+rs).^2 )./rs;
                    let rs = (1.0 - xs_i).sqrt();
                    let ep = (-(hk / 2.0) * xs_i / ((1.0 + rs) * (1.0 + rs))).exp() / rs;
                    tmp += w_i * asr_i.exp() * (sp - ep);
                }
            }

            // bvn = ( a*( (exp(asr(ix)).*(sp-ep))*w(ix)' ) - bvn )/tp;
            bvn = (a * tmp - bvn) / tp;
        }
        // end
        // if r > 0, bvn =  bvn + phid( -max( h, k ) );
        if r > 0.0 {
            bvn += ndtr(-(if h > k { h } else { k }));
        } else if h >= k {
            // elseif h >= k, bvn = -bvn;
            bvn = -bvn;
        } else {
            // else, if h < 0, L = phid(k)-phid(h); else, L = phid(-h)-phid(-k); end
            let l = if h < 0.0 {
                ndtr(k) - ndtr(h)
            } else {
                ndtr(-h) - ndtr(-k)
            };
            // bvn =  L - bvn;
            bvn = l - bvn;
        }
    }
    // end, p = max( 0, min( 1, bvn ) );
    bvn.clamp(0.0, 1.0)
}

/// Computes `sin(x)` and `cos(x)^2` with series approximation for `|x|` near `pi/2`
fn sin_cos2(x: f64) -> (f64, f64) {
    let ee = (FRAC_PI_2 - x.abs()) * (FRAC_PI_2 - x.abs());
    if ee < 5e-5 {
        let cs = ee * (1.0 - ee * (1.0 - 2.0 * ee / 15.0) / 3.0);
        let sgn = f64::from(u8::from(x > 0.0)) - f64::from(u8::from(x < 0.0));
        let sx = (1.0 - ee * (1.0 - ee / 12.0) / 2.0) * sgn;
        (sx, cs)
    } else {
        let sx = x.sin();
        (sx, 1.0 - sx * sx)
    }
}

/// Evaluate one Plackett formula integrand term.
fn plackett_term(ba: f64, bb: f64, bc: f64, ra: f64, rb: f64, r: f64, rr: f64) -> f64 {
    let mut f = 0.0;
    let dt = rr * (rr - (ra - rb) * (ra - rb) - 2.0 * ra * rb * (1.0 - r));
    if dt > 0.0 {
        let bt = (bc * rr + ba * (r * rb - ra) + bb * (r * ra - rb)) / dt.sqrt();
        let ft = (ba - r * bb) * (ba - r * bb) / rr + bb * bb;
        if bt > -10.0 && ft < 100.0 {
            f = (-ft / 2.0).exp();
            if bt < 10.0 {
                f *= ndtr(bt);
            }
        }
    }
    f
}

/// Combine the two Plackett integrand terms at integration point `x`.
fn tvn_integrand(x: f64, h1: f64, h2: f64, h3: f64, r23: f64, a12: f64, a13: f64) -> f64 {
    let mut f = 0.0;
    let (r12, rr2) = sin_cos2(a12 * x);
    let (r13, rr3) = sin_cos2(a13 * x);
    if a12.abs() > 0.0 {
        f += a12 * plackett_term(h1, h2, h3, r13, r23, r12, rr2);
    }
    if a13.abs() > 0.0 {
        f += a13 * plackett_term(h1, h3, h2, r12, r23, r13, rr3);
    }
    f
}

/// Apply a Gauss-Kronrod rule on one adaptive subinterval, returning `(resk, err)`.
#[allow(clippy::too_many_arguments)]
fn tvn_gauss_kronrod(
    a: f64,
    b: f64,
    h1: f64,
    h2: f64,
    h3: f64,
    r23: f64,
    a12: f64,
    a13: f64,
) -> (f64, f64) {
    const WG0: f64 = 0.272_925_086_777_900_7;
    const WG: [f64; 5] = [
        0.055_668_567_116_174_49,
        0.125_580_369_464_904_8,
        0.186_290_210_927_735_2,
        0.233_193_764_591_991_4,
        0.262_804_544_510_247_8,
    ];
    const XGK: [f64; 11] = [
        0.996_369_613_889_542_7,
        0.978_228_658_146_057,
        0.941_677_108_578_068_1,
        0.887_062_599_768_095_3,
        0.816_057_456_656_221_1,
        0.730_152_005_574_049_2,
        0.630_599_520_161_965_1,
        0.519_096_129_206_811_8,
        0.397_944_140_952_377_6,
        0.269_543_155_952_345,
        0.136_113_000_799_361_7,
    ];
    const WGK0: f64 = 0.136_577_794_711_118_3;
    const WGK: [f64; 11] = [
        0.009_765_441_045_961_29,
        0.027_156_554_682_104_43,
        0.045_829_378_564_426_71,
        0.063_097_424_750_374_84,
        0.078_664_571_932_227_64,
        0.092_953_098_596_900_74,
        0.105_872_074_481_389_4,
        0.116_739_502_461_047_2,
        0.125_158_799_100_319_5,
        0.131_280_684_229_805_7,
        0.135_193_572_799_884_5,
    ];

    let integrand = |x| tvn_integrand(x, h1, h2, h3, r23, a12, a13);

    let wid = (b - a) / 2.0;
    let cen = f64::midpoint(b, a);
    let mut fc = integrand(cen);
    let mut resg = fc * WG0;
    let mut resk = fc * WGK0;

    for j in 0..5 {
        let mut t = wid * XGK[2 * j];
        fc = integrand(cen - t) + integrand(cen + t);
        resk += WGK[2 * j] * fc;

        t = wid * XGK[2 * j + 1];
        fc = integrand(cen - t) + integrand(cen + t);
        resk += WGK[2 * j + 1] * fc;
        resg += WG[j] * fc;
    }

    let t = wid * XGK[10];
    fc = integrand(cen - t) + integrand(cen + t);
    resk = wid * (resk + WGK[10] * fc);
    let err = (resk - wid * resg).abs();
    (resk, err)
}

/// Adaptively integrate the Plackett integrand over `[0, 1]`.
fn tvn_adaptive_integral(h1: f64, h2: f64, h3: f64, r23: f64, a12: f64, a13: f64, tol: f64) -> f64 {
    const NL: usize = 100;
    let mut ai = [0.0; NL];
    let mut bi = [0.0; NL];
    let mut fi = [0.0; NL];
    let mut ei = [0.0; NL];
    ai[0] = 0.0;
    bi[0] = 1.0;
    fi[0] = 0.0;
    ei[0] = 0.0;

    let mut ip = 0;
    let mut im = 0;
    let mut err = 1.0;
    let mut fin = 0.0;
    while 4.0 * err > tol && im < NL - 1 {
        im += 1;
        bi[im] = bi[ip];
        ai[im] = f64::midpoint(ai[ip], bi[ip]);
        bi[ip] = ai[im];
        (fi[ip], ei[ip]) = tvn_gauss_kronrod(ai[ip], bi[ip], h1, h2, h3, r23, a12, a13);
        (fi[im], ei[im]) = tvn_gauss_kronrod(ai[im], bi[im], h1, h2, h3, r23, a12, a13);

        fin = 0.0;
        let mut err2 = 0.0;
        ip = 0;
        for i in 0..=im {
            fin += fi[i];
            err2 += ei[i] * ei[i];
            if ei[i] > ei[ip] {
                ip = i;
            }
        }
        err = err2.sqrt();
    }
    fin
}

/// Cumulative distribution function of the standard trivariate normal distribution
///
/// Returns $P(X \le h_1, Y \le h_2, Z \le h_3)$ for a standard trivariate normal vector
/// $(X, Y, Z)$ with correlation matrix
///
/// $$
/// \begin{bmatrix}
///     1      & r_{12} & r_{13} \\\\
///     r_{12} & 1      & r_{23} \\\\
///     r_{13} & r_{23} & 1
/// \end{bmatrix} ,
/// $$
///
/// which must be positive semidefinite. `epsi` is the requested absolute accuracy (`1e-7` by
/// default in xsf), with a lower bound of `1e-14`.
///
/// Corresponds to `xsf::trivariate_normal_cdf` in xsf main, which is used by
/// [`scipy.stats.multivariate_normal`][multivariate_normal].
///
/// [multivariate_normal]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multivariate_normal.html
#[must_use]
#[allow(clippy::float_cmp, clippy::too_many_arguments, clippy::too_many_lines)]
pub fn trivariate_normal_cdf(
    h1: f64,
    h2: f64,
    h3: f64,
    r12: f64,
    r13: f64,
    r23: f64,
    epsi: f64,
) -> f64 {
    // Cumulative distribution function of the standard trivariate normal distribution.
    //
    // Return P(X <= h1, Y <= h2, Z <= h3) for a standard trivariate normal
    // vector (X, Y, Z) with correlation matrix
    //
    //     [ 1   r12  r13 ]
    //     [ r12  1   r23 ]
    //     [ r13 r23   1  ].
    //
    // h1, h2, h3 are the upper integration limits. The correlations must
    // describe a positive semidefinite correlation matrix. epsi is the requested
    // absolute accuracy (default 1e-7), with a lower bound of 1e-14.
    //
    // Port of tvnl in Alan Genz's MATLAB tvn.m. See LICENSES_bundled.txt for the license.
    // Original source:
    // https://web.archive.org/web/20200205123040/http://www.math.wsu.edu/faculty/genz/software/matlab/tvn.m
    let (mut h1, mut h2, mut h3) = (h1, h2, h3);
    let (mut r12, mut r13, mut r23) = (r12, r13, r23);
    let epst = if epsi > 1e-14 { epsi } else { 1e-14 };

    if h1.is_nan()
        || h2.is_nan()
        || h3.is_nan()
        || r12.is_nan()
        || r13.is_nan()
        || r23.is_nan()
        || epsi.is_nan()
    {
        return f64::NAN;
    }
    let det = 1.0 + 2.0 * r12 * r13 * r23 - r12 * r12 - r13 * r13 - r23 * r23;
    // Tolerate tiny negative determinants caused by floating-point roundoff
    // when the correlation matrix is singular.
    let det_tol = 8.0 * f64::EPSILON;
    if r12.abs() > 1.0 || r13.abs() > 1.0 || r23.abs() > 1.0 || det < -det_tol {
        return f64::NAN;
    }
    if h1 == f64::NEG_INFINITY || h2 == f64::NEG_INFINITY || h3 == f64::NEG_INFINITY {
        return 0.0;
    }
    if h1 == f64::INFINITY {
        if h2 == f64::INFINITY {
            return if h3 == f64::INFINITY { 1.0 } else { ndtr(h3) };
        }
        if h3 == f64::INFINITY {
            return ndtr(h2);
        }
        return bivariate_normal_cdf(h2, h3, r23);
    }
    if h2 == f64::INFINITY {
        if h3 == f64::INFINITY {
            return ndtr(h1);
        }
        return bivariate_normal_cdf(h1, h3, r13);
    }
    if h3 == f64::INFINITY {
        return bivariate_normal_cdf(h1, h2, r12);
    }

    // Sort correlations and check special cases.
    if r12.abs() > r13.abs() {
        core::mem::swap(&mut h2, &mut h3);
        core::mem::swap(&mut r12, &mut r13);
    }
    if r13.abs() > r23.abs() {
        core::mem::swap(&mut h1, &mut h2);
        core::mem::swap(&mut r23, &mut r13);
    }

    let tvn = if h1.abs() + h2.abs() + h3.abs() < epst {
        (1.0 + 2.0 * (r12.asin() + r13.asin() + r23.asin()) / PI) / 8.0
    } else if r12.abs() + r13.abs() < epst {
        ndtr(h1) * bivariate_normal_cdf(h2, h3, r23)
    } else if r13.abs() + r23.abs() < epst {
        ndtr(h3) * bivariate_normal_cdf(h1, h2, r12)
    } else if r12.abs() + r23.abs() < epst {
        ndtr(h2) * bivariate_normal_cdf(h1, h3, r13)
    } else if 1.0 - r23 < epst {
        bivariate_normal_cdf(h1, if h2 < h3 { h2 } else { h3 }, r12)
    } else if r23 == -1.0 && h2 <= -h3 {
        // When Z = -Y, the constraints require -h3 <= Y <= h2.
        // For h2 <= -h3, this event has probability zero; avoid quadrature.
        0.0
    } else if r23 + 1.0 < epst && h2 > -h3 {
        bivariate_normal_cdf(h1, h2, r12) - bivariate_normal_cdf(h1, -h3, r12)
    } else {
        // Use numerical integration to compute probability and add
        // singular values from the Plackett formula.
        let a12 = r12.asin();
        let a13 = r13.asin();
        let tvn = tvn_adaptive_integral(h1, h2, h3, r23, a12, a13, epst) / (2.0 * PI);
        tvn + bivariate_normal_cdf(h2, h3, r23) * ndtr(h1)
    };

    tvn.clamp(0.0, 1.0)
}
