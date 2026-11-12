use num_complex::c64;
use xsref::np_assert_allclose;

xsref::tests! {
    test_spence_xsref_f64: "spence", "d-d", |x| xsf::spence(x[0]);
    // the complex result isn't real for x < 0, where the real one is NaN
    test_spence_xsref_c64: "spence", "d-d", |x| {
        if x[0] < 0.0 { f64::NAN } else { xsf::spence(c64(x[0], 0.0)).re }
    };
}

#[test]
fn test_spence_manual_c64() {
    // values calculated using Wolfram|Alpha with PolyLog[2, 1-z]
    let zs = [
        c64(-1.0, -1.0),
        c64(-1.0, 1.0),
        c64(0.0, -1.0),
        c64(0.0, 1.0),
        c64(1.0, -1.0),
        c64(1.0, 1.0),
    ];
    let expect = [
        c64(1.186_688_537_000_058, 2.407_740_769_345_772),
        c64(1.186_688_537_000_058, -2.407_740_769_345_772),
        c64(0.616_850_275_068_085, 1.460_362_116_753_12),
        c64(0.616_850_275_068_085, -1.460_362_116_753_12),
        c64(-0.205_616_758_356_028_3, 0.915_965_594_177_219),
        c64(-0.205_616_758_356_028_3, -0.915_965_594_177_219),
    ];
    np_assert_allclose!(zs.map(xsf::spence), expect, rtol = 1e-13);
}

#[test]
fn test_spence_near_1_c64() {
    // mpmath: polylog(2, 1 - z)
    let z = c64(1.0 + 1e-9, 1e-9);
    let expect = c64(-1.000_000_082_740_371e-9, -9.999_999_995e-10);
    np_assert_allclose!([xsf::spence(z)], [expect], rtol = 1e-15);
}

#[test]
fn test_spence_huge_c64() {
    // mpmath: polylog(2, 1 - z), where |z|^2 or |z| overflows
    let zs = [
        c64(1e200, 1e200),
        c64(-3e160, 2e155),
        c64(0.0, -1e300),
        c64(f64::MAX, f64::MAX),
    ];
    let expect = [
        c64(-106_198.961_811_644_83, -361.961_418_882_061_2),
        c64(-68_266.353_166_041_1, -1_160.854_434_866_919_6),
        c64(-238_585.826_205_044_6, 1_085.067_661_862_319_7),
        c64(-252_143.138_269_828_16, -557.734_237_379_010_1),
    ];
    np_assert_allclose!(zs.map(xsf::spence), expect, rtol = 1e-15);
}
