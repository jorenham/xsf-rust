use num_complex::c64;
use xsref::np_assert_allclose;

#[test]
fn test_spence_xsref_f64() {
    xsref::test("spence", "d-d", |x| xsf::spence(x[0]));
}

#[test]
fn test_spence_xsref_c64() {
    xsref::test("spence", "d-d", |x| xsf::spence(c64(x[0], 0.0)).re);
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
