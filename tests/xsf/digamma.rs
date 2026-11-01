#[test]
fn test_digamma_f64() {
    xsref::test("digamma", "d-d", |x| xsf::digamma(x[0]));
}

#[test]
fn test_digamma_near_negative_pole() {
    // just above a negative integer (reference values from mpmath)
    for (x, expected) in [
        (-0.999_999, -999_999.577_184_264_3),
        (-2.999_999_9, -9_999_998.760_247_82),
    ] {
        let rel = ((xsf::digamma(x) - expected) / expected).abs();
        assert!(rel < 1e-15, "x = {x}: rel = {rel:e}");
    }
}

#[test]
fn test_digamma_c64_abs16() {
    // the forward recurrence boundary at |z| = 16 (reference values from mpmath)
    let w = xsf::digamma(num_complex::c64(-1e-300, 16.0));
    let expected = num_complex::c64(2.772_914_370_467_193_5, 1.602_046_326_794_896_6);
    let rel = (w - expected).norm() / expected.norm();
    assert!(rel < 1e-15, "rel = {rel:e}");
}

#[test]
fn test_digamma_c64() {
    use num_complex::c64;

    xsref::test("digamma", "cd-cd", |x| xsf::digamma(c64(x[0], x[1])));
}
