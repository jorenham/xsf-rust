#[test]
fn test_lgam1p() {
    xsref::test("lgam1p", "d-d", |x| xsf::lgam1p(x[0]));
}

#[test]
fn test_lgam1p_half() {
    // mpmath: loggamma(1.5), loggamma(0.5)
    for (x, expected) in [
        (0.5, -0.120_782_237_635_245_22),
        (-0.5, 0.572_364_942_924_700_1),
    ] {
        let rel = (xsf::lgam1p(x) / expected - 1.0).abs();
        assert!(rel < 1e-15, "x = {x}, rel = {rel:e}");
    }
}
