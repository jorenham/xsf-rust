use num_complex::c64;

#[test]
fn test_hyp2f1_c64_transform2_euler() {
    // unlike in xsf: the 1/z transformation of the Euler transform (mpmath references)
    for (a, b, c, z, expected) in [
        // a - b an integer: the limiting case
        (
            8.5,
            6.5,
            -4.25,
            c64(-1.0, 0.5),
            c64(-14.914_154_962_801_387, 2.508_617_290_210_916_7),
        ),
        (
            4.3,
            4.4,
            -4.2,
            c64(-0.9, 0.8),
            c64(-13.468_832_852_470_712, -27.983_525_735_799_486),
        ),
        // the limiting case with b within 1e-15 of an integer
        (
            8.881_784_197_001_252e-16,
            8.881_784_197_001_252e-16,
            -1.25,
            c64(-1.0, 0.5),
            c64(1.0, -4.661_160_101_737_811e-31),
        ),
        // (1 - z)^(c - a - b) is subnormal
        (
            9.48,
            9.48,
            -7.69,
            c64(2.97e11, 8.92e11),
            c64(3.884_057_170_000_241_6e-108, -6.451_235_141_304_403e-108),
        ),
        // c - (c - a) loses a
        (
            5e-6,
            -3.125,
            -9.0625,
            c64(133.5, -390.75),
            c64(3.466_873_805_327_974, -0.823_490_927_100_644_6),
        ),
    ] {
        let actual = xsf::hyp2f1(a, b, c, z);
        assert!(
            (actual - expected).norm() < 1e-13 * expected.norm(),
            "hyp2f1({a}, {b}, {c}, {z}) = {actual} != {expected}"
        );
    }
}
