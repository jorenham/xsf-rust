use core::f64::consts::PI;

use num_complex::c64;

#[test]
fn test_trigamma_f64() {
    // reference values from xsf's test_trigamma.cpp
    for (x, expected) in [
        (0.1, 101.433_299_150_792_76),
        (0.5, PI * PI / 2.0),
        (2.0, PI * PI / 6.0 - 1.0),
        (-7.5, 9.744_766_282_170_433),
        (-3.2, 28.298_186_402_194_08),
    ] {
        let rel = xsf::extended_relative_error(xsf::trigamma(x), expected);
        assert!(rel <= 1e-14, "x = {x}: rel = {rel:e}");
    }
}

#[test]
fn test_trigamma_f64_poles() {
    assert_eq!(xsf::trigamma(0.0), f64::INFINITY);
    assert_eq!(xsf::trigamma(-1.0), f64::INFINITY);
}

#[test]
fn test_trigamma_c64() {
    // reference values from xsf's test_trigamma.cpp
    for (z, expected) in [
        (c64(0.0, 8.0), c64(-0.007_812_5, -0.124_673_450_303_127_63)),
        (
            c64(-3.2, 0.1),
            c64(15.207_350_644_973_363, 15.708_103_885_511_357),
        ),
    ] {
        let rel = xsf::extended_relative_error(xsf::trigamma(z), expected);
        assert!(rel <= 1e-14, "z = {z}: rel = {rel:e}");
    }
}
