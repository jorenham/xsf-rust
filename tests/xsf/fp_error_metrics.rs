#![allow(clippy::float_cmp)]

use num_complex::c64;

#[test]
fn test_extended_absolute_error_f64() {
    assert_eq!(xsf::extended_absolute_error(0.0, 0.0), 0.0);
    assert_eq!(xsf::extended_absolute_error(1.0, 0.0), 1.0);
    assert_eq!(xsf::extended_absolute_error(1.0, 2.0), 1.0);
    assert_eq!(xsf::extended_absolute_error(2.0, 1.0), 1.0);
    assert_eq!(xsf::extended_absolute_error(3.0, 1.0), 2.0);
}

#[test]
fn test_extended_absolute_error_c64() {
    assert_eq!(
        xsf::extended_absolute_error(c64(1.0, 1.0), c64(1.0, 1.0)),
        0.0
    );
    assert_eq!(
        xsf::extended_absolute_error(c64(0.0, 0.0), c64(3.0, 4.0)),
        5.0
    );
}

#[test]
fn test_extended_relative_error_f64() {
    assert_eq!(xsf::extended_relative_error(0.0, 0.0), 0.0);
    assert_eq!(xsf::extended_relative_error(1.0, 0.0), f64::INFINITY);
    assert_eq!(xsf::extended_relative_error(1.0, 2.0), 0.5);
    assert_eq!(xsf::extended_relative_error(2.0, 1.0), 1.0);
    assert_eq!(xsf::extended_relative_error(3.0, 1.0), 2.0);
}

#[test]
fn test_extended_relative_error_c64() {
    assert_eq!(
        xsf::extended_relative_error(c64(1.0, 1.0), c64(1.0, 1.0)),
        0.0
    );
    assert_eq!(
        xsf::extended_relative_error(c64(0.0, 0.0), c64(3.0, 4.0)),
        1.0
    );
}
