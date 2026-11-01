use num_complex::c64;

#[test]
fn test_cevalpoly_0() {
    // p(z) = 0
    let y = xsf::cevalpoly(&[], c64(2.0, 3.0));
    assert_eq!(y, c64(0.0, 0.0));
}

#[test]
fn test_cevalpoly_1() {
    // p(z) = 5
    let y = xsf::cevalpoly(&[5.0], c64(2.0, 3.0));
    // p(2+3i) = 5
    assert_eq!(y, c64(5.0, 0.0));
}

#[test]
fn test_cevalpoly_2() {
    // p(z) = 2z + 3
    let y = xsf::cevalpoly(&[2.0, 3.0], c64(1.0, 1.0));
    // p(1+i) = 5 + 2i
    assert_eq!(y, c64(5.0, 2.0));
}

#[test]
fn test_cevalpoly_4() {
    // p(z) = z^3 + 2z^2 + 3z + 4
    let y = xsf::cevalpoly(&[1.0, 2.0, 3.0, 4.0], c64(1.0, 1.0));
    // p(1+i) = 5 + 9i
    assert_eq!(y, c64(5.0, 9.0));
}
