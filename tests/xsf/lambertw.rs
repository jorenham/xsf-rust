use num_complex::c64;
use num_traits::ToPrimitive;

#[test]
fn test_lambertw_c64() {
    xsref::test("lambertw", "cd_p_d-cd", |x| {
        xsf::lambertw(c64(x[0], x[1]), x[2].to_isize().unwrap(), x[3])
    });
}

#[test]
fn test_lambertw_real() {
    // mpmath: lambertw(z, k)
    for (z, k, expected) in [
        (10.0, 0, 1.745_528_002_740_699_4),
        (1e100, 0, 224.843_106_445_118_5),
        (-0.2, 0, -0.259_171_101_819_073_77),
        (0.01, 0, 0.009_901_473_843_595_012),
        (-0.2, -1, -2.542_641_357_773_526_5),
        (-1e-300, -1, -697.322_776_295_460_1),
        (-1e-320, -1, -743.438_526_972_854_4),
    ] {
        let w = xsf::lambertw(c64(z, 0.0), k, 1e-8);
        assert_eq!(w.im, 0.0);
        xsref::np_assert_allclose!([w.re], [expected], rtol = 1e-15);
    }
}
