#[test]
fn test_gdtrib() {
    xsref::test("gdtrib", "d_d_d-d", |x| xsf::gdtrib(x[0], x[1], x[2]));
}

#[test]
fn test_gdtria_gdtrix() {
    // reference values from scipy.special
    xsref::np_assert_allclose!([xsf::gdtria(0.5, 2.0, 3.0)], [0.559_448_996_672_220_4]);
    xsref::np_assert_allclose!([xsf::gdtrix(1.5, 2.5, 0.25)], [0.891_534_269_810_721_3]);
    assert_eq!(xsf::gdtria(0.0, 0.0, 1.0), 0.0);
    assert!(xsf::gdtria(0.0, 0.0, f64::INFINITY).is_nan());
    assert!(xsf::gdtria(0.5, 1.0, 0.0).is_nan());
    assert!(xsf::gdtrix(0.0, 0.0, 0.5).is_nan());
    assert!(xsf::gdtrix(f64::INFINITY, 1.0, 0.5).is_nan());
}
