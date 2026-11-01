use core::f64::consts::TAU;

const ATOL: f64 = 1.5e-15;

/// Translated from `scipy.special.tests.test_basic.TestSystematic.test_diric`
#[test]
fn test_diric() {
    // Test behavior near multiples of tau. Regression test for scipy/scipy#4001.

    // n_odd = [1, 5, 25]
    let n_odd = [1, 5, 25];
    // x = np.array(2*np.pi + 1e-9).astype(np.float64)
    let x = TAU + 1e-9;
    // assert_allclose(special.diric(x, n_odd), 1.0, atol=ATOL, rtol=0)
    xsref::np_assert_allclose!(n_odd.map(|n| xsf::diric(n, x)), &[1.0; 3], atol = ATOL);

    // x = np.array(2*np.pi + 1e-15).astype(np.float64)
    let x = TAU + 1e-15;
    // assert_allclose(special.diric(x, n_odd), 1.0, atol=ATOL, rtol=0)
    xsref::np_assert_allclose!(n_odd.map(|n| xsf::diric(n, x)), &[1.0; 3], atol = ATOL);

    // n_even = [2, 4, 24]
    let n_even = [2, 4, 24];
    // x = np.array(2*np.pi + 1e-9).astype(np.float64)
    let x = TAU + 1e-9;
    // assert_allclose(special.diric(x, n_even), -1.0, atol=ATOL, rtol=0)
    xsref::np_assert_allclose!(n_even.map(|n| xsf::diric(n, x)), &[-1.0; 3], atol = ATOL);

    // Test at some values not near a multiple of pi
    // x = np.arange(0.2*np.pi, 1.0*np.pi, 0.2*np.pi)
    let x = [0.1 * TAU, 0.2 * TAU, 0.3 * TAU, 0.4 * TAU];
    // octave_result = [0.872677996249965, 0.539344662916632,
    //                  0.127322003750035, -0.206011329583298]
    let octave_result = [
        0.872_677_996_249_965,
        0.539_344_662_916_632,
        0.127_322_003_750_035,
        -0.206_011_329_583_298,
    ];
    // assert_allclose(special.diric(x, 3), octave_result, atol=ATOL, rtol=0)
    xsref::np_assert_allclose!(x.map(|x| xsf::diric(3, x)), &octave_result, atol = ATOL);
}
