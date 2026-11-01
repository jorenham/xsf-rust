//! Translated from `scipy.special.tests.test_logsumexp` at
//! <https://github.com/scipy/scipy/blob/5a7df53/scipy/special/tests/test_logsumexp.py>

use core::f64::consts::{E, LN_2};
use xsref::np_assert_allclose;

/// Corresponds to `TestLogSumExp.test_logsumexp`.
#[test]
#[allow(clippy::float_cmp, clippy::cast_precision_loss)]
fn test_logsumexp() {
    // default rtol of `xp_assert_close`
    let rtol: f64 = f64::EPSILON.sqrt() * 4.0;

    // Test with zero-size array
    assert_eq!(xsf::logsumexp(&[]), f64::NEG_INFINITY);

    // Test whether logsumexp() function correctly handles large inputs.
    let xs = (0..200).map(f64::from).collect::<Vec<f64>>();
    let desired = xs.iter().copied().map(f64::exp).sum::<f64>().ln();
    np_assert_allclose!([xsf::logsumexp(&xs)], [desired], rtol = rtol);

    // Now test with large numbers
    let xs = [1000.0, 1000.0];
    let desired = 1000.0 + LN_2;
    np_assert_allclose!([xsf::logsumexp(&xs)], [desired], rtol = rtol);

    let (x, n) = (10_000.0, 1_000);
    let xs = vec![x; n];
    let desired = x + (n as f64).ln();
    np_assert_allclose!([xsf::logsumexp(&xs)], [desired], rtol = rtol);

    let ys = vec![1e-40; 1_000_000];
    let xs = ys.iter().copied().map(f64::ln).collect::<Vec<f64>>();
    np_assert_allclose!(
        [xsf::logsumexp(&xs).exp()],
        [ys.iter().sum::<f64>()],
        rtol = rtol
    );

    // Handling special values properly
    assert!(xsf::logsumexp(&[f64::NAN]).is_nan());
    assert_eq!(xsf::logsumexp(&[f64::INFINITY]), f64::INFINITY);
    assert_eq!(xsf::logsumexp(&[f64::NEG_INFINITY]), f64::NEG_INFINITY);
    assert_eq!(
        xsf::logsumexp(&[f64::NEG_INFINITY, f64::NEG_INFINITY]),
        f64::NEG_INFINITY
    );
}

/// Corresponds to `TestSoftmax.test_softmax_fixtures`.
#[test]
fn test_softmax() {
    np_assert_allclose!(
        xsf::softmax(&[1000.0, 0.0, 0.0, 0.0]),
        [1.0, 0.0, 0.0, 0.0],
        rtol = 1e-13
    );
    np_assert_allclose!(xsf::softmax(&[1.0, 1.0]), [0.5, 0.5], rtol = 1e-13);
    np_assert_allclose!(
        xsf::softmax(&[0.0, 1.0]),
        [1.0 / (1.0 + E), E / (1.0 + E)],
        rtol = 1e-13
    );

    // Expected value computed using mpmath (with mpmath.mp.dps = 200) and then
    // converted to float.
    let x = [0.0, 1.0, 2.0, 3.0];
    let expected = [
        0.032_058_603_280_084_99,
        0.087_144_318_742_032_56,
        0.236_882_818_089_910_13,
        0.643_914_259_887_972_2,
    ];
    np_assert_allclose!(xsf::softmax(&x), expected, rtol = 1e-13);

    // Translation property.  If all the values are changed by the same amount,
    // the softmax result does not change.
    np_assert_allclose!(xsf::softmax(&x.map(|x| x + 100.0)), expected, rtol = 1e-13);
}

/// Corresponds to `TestLogSoftmax.*`
#[test]
fn test_log_softmax() {
    // test_log_softmax_basic
    np_assert_allclose!(xsf::log_softmax(&[1000.0, 1.0]), [0., -999.], rtol = 1e-13);

    // test_log_softmax_scalar
    np_assert_allclose!(xsf::log_softmax(&[1.0]), [0.], rtol = 1e-13);

    // test_log_softmax_translation
    let x = [0.0, 1.0, 2.0, 3.0];
    // Expected value computed using mpmath (with mpmath.mp.dps = 200)
    let expect = [
        -3.440_189_698_561_195_3,
        -2.440_189_698_561_195_3,
        -1.440_189_698_561_195_3,
        -0.440_189_698_561_195_33,
    ];
    np_assert_allclose!(xsf::log_softmax(&x), expect, rtol = 1e-13);

    // Translation property.  If all the values are changed by the same amount,
    // the softmax result does not change.
    np_assert_allclose!(
        xsf::log_softmax(&x.map(|x| x + 100.0)),
        expect,
        rtol = 1e-13
    );
}
