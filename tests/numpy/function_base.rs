use xsref::np_assert_allclose;

#[test]
fn test_sinc() {
    // based on `numpy.lib.tests.test_function_base.TestSinc.test_simple`
    let w = (-1..=1)
        .map(|i| f64::from(i) * 0.01)
        .map(xsf::sinc)
        .collect::<Vec<f64>>();
    np_assert_allclose!(
        w.clone(),
        w.iter().copied().rev().collect::<Vec<_>>(),
        atol = 1e-15
    );
}

#[test]
fn test_sinc_large_or_near_int() {
    // mpmath reference values
    let x = [
        -9_781_324_477.118_902,
        1.000_000_000_1,
        2.999_999_75,
        1e15 + 0.5,
    ];
    let y = [
        -1.187_531_271_754_768_9e-11,
        -1.000_000_082_640_371e-10,
        8.333_334_021_540_307e-8,
        3.183_098_861_837_905e-16,
    ];
    np_assert_allclose!(x.map(xsf::sinc), y, rtol = 1e-15);
}
