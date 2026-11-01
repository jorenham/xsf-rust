const ATOL: f64 = 1.5e-10;

/// Translated from `scipy.special.tests.test_basic.TestBessel.test_jvp`
#[test]
fn test_bessel_j_prime() {
    // jvprim = special.jvp(2,2)
    let jvprim = xsf::bessel_j_prime(2.0, 2.0, 1);
    // jv0 = (special.jv(1,2)-special.jv(3,2))/2
    let jv0 = (xsf::bessel_j(1.0, 2.0) - xsf::bessel_j(3.0, 2.0)) * 0.5;
    // assert_allclose(jvprim, jv0, atol=1.5e-10, rtol=0)
    xsref::np_assert_allclose!(&[jvprim], &[jv0], atol = ATOL);
}

/// Translated from `scipy.special.tests.test_basic.TestBessel.test_yvp`
#[test]
#[allow(clippy::similar_names)]
fn test_bessel_y_prime() {
    // yvpr = (special.yv(1,.2) - special.yv(3,.2))/2.0
    let yvpr = (xsf::bessel_y(1.0, 0.2) - xsf::bessel_y(3.0, 0.2)) * 0.5;
    // yvp1 = special.yvp(2,.2)
    let yvp1 = xsf::bessel_y_prime(2.0, 0.2, 1);
    // assert_allclose(yvp1, yvpr, atol=1.5e-10, rtol=0)
    xsref::np_assert_allclose!(&[yvp1], &[yvpr], atol = ATOL);
}

/// Translated from `scipy.special.tests.test_basic.TestBessel.test_ivp0`
#[test]
fn test_bessel_i_prime_0() {
    // assert_allclose(special.iv(1, 2), special.ivp(0, 2), atol=1.5e-10, rtol=0)
    xsref::np_assert_allclose!(
        &[xsf::bessel_i(1.0, 2.0)],
        &[xsf::bessel_i_prime(0.0, 2.0, 1)],
        atol = ATOL
    );
}

/// Translated from `scipy.special.tests.test_basic.TestBessel.test_ivp`
#[test]
fn test_bessel_i_prime() {
    // y = (special.iv(0,2) + special.iv(2,2))/2
    let y = f64::midpoint(xsf::bessel_i(0.0, 2.0), xsf::bessel_i(2.0, 2.0));
    // x = special.ivp(1,2)
    let x = xsf::bessel_i_prime(1.0, 2.0, 1);
    // assert_allclose(x, y, atol=1.5e-10, rtol=0)
    xsref::np_assert_allclose!(&[x], &[y], atol = ATOL);
}

/// Translated from `scipy.special.tests.test_basic.TestBessel.test_kvp_v0n1`
#[test]
fn test_bessel_k_prime_v0n1() {
    // z = 2.2
    const Z: f64 = 2.2;
    // assert_allclose(-special.kv(1, z), special.kvp(0, z, n=1), atol=1.5e-10, rtol=0)
    xsref::np_assert_allclose!(
        &[-xsf::bessel_k(1.0, Z)],
        &[xsf::bessel_k_prime(0.0, Z, 1)],
        atol = ATOL
    );
}

/// Translated from `scipy.special.tests.test_basic.TestBessel.test_kvp_n1`
#[test]
fn test_bessel_k_prime_n1() {
    // v = 3.
    const V: f64 = 3.0;
    // z = 2.2
    const Z: f64 = 2.2;
    // xc = -special.kv(v+1,z) + v/z*special.kv(v,z)
    let xc = -xsf::bessel_k(V + 1.0, Z) + V / Z * xsf::bessel_k(V, Z);
    // x = special.kvp(v,z, n=1)
    let x = xsf::bessel_k_prime(V, Z, 1);
    // assert_allclose(xc, x, atol=1.5e-10, rtol=0)
    xsref::np_assert_allclose!(&[x], &[xc], atol = ATOL);
}

/// Translated from `scipy.special.tests.test_basic.TestBessel.test_kvp_n2`
#[test]
fn test_bessel_k_prime_n2() {
    // v = 3.
    const V: f64 = 3.0;
    // z = 2.2
    const Z: f64 = 2.2;
    // xc = (z**2 + v**2 - v) / z**2 * special.kv(v, z) + special.kv(v+1, z) / z
    const Z2: f64 = Z * Z;
    let xc = (Z2 + V * V - V) / Z2 * xsf::bessel_k(V, Z) + xsf::bessel_k(V + 1.0, Z) / Z;
    // x = special.kvp(v, z, n=2)
    let x = xsf::bessel_k_prime(V, Z, 2);
    // assert_allclose(xc, x, atol=1.5e-10, rtol=0)
    xsref::np_assert_allclose!(&[x], &[xc], atol = ATOL);
}
