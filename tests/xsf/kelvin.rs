#[test]
fn test_ber_f64() {
    xsref::test("ber", "d-d", |x| xsf::ber(x[0]));
}

#[test]
fn test_bei_f64() {
    xsref::test("bei", "d-d", |x| xsf::bei(x[0]));
}

#[test]
fn test_ker_f64() {
    xsref::test("ker", "d-d", |x| xsf::ker(x[0]));
}

#[test]
fn test_kei_f64() {
    xsref::test("kei", "d-d", |x| xsf::kei(x[0]));
}

#[test]
fn test_berp_f64() {
    xsref::test("berp", "d-d", |x| xsf::berp(x[0]));
}

#[test]
fn test_beip_f64() {
    xsref::test("beip", "d-d", |x| xsf::beip(x[0]));
}

#[test]
fn test_kerp_f64() {
    xsref::test("kerp", "d-d", |x| xsf::kerp(x[0]));
}

#[test]
fn test_keip_f64() {
    xsref::test("keip", "d-d", |x| xsf::keip(x[0]));
}

#[test]
fn test_kelvin_c64() {
    xsref::test("kelvin", "d-cd_cd_cd_cd", |x| {
        let [a, b, c, d] = xsf::kelvin(x[0]);
        (a, b, c, d)
    });
}

// Kelvin function zeros tests

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_ber_zeros`
#[test]
fn test_ber_zeros() {
    let ber = xsf::ber_zeros(5);
    xsref::np_assert_allclose!(
        &ber,
        &[2.84892, 7.23883, 11.67396, 16.11356, 20.55463],
        atol = 1.5e-4
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_bei_zeros`
#[test]
fn test_bei_zeros() {
    // Abramowitz & Stegun, Table 9.12
    let bi = xsf::bei_zeros(5);
    xsref::np_assert_allclose!(
        &bi,
        &[5.02622, 9.45541, 13.89349, 18.33398, 22.77544],
        atol = 1.5e-4
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_ker_zeros`
#[test]
fn test_ker_zeros() {
    let ker = xsf::ker_zeros(5);
    xsref::np_assert_allclose!(
        &ker,
        &[1.71854, 6.12728, 10.56294, 15.00269, 19.44381],
        atol = 1.5e-4
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_kei_zeros`
#[test]
fn test_kei_zeros() {
    let kei = xsf::kei_zeros(5);
    xsref::np_assert_allclose!(
        &kei,
        &[3.91467, 8.34422, 12.78256, 17.22314, 21.66464],
        atol = 1.5e-4
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_berp_zeros`
#[test]
fn test_berp_zeros() {
    let brp = xsf::berp_zeros(5);
    xsref::np_assert_allclose!(
        &brp,
        &[6.03871, 10.51364, 14.96844, 19.41758, 23.86430],
        atol = 1.5e-4
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_beip_zeros`
#[test]
fn test_beip_zeros() {
    let bip = xsf::beip_zeros(5);
    xsref::np_assert_allclose!(
        &bip,
        &[
            3.772_673_304_934_953,
            8.280_987_849_760_042,
            12.742_147_523_633_703,
            17.193_431_752_512_54,
            21.641_143_941_167_325,
        ],
        atol = 1.5e-8
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_kerp_zeros`
#[test]
fn test_kerp_zeros() {
    let kerp = xsf::kerp_zeros(5);
    xsref::np_assert_allclose!(
        &kerp,
        &[2.66584, 7.17212, 11.63218, 16.08312, 20.53068],
        atol = 1.5e-4
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_keip_zeros`
#[test]
fn test_keip_zeros() {
    let keip = xsf::keip_zeros(5);
    xsref::np_assert_allclose!(
        &keip,
        &[4.93181, 9.40405, 13.85827, 18.30717, 22.75379],
        atol = 1.5e-4
    );
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_kelvin_zeros`
#[test]
fn test_kelvin_zeros() {
    // numbers come from 9.9 of A&S pg. 381
    let tmp = xsf::kelvin_zeros(5);
    let [berz, beiz, kerz, keiz, berpz, beipz, kerpz, keipz] = tmp;

    xsref::np_assert_allclose!(
        &berz,
        &[2.84892, 7.23883, 11.67396, 16.11356, 20.55463],
        atol = 1.5e-4
    );
    xsref::np_assert_allclose!(
        &beiz,
        &[5.02622, 9.45541, 13.89349, 18.33398, 22.77544],
        atol = 1.5e-4
    );
    xsref::np_assert_allclose!(
        &kerz,
        &[1.71854, 6.12728, 10.56294, 15.00269, 19.44382],
        atol = 1.5e-4
    );
    xsref::np_assert_allclose!(
        &keiz,
        &[3.91467, 8.34422, 12.78256, 17.22314, 21.66464],
        atol = 1.5e-4
    );
    xsref::np_assert_allclose!(
        &berpz,
        &[6.03871, 10.51364, 14.96844, 19.41758, 23.86430],
        atol = 1.5e-4
    );
    xsref::np_assert_allclose!(
        &beipz,
        // table from 1927 had 3.77320 but this is more accurate
        &[3.77267, 8.28099, 12.74215, 17.19343, 21.64114],
        atol = 1.5e-4
    );
    xsref::np_assert_allclose!(
        &kerpz,
        &[2.66584, 7.17212, 11.63218, 16.08312, 20.53068],
        atol = 1.5e-4
    );
    xsref::np_assert_allclose!(
        &keipz,
        &[4.93181, 9.40405, 13.85827, 18.30717, 22.75379],
        atol = 1.5e-4
    );
}
