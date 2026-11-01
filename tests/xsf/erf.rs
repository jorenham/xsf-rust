use num_complex::c64;

// erf

#[test]
fn test_erf_f64() {
    xsref::test("erf", "d-d", |x| xsf::erf(x[0]));
}

#[test]
fn test_erf_c64() {
    xsref::test("erf", "cd-cd", |x| xsf::erf(c64(x[0], x[1])));
}

// erfc

#[test]
fn test_erfc_f64() {
    xsref::test("erfc", "d-d", |x| xsf::erfc(x[0]));
}

#[test]
fn test_erfc_c64() {
    xsref::test("erfc", "cd-cd", |x| xsf::erfc(c64(x[0], x[1])));
}

// erfcx

#[test]
fn test_erfcx_f64() {
    xsref::test("erfcx", "d-d", |x| xsf::erfcx(x[0]));
}

#[test]
fn test_erfcx_c64() {
    xsref::test("erfcx", "cd-cd", |x| xsf::erfcx(c64(x[0], x[1])));
}

// erfi

#[test]
fn test_erfi_f64() {
    xsref::test("erfi", "d-d", |x| xsf::erfi(x[0]));
}

#[test]
fn test_erfi_c64() {
    xsref::test("erfi", "cd-cd", |x| xsf::erfi(c64(x[0], x[1])));
}

// erf_zeros (from `scipy.speceial.tests.test_basic.TestErf.test_erf_zeros`)

#[test]
fn test_erf_zeros() {
    let erz = xsf::erf_zeros::<5>();
    let erzr = [
        c64(1.450_616_16, 1.880_943_00),
        c64(2.244_659_28, 2.616_575_14),
        c64(2.839_741_05, 3.175_628_10),
        c64(3.335_460_74, 3.646_174_38),
        c64(3.769_005_57, 4.060_697_23),
    ];
    xsref::np_assert_allclose!(&erz, &erzr, atol = 1.5e-4);
}

// dawsn

#[test]
fn test_dawsn_f64() {
    xsref::test("dawsn", "d-d", |x| xsf::dawsn(x[0]));
}

#[test]
fn test_dawsn_c64() {
    xsref::test("dawsn", "cd-cd", |x| xsf::dawsn(c64(x[0], x[1])));
}

// wofz

#[test]
fn test_wofz() {
    xsref::test("wofz", "cd-cd", |x| xsf::wofz(c64(x[0], x[1])));
}

// voigt_profile

#[test]
fn test_voigt_profile() {
    xsref::test("voigt_profile", "d_d_d-d", |x| {
        xsf::voigt_profile(x[0], x[1], x[2])
    });
}
