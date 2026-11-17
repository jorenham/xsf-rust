use num_complex::c64;

xsref::tests! {
    test_ber_f64: "ber", "d-d", |x| xsf::ber(x[0]);
    test_bei_f64: "bei", "d-d", |x| xsf::bei(x[0]);
    test_ker_f64: "ker", "d-d", |x| xsf::ker(x[0]);
    test_kei_f64: "kei", "d-d", |x| xsf::kei(x[0]);
    test_berp_f64: "berp", "d-d", |x| xsf::berp(x[0]);
    test_beip_f64: "beip", "d-d", |x| xsf::beip(x[0]);
    test_kerp_f64: "kerp", "d-d", |x| xsf::kerp(x[0]);
    test_keip_f64: "keip", "d-d", |x| xsf::keip(x[0]);
    test_kelvin_c64: "kelvin", "d-cd_cd_cd_cd", |x| {
        let [a, b, c, d] = xsf::kelvin(x[0]);
        (a, b, c, d)
    };
}

/// Ported from xsf's `test_kelvin.cpp`
#[test]
fn test_kelvin_wronskian() {
    // x (Be Ke' - Be' Ke) = -1 checks all eight functions at once.
    //
    // Unlike in xsf (scipy/xsf#232), it holds to 1e-14 (xsf: 2.5e-10 at x = 9.9).
    for x in [0.01, 1.0, 9.9, 10.1, 30.0] {
        let [be, ke, bep, kep] = xsf::kelvin(x);
        let residual = (x * (be * kep - bep * ke) + 1.0).norm();
        assert!(residual <= 1e-14, "x = {x}: residual = {residual:e}");
    }
}

#[test]
fn test_kelvin_be_mid() {
    // ber + i bei and ber' + i bei' for 10 <= x < 18, computed with mpmath (40 digits)
    let cases = [
        (
            10.0,
            (138.840_465_941_632_66, 56.370_458_553_906_64),
            (51.195_258_393_599_61, 135.309_301_715_964_6),
        ),
        (
            14.0,
            (-2_131.281_170_574_877_7, -160.937_689_652_486_43),
            (-1_316.093_925_831_294_4, -1_616.089_442_962_944),
        ),
        (
            17.5,
            (18_848.671_223_081_736, -12_619.261_010_554_887),
            (21_711.202_282_744_056, 4_774.935_372_065_795),
        ),
    ];
    for (x, (ber, bei), (berp, beip)) in cases {
        let [be, _, bep, _] = xsf::kelvin(x);
        for (y, r) in [(be, c64(ber, bei)), (bep, c64(berp, beip))] {
            let err = (y - r).norm() / r.norm();
            assert!(err <= 1e-14, "x = {x}: {y} != {r} (err = {err:e})");
        }
    }
}

#[test]
fn test_kelvin_ke_mid() {
    // ker + i kei and ker' + i kei' for 3 <= x < 18, computed with mpmath (40 digits)
    let cases = [
        (
            4.0,
            (-0.036_178_847_899_547_61, 0.002_198_399_294_972_519_7),
            (0.031_478_489_822_108_15, 0.023_910_613_791_704_644),
        ),
        (
            9.5,
            (0.000_333_028_933_957_937_05, -0.000_355_743_339_928_496_86),
            (-0.000_504_543_773_915_663_6, 3.415_790_306_469_979e-5),
        ),
        (
            12.5,
            (-5.002_032_136_217_098_3e-5, -1.013_641_248_284_213_7e-5),
            (3.016_978_394_412_561_7e-5, 4.296_237_309_176_477e-5),
        ),
    ];
    for (x, (ker, kei), (kerp, keip)) in cases {
        let [_, ke, _, kep] = xsf::kelvin(x);
        for (y, r) in [(ke, c64(ker, kei)), (kep, c64(kerp, keip))] {
            let err = (y - r).norm() / r.norm();
            assert!(err <= 1e-14, "x = {x}: {y} != {r} (err = {err:e})");
        }
    }
}

/// Ported from `scipy.special.tests.test_basic.TestKelvin.test_*_zeros`
#[test]
fn test_kelvin_zeros() {
    let funcs: [fn(usize) -> Vec<f64>; 8] = [
        xsf::ber_zeros,
        xsf::bei_zeros,
        xsf::ker_zeros,
        xsf::kei_zeros,
        xsf::berp_zeros,
        xsf::beip_zeros,
        xsf::kerp_zeros,
        xsf::keip_zeros,
    ];
    // A&S 9.9, p. 381; the precise beip row is more accurate than the 3.77320 of the 1927 table
    let zeros = [
        [2.84892, 7.23883, 11.67396, 16.11356, 20.55463],
        [5.02622, 9.45541, 13.89349, 18.33398, 22.77544],
        [1.71854, 6.12728, 10.56294, 15.00269, 19.44382],
        [3.91467, 8.34422, 12.78256, 17.22314, 21.66464],
        [6.03871, 10.51364, 14.96844, 19.41758, 23.86430],
        [
            3.772_673_304_934_953,
            8.280_987_849_760_042,
            12.742_147_523_633_703,
            17.193_431_752_512_54,
            21.641_143_941_167_325,
        ],
        [2.66584, 7.17212, 11.63218, 16.08312, 20.53068],
        [4.93181, 9.40405, 13.85827, 18.30717, 22.75379],
    ];
    for ((f, expected), all) in funcs.iter().zip(&zeros).zip(xsf::kelvin_zeros(5)) {
        xsref::np_assert_allclose!(&f(5), expected, atol = 1.5e-4);
        xsref::np_assert_allclose!(&all, expected, atol = 1.5e-4);
    }
    xsref::np_assert_allclose!(&xsf::beip_zeros(5), &zeros[5], atol = 1.5e-8);
}
