// Ported from `scipy.special.tests.test_gamma` and `test_basic`

const INF: f64 = f64::INFINITY;

#[test]
fn test_log_gammainc_gammaincc() {
    // consistent with `gammainc` and `gammaincc` away from the tails
    for (a, x) in [(1.0, 2.0), (2.0, 1.0)] {
        xsref::np_assert_allclose!(
            &[xsf::log_gammainc(a, x), xsf::log_gammaincc(a, x)],
            &[xsf::gammainc(a, x).ln(), xsf::gammaincc(a, x).ln()],
            rtol = 1e-15
        );
    }
}

#[test]
fn test_log_gammainc_edge_cases() {
    for (a, x) in [(0.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (0.0, 0.0), (INF, INF)] {
        assert!(xsf::log_gammainc(a, x).is_nan());
        assert!(xsf::log_gammaincc(a, x).is_nan());
    }
    for a in [0.5, 1.0, 2.0, 100.0] {
        assert_eq!(xsf::log_gammainc(a, 0.0), -INF);
        assert_eq!(xsf::log_gammaincc(a, 0.0), 0.0);
        assert_eq!(xsf::log_gammainc(a, INF), 0.0);
        assert_eq!(xsf::log_gammaincc(a, INF), -INF);
    }
    for x in [0.5, 1.0, 2.0, 100.0] {
        assert_eq!(xsf::log_gammainc(INF, x), -INF);
        assert_eq!(xsf::log_gammaincc(INF, x), 0.0);
    }
}

#[test]
fn test_log_gammainc_tails() {
    // `(a, x, ln P, ln Q)` from Arb, also where `P` or `Q` underflows, or rounds to 1
    for (a, x, ln_p, ln_q) in [
        (0.5, 5e-101, -115.355_046_002_347_01, -7.978_845_608_028_654e-51),
        (0.5, 750.0, -0.0, -753.883_067_105_382_5),
        (1e-300, 700.0, -0.0, -1_397.328_033_756_021_3),
        (1e-310, 1e-10, -2.244_863_526_513_885e-309, -710.690_149_005_785_2),
        (1e-310, 700.0, -0.0, -1_420.353_884_685_962),
        (1e-320, 1e-100, -2.296_79e-318, -731.390_548_222_842_9),
        (1.0, 1e200, -0.0, -1e200),
        (5.0, 1e-300, -3_458.665_131_233_850_4, -0.0),
        (10.0, 1e-100, -2_317.689_505_567_121, -0.0),
        (10.0, 700.0, -1.110_692_465_019_363_7e-284, -653.829_182_747_133_7),
        (10.0, 1e5, -0.0, -99_909.185_408_292_2),
        (100.0, 10.0, -143.376_723_100_618_88, -5.398_589_728_139_581e-63),
        (1e6, 9e5, -5_366.039_856_607_687, -0.0),
        (1e9, 9.9e8, -50_342.528_912_606_71, -0.0),
        (1e9, 1.01e9, -0.0, -49_675.822_243_282_67),
        (6.666_666_666_666_667e199, 1e200, -0.0, -6.302_326_126_122_374e198),
    ] {
        xsref::np_assert_allclose!(
            &[xsf::log_gammainc(a, x), xsf::log_gammaincc(a, x)],
            &[ln_p, ln_q],
            rtol = 1e-15
        );
    }
}
