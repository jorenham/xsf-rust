xsref::tests! {
    test_prolate_segv: "prolate_segv", "d_d_d-d", |x| {
        xsref::with_int(x[0], |m| xsref::with_int(x[1], |n| xsf::prolate_segv(m, n, x[2])))
    };

    // https://github.com/scipy/xsref/issues/9
    // test_oblate_segv: "oblate_segv", "d_d_d-d", |x| {
    //     xsf::oblate_segv(x[0] as u64, x[1] as u64, x[2])
    // };

    test_prolate_aswfa_nocv: "prolate_aswfa_nocv", "d_d_d_d-d_d", |x| {
        xsf::prolate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    };
    test_oblate_aswfa_nocv: "oblate_aswfa_nocv", "d_d_d_d-d_d", |x| {
        xsf::oblate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    };
    test_prolate_radial1_nocv: "prolate_radial1_nocv", "d_d_d_d-d_d", |x| {
        xsf::prolate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    };
    test_oblate_radial1_nocv: "oblate_radial1_nocv", "d_d_d_d-d_d", |x| {
        xsf::oblate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    };
    test_prolate_radial2_nocv: "prolate_radial2_nocv", "d_d_d_d-d_d", |x| {
        xsf::prolate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    };
    test_oblate_radial2_nocv: "oblate_radial2_nocv", "d_d_d_d-d_d", |x| {
        xsf::oblate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    };
    test_prolate_aswfa: "prolate_aswfa", "d_d_d_d_d-d_d", |x| {
        xsf::prolate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    };
    test_oblate_aswfa: "oblate_aswfa", "d_d_d_d_d-d_d", |x| {
        xsf::oblate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    };
    test_prolate_radial1: "prolate_radial1", "d_d_d_d_d-d_d", |x| {
        xsf::prolate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    };
    test_oblate_radial1: "oblate_radial1", "d_d_d_d_d-d_d", |x| {
        xsf::oblate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    };
    test_prolate_radial2: "prolate_radial2", "d_d_d_d_d-d_d", |x| {
        xsf::prolate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    };
    test_oblate_radial2: "oblate_radial2", "d_d_d_d_d-d_d", |x| {
        xsf::oblate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    };
}

#[test]
fn test_prolate_radial1_nocv_out_of_bounds() {
    // the 252-element buffers of `rmn1` would overflow
    let (r1f, r1d) = xsf::prolate_radial1_nocv(15, 174, 30.2, 1.077);
    assert!(r1f.is_nan() && r1d.is_nan());
}

#[test]
fn test_oblate_radial2_wronskian() {
    // `R1 R2' - R1' R2 = 1 / (c (x^2 + 1))`
    for (m, n, c, x) in [
        (0, 0, 1.0, 0.0),
        (0, 0, 1.0, 0.5),
        (1, 1, 1.0, 0.0),
        (0, 1, 1.0, 0.5),
        (1, 1, 1.0, 0.5),
        (1, 2, 1.0, 0.5),
        (2, 3, 2.0, 0.2),
        (3, 6, 0.7, 0.05),
        (5, 8, 3.0, 0.4),
    ] {
        let (r1f, r1d) = xsf::oblate_radial1_nocv(m, n, c, x);
        let (r2f, r2d) = xsf::oblate_radial2_nocv(m, n, c, x);
        let w = (r1f * r2d - r1d * r2f) * c * (x * x + 1.0);
        assert!((w - 1.0).abs() < 1e-13, "{m} {n} {c} {x}: {w}");

        let cv = xsf::oblate_segv(m, n, c);
        assert_eq!(xsf::oblate_radial2(m, n, c, cv, x), (r2f, r2d));
    }
}

#[test]
fn test_prolate_segv_c0_large_n() {
    // `n (n + 1)`, without `int` overflow
    assert_eq!(xsf::prolate_segv(46341, 46341, 0.0), 46341.0 * 46342.0);
}

type Case = (u64, u64, f64, &'static [f64], &'static [f64]);

/// Compares with `scipy.special.{pro,obl}_cv_seq` (SciPy 1.18.1)
#[test]
fn test_segv_seq_scipy() {
    #[rustfmt::skip]
    let cases: [Case; 4] = [
        (0, 5, 1.0,
            &[0.319_000_055_146_892_84, 2.593_084_579_977_15, 6.533_471_800_523_812,
              12.514_462_145_094_072, 20.508_274_362_570_884, 30.505_404_625_322_107],
            &[-0.348_602_399_470_268_2, 1.393_206_310_448_417, 5.486_800_053_818_691,
              11.492_120_902_297_92, 19.495_276_746_591_117, 29.496_855_283_297_478]),
        (2, 6, 3.5,
            &[7.513_152_060_394_042_5, 15.770_798_673_512_214, 24.854_790_768_938_39,
              35.358_649_669_244_656, 47.612_636_237_010_71],
            &[3.877_796_220_360_97, 7.558_328_003_098_823, 15.107_696_662_176_464,
              24.716_458_501_131_573, 36.490_828_227_858_87]),
        (1, 3, 50.0,
            &[50.256_464_865_030_09, 149.262_155_810_652_98, 247.245_314_587_757_38],
            &[-2_302.020_626_337_337_5, -2_302.020_626_337_342, -2_108.169_051_280_062_4]),
        (0, 3, 0.0, &[0.0, 2.0, 6.0, 12.0], &[0.0, 2.0, 6.0, 12.0]),
    ];
    for (m, n, c, pro, obl) in cases {
        for (name, got, want, last) in [
            (
                "pro",
                xsf::prolate_segv_seq(m, n, c),
                pro,
                xsf::prolate_segv(m, n, c),
            ),
            (
                "obl",
                xsf::oblate_segv_seq(m, n, c),
                obl,
                xsf::oblate_segv(m, n, c),
            ),
        ] {
            assert_eq!(got.len(), want.len());
            assert_eq!(got.last(), Some(&last));
            for (i, (g, w)) in got.iter().zip(want).enumerate() {
                let tol = 1e-12 * w.abs().max(1e-3);
                assert!((g - w).abs() <= tol, "{name} {m} {n} {c} [{i}]: {g} vs {w}");
            }
        }
    }
}

#[test]
fn test_segv_seq_longest() {
    // `n - m = 199` is the largest difference accepted by SciPy
    for (kd, seq) in [
        (1, xsf::prolate_segv_seq(0, 199, 1.0)),
        (-1, xsf::oblate_segv_seq(0, 199, 1.0)),
    ] {
        assert_eq!(seq.len(), 200);
        let want = if kd == 1 {
            [
                0.319_000_055_146_892_56,
                10_100.500_015_471_911,
                39_800.500_003_925_874,
            ]
        } else {
            [
                -0.348_602_399_470_269_7,
                10_099.499_990_717_593,
                39_799.499_997_644_336,
            ]
        };
        for (g, w) in [seq[0], seq[100], seq[199]].into_iter().zip(want) {
            assert!((g - w).abs() <= 1e-12 * w.abs(), "{kd}: {g} vs {w}");
        }
    }
}

#[test]
fn test_segv_seq_nan() {
    // m and n round to the same f64
    let seq = xsf::prolate_segv_seq(1 << 53, (1 << 53) + 1, 0.0);
    assert!(seq.len() == 2 && seq.iter().all(|x| x.is_nan()));
    // SciPy segfaults for `c` this large
    assert!(xsf::oblate_segv_seq(1, 3, 400.0).iter().all(|x| x.is_nan()));
}

#[test]
#[should_panic(expected = "n - m must be in 0..=199")]
fn test_segv_seq_too_long() {
    let _ = xsf::prolate_segv_seq(0, 200, 1.0);
}

#[test]
#[should_panic(expected = "n - m must be in 0..=199")]
fn test_segv_seq_n_lt_m() {
    let _ = xsf::oblate_segv_seq(3, 1, 1.0);
}
