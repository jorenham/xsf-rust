#![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

#[test]
fn test_prolate_segv() {
    xsref::test("prolate_segv", "d_d_d-d", |x| {
        xsf::prolate_segv(x[0] as u64, x[1] as u64, x[2])
    });
}

// https://github.com/scipy/xsref/issues/9
// #[test]
// fn test_oblate_segv() {
//     testing::test("oblate_segv", "d_d_d-d", |x| {
//         xsf::oblate_segv(x[0] as u64, x[1] as u64, x[2])
//     });
// }

#[test]
fn test_prolate_aswfa_nocv() {
    xsref::test("prolate_aswfa_nocv", "d_d_d_d-d_d", |x| {
        xsf::prolate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_oblate_aswfa_nocv() {
    xsref::test("oblate_aswfa_nocv", "d_d_d_d-d_d", |x| {
        xsf::oblate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_prolate_radial1_nocv() {
    xsref::test("prolate_radial1_nocv", "d_d_d_d-d_d", |x| {
        xsf::prolate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_oblate_radial1_nocv() {
    xsref::test("oblate_radial1_nocv", "d_d_d_d-d_d", |x| {
        xsf::oblate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_prolate_radial2_nocv() {
    xsref::test("prolate_radial2_nocv", "d_d_d_d-d_d", |x| {
        xsf::prolate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_prolate_aswfa() {
    xsref::test("prolate_aswfa", "d_d_d_d_d-d_d", |x| {
        xsf::prolate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_oblate_aswfa() {
    xsref::test("oblate_aswfa", "d_d_d_d_d-d_d", |x| {
        xsf::oblate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_prolate_radial1() {
    xsref::test("prolate_radial1", "d_d_d_d_d-d_d", |x| {
        xsf::prolate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_oblate_radial1() {
    xsref::test("oblate_radial1", "d_d_d_d_d-d_d", |x| {
        xsf::oblate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_prolate_radial2() {
    xsref::test("prolate_radial2", "d_d_d_d_d-d_d", |x| {
        xsf::prolate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_prolate_radial1_nocv_out_of_bounds() {
    // the 252-element buffers of `rmn1` would overflow
    let (r1f, r1d) = xsf::prolate_radial1_nocv(15, 174, 30.2, 1.077);
    assert!(r1f.is_nan() && r1d.is_nan());
}

// The xsref `oblate_radial2` cases are computed with scipy, which has the `qstar` and `gmn`
// bugs that are fixed in xsf-rust.
#[test]
fn test_oblate_radial2_wronskian() {
    // `R1 R2' - R1' R2 = 1 / (c (x^2 + 1))`
    for (m, n, c, x) in [
        (0, 0, 1.0, 0.0),
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
