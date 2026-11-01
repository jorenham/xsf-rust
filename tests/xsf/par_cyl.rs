#[test]
fn test_pbwa() {
    xsref::test("pbwa", "d_d-d_d", |x| xsf::pbwa(x[0], x[1]));
}

#[test]
fn test_pbdv() {
    xsref::test("pbdv", "d_d-d_d", |x| xsf::pbdv(x[0], x[1]));
}

#[test]
fn test_pbdv_negative_x() {
    // non-integer v > 0 and -10 <= x < -5.8 (reference values from mpmath)
    for (v, x, d, dp) in [
        (
            7.821_919_754_766_831_5,
            -6.126_185_993_244_382,
            142.563_617_884_111_64,
            -29.698_107_853_800_717,
        ),
        (
            0.373_152_185_433_848_87,
            -5.814_209_871_190_306,
            -287.813_586_230_075_5,
            762.841_520_864_751_4,
        ),
        (9.5, -8.0, 9_081.000_904_297_447, -20_350.620_138_313_516),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        let rel = ((y - d) / d).abs().max(((yp - dp) / dp).abs());
        assert!(rel < 1e-13, "v = {v}, x = {x}: rel = {rel:e}");
    }

    // near an integer v, or v >= x^2 / 3, where the series is not used
    for (v, x, d, dp, tol) in [
        (
            1e-20,
            -8.0,
            1.125_351_464_193_068_7e-7,
            4.501_408_120_768_454_4e-7,
            1e-6,
        ),
        (
            37.999_998_999_999,
            -9.0,
            -9.676_543_736_549_405e21,
            -7.936_223_618_237_496e21,
            1e-2,
        ),
    ] {
        let (y, yp) = xsf::pbdv(v, x);
        let rel = ((y - d) / d).abs().max(((yp - dp) / dp).abs());
        assert!(rel < tol, "v = {v}, x = {x}: rel = {rel:e}");
    }
}

#[test]
fn test_pbvv() {
    xsref::test("pbvv", "d_d-d_d", |x| xsf::pbvv(x[0], x[1]));
}

#[test]
fn test_pbvv_extreme() {
    // v >= 0 and 7.5 < |x| <= 10, and v near an integer (reference values from mpmath)
    for (v, x, w, wp) in [
        (
            8.463_751_277_865_178,
            7.641_451_046_470_024,
            0.023_227_491_821_967_712,
            0.050_155_659_621_339_22,
        ),
        (
            9.954_321_571_892_986,
            -7.634_719_796_831_582,
            -0.001_729_313_216_016_740_6,
            0.002_985_548_082_122_786,
        ),
        (6.0, -8.0, -5.663_320_431_230_510_5, 16.810_707_783_832_23),
        // near an integer v
        (
            2.000_000_01,
            7.025_168_798_105_803_5,
            601.282_319_691_077_8,
            1_829.270_702_453_187_4,
        ),
    ] {
        let (y, yp) = xsf::pbvv(v, x);
        let rel = ((y - w) / w).abs().max(((yp - wp) / wp).abs());
        assert!(rel < 1e-12, "v = {v}, x = {x}: rel = {rel:e}");
    }
}
