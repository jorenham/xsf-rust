const PASCAL_L8: [[usize; 8]; 8] = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 0, 0, 0],
    [1, 2, 1, 0, 0, 0, 0, 0],
    [1, 3, 3, 1, 0, 0, 0, 0],
    [1, 4, 6, 4, 1, 0, 0, 0],
    [1, 5, 10, 10, 5, 1, 0, 0],
    [1, 6, 15, 20, 15, 6, 1, 0],
    [1, 7, 21, 35, 35, 21, 7, 1],
];
const PASCAL_S08: [[usize; 8]; 9] = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 2, 3, 4, 5, 6, 7, 8],
    [1, 3, 6, 10, 15, 21, 28, 36],
    [1, 4, 10, 20, 35, 56, 84, 120],
    [1, 5, 15, 35, 70, 126, 210, 330],
    [1, 6, 21, 56, 126, 252, 462, 792],
    [1, 7, 28, 84, 210, 462, 924, 1716],
    [1, 8, 36, 120, 330, 792, 1716, 3432],
];

#[test]
fn test_comb() {
    for (n, nc) in PASCAL_L8.iter().enumerate() {
        for (k, &nck) in nc.iter().enumerate() {
            assert_eq!(xsf::comb(n, k), nck);
        }
    }
}

#[test]
fn test_comb_rep() {
    for (n, nr) in PASCAL_S08.iter().enumerate() {
        for (k, &nrk) in nr.iter().enumerate() {
            assert_eq!(xsf::comb_rep(n, k), nrk);
        }
    }
}

#[test]
fn test_comb_no_overflow() {
    // large results that fit in a u64, where the intermediate products used to overflow
    assert_eq!(xsf::comb(63_u64, 29), 759_510_004_936_100_355);
    assert_eq!(xsf::comb(67_u64, 33), 14_226_520_737_620_288_370);
    assert_eq!(xsf::comb_rep(34_u64, 34), 14_226_520_737_620_288_370);
}

#[test]
fn test_comb_max_n() {
    // results that fit, at or near `N::MAX`
    assert_eq!(xsf::comb(u64::MAX, 1), u64::MAX);
    assert_eq!(xsf::comb_rep(1_u8, 255), 1);
    assert_eq!(xsf::comb_rep(u8::MAX, 1), 255);
    assert_eq!(xsf::comb_rep(2_u8, 254), 255);
}

#[test]
fn test_comb_u8_exhaustive() {
    // Pascal's triangle, saturated at u16::MAX, against every `u8` result that fits
    let mut row = vec![1_u16];
    for n in 0..=u8::MAX {
        for (k, &c) in (0..=u8::MAX).zip(&row) {
            if let Ok(c) = u8::try_from(c) {
                assert_eq!(xsf::comb(n, k), c, "comb({n}, {k})");
            }
        }
        let mut next = vec![1];
        next.extend(row.windows(2).map(|w| w[0].saturating_add(w[1])));
        next.push(1);
        row = next;
    }
}
