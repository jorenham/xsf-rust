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
