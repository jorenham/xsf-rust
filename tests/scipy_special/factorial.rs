#![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

const FACT1_5: u8 = 120;
const FACT1_8: u16 = 40_320;
const FACT1_12: u32 = 479_001_600;
const FACT1_20: u64 = 2_432_902_008_176_640_000;
const FACT1_34: u128 = 295_232_799_039_604_140_847_618_609_643_520_000_000;

#[test]
fn test_factorial_one() {
    for n in [0, 1] {
        assert_eq!(xsf::factorial(n as u8), 1);
        assert_eq!(xsf::factorial(n as u16), 1);
        assert_eq!(xsf::factorial(n as u32), 1);
        assert_eq!(xsf::factorial(n as u64), 1);
        assert_eq!(xsf::factorial(n as u128), 1);
    }
}

#[test]
fn test_factorial_max() {
    // anything larger will overflow
    assert_eq!(xsf::factorial(5_u8), FACT1_5);
    assert_eq!(xsf::factorial(8_u16), FACT1_8);
    assert_eq!(xsf::factorial(12_u32), FACT1_12);
    assert_eq!(xsf::factorial(20_u64), FACT1_20);
    assert_eq!(xsf::factorial(34_u128), FACT1_34);
}

#[test]
fn test_factorial_checked_max() {
    assert_eq!(xsf::factorial_checked(5_u8), Some(FACT1_5));
    assert_eq!(xsf::factorial_checked(8_u16), Some(FACT1_8));
    assert_eq!(xsf::factorial_checked(12_u32), Some(FACT1_12));
    assert_eq!(xsf::factorial_checked(20_u64), Some(FACT1_20));
    assert_eq!(xsf::factorial_checked(34_u128), Some(FACT1_34));
}

#[test]
fn test_factorial_checked_overflow() {
    assert_eq!(xsf::factorial_checked(6_u8), None);
    assert_eq!(xsf::factorial_checked(9_u16), None);
    assert_eq!(xsf::factorial_checked(13_u32), None);
    assert_eq!(xsf::factorial_checked(21_u64), None);
    assert_eq!(xsf::factorial_checked(35_u128), None);
}

const FACT2_7: u8 = 105;
const FACT2_12: u16 = 46_080;
const FACT2_20: u32 = 3_715_891_200;
const FACT2_33: u64 = 6_332_659_870_762_850_625;
const FACT2_56: u128 = 81_842_841_814_930_553_085_241_614_925_824_000_000;

#[test]
fn test_multifactorial_k2_one() {
    for n in 0..=2 {
        let y = if n == 0 { 1 } else { n };
        assert_eq!(xsf::multifactorial(n as u8, 2), y as u8, "{n}");
        assert_eq!(xsf::multifactorial(n as u16, 2), y as u16, "{n}");
        assert_eq!(xsf::multifactorial(n as u32, 2), y as u32, "{n}");
        assert_eq!(xsf::multifactorial(n as u64, 2), y as u64, "{n}");
        assert_eq!(xsf::multifactorial(n as u128, 2), y as u128, "{n}");
    }
}

#[test]
fn test_multifactorial_k2_max() {
    assert_eq!(xsf::multifactorial(7_u8, 2), FACT2_7);
    assert_eq!(xsf::multifactorial(12_u16, 2), FACT2_12);
    assert_eq!(xsf::multifactorial(20_u32, 2), FACT2_20);
    assert_eq!(xsf::multifactorial(33_u64, 2), FACT2_33);
    assert_eq!(xsf::multifactorial(56_u128, 2), FACT2_56);
}

#[test]
fn test_multifactorial_k() {
    const N: [u64; 9] = [0, 1, 2, 3, 4, 5, 6, 7, 8];
    assert_eq!(
        N.map(|n| xsf::multifactorial(n, 1)),
        [1, 1, 2, 6, 24, 120, 720, 5040, 40320],
    );
    assert_eq!(
        N.map(|n| xsf::multifactorial(n, 2)),
        [1, 1, 2, 3, 8, 15, 48, 105, 384],
    );
    assert_eq!(
        N.map(|n| xsf::multifactorial(n, 3)),
        [1, 1, 2, 3, 4, 10, 18, 28, 80],
    );
    assert_eq!(
        N.map(|n| xsf::multifactorial(n, 4)),
        [1, 1, 2, 3, 4, 5, 12, 21, 32],
    );
}
