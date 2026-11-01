use num_complex::c64;

fn as_tuple<T: Clone>(x: &[T; 4]) -> (T, T, T, T) {
    (x[0].clone(), x[1].clone(), x[2].clone(), x[3].clone())
}

#[test]
fn test_airy_f64() {
    xsref::test("airy", "d-d_d_d_d", |x| as_tuple(&xsf::airy(x[0])));
}

#[test]
fn test_airy_c64() {
    xsref::test("airy", "cd-cd_cd_cd_cd", |x| {
        as_tuple(&xsf::airy(c64(x[0], x[1])))
    });
}

#[test]
fn test_airye_f64() {
    xsref::test("airye", "d-d_d_d_d", |x| as_tuple(&xsf::airy_scaled(x[0])));
}

#[test]
fn test_airye_c64() {
    xsref::test("airye", "cd-cd_cd_cd_cd", |x| {
        as_tuple(&xsf::airy_scaled(c64(x[0], x[1])))
    });
}

#[test]
fn test_itairy() {
    xsref::test("itairy", "d-d_d_d_d", |x| {
        as_tuple(&xsf::airy_integrals(x[0]))
    });
}

#[test]
fn test_ai_zeros_0() {
    assert!(xsf::airy_ai_zeros::<0>().iter().all(|a| a.is_empty()));
}

#[test]
fn test_ai_zeros() {
    // >>> from scipy import special
    // >>> a, ap, ai, aip = special.ai_zeros(3)
    let [a, ap, ai, aip] = xsf::airy_ai_zeros::<3>();
    // >>> a
    // array([-2.33810741, -4.08794944, -5.52055983])
    xsref::np_assert_allclose!(&a, &[-2.338_107_41, -4.087_949_44, -5.520_559_83]);
    // >>> ap
    // array([-1.01879297, -3.24819758, -4.82009921])
    xsref::np_assert_allclose!(&ap, &[-1.018_792_97, -3.248_197_58, -4.820_099_21]);
    // >>> ai
    // array([ 0.53565666, -0.41901548,  0.38040647])
    xsref::np_assert_allclose!(&ai, &[0.535_656_66, -0.419_015_48, 0.380_406_47]);
    // >>> aip
    // array([ 0.70121082, -0.80311137,  0.86520403])
    xsref::np_assert_allclose!(&aip, &[0.701_210_82, -0.803_111_37, 0.865_204_03]);
}

#[test]
fn test_bi_zeros_0() {
    assert!(xsf::airy_bi_zeros::<0>().iter().all(|a| a.is_empty()));
}

#[test]
fn test_bi_zeros() {
    // >>> from scipy import special
    // >>> b, bp, bi, bip = special.bi_zeros(3)
    let [b, bp, bi, bip] = xsf::airy_bi_zeros::<3>();
    // >>> b
    // array([-1.17371322, -3.2710933 , -4.83073784])
    xsref::np_assert_allclose!(&b, &[-1.173_713_22, -3.271_093_3, -4.830_737_84]);
    // >>> bp
    // array([-2.29443968, -4.07315509, -5.51239573])
    xsref::np_assert_allclose!(&bp, &[-2.294_439_68, -4.073_155_09, -5.512_395_73]);
    // >>> bi
    // array([-0.45494438,  0.39652284, -0.36796916])
    xsref::np_assert_allclose!(&bi, &[-0.454_944_38, 0.396_522_84, -0.367_969_16]);
    // >>> bip
    // array([ 0.60195789, -0.76031014,  0.83699101])
    xsref::np_assert_allclose!(&bip, &[0.601_957_89, -0.760_310_14, 0.836_991_01]);
}
