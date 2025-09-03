use crate::bindings::xsf_impl;

xsf_impl!(
    iv_ratio,
    (v: f64, x: f64),
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind"
);
xsf_impl!(
    iv_ratio_c,
    (v: f64, x: f64),
    "Compute `iv(v,x)/iv(v-1,x)` of the modified Bessel function of the first kind"
);
