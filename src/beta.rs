use crate::bindings::xsf_impl;

xsf_impl!(beta, (a: f64, b: f64), "Beta function");
xsf_impl!(betaln, (a: f64, b: f64), "Natural log of `|beta|`");
