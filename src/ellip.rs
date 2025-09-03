use crate::bindings::xsf_impl;

xsf_impl!(ellipk, (m: f64), "Complete elliptic integral of the first kind");
xsf_impl!(ellipkm1, (p: f64), "Complete elliptic integral of the first kind around `m = 1`");
xsf_impl!(ellipkinc, (phi: f64, m: f64), "Incomplete elliptic integral of the first kind");

xsf_impl!(ellipe, (m: f64), "Complete elliptic integral of the second kind");
xsf_impl!(ellipeinc, (phi: f64, m: f64), "Incomplete elliptic integral of the second kind");
