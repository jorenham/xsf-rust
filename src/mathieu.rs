use crate::bindings::xsf_impl;

xsf_impl!(cem_cva, (m: f64, q: f64), "Characteristic value of even Mathieu functions");
xsf_impl!(sem_cva, (m: f64, q: f64), "Characteristic value of odd Mathieu functions");

// TODO: cen, sem, mcm1, msm1, mcm2, msm2
