//! Code that is not translated from SciPy or xsf: functions that SciPy doesn't have, and
//! algorithms that the translations use.

mod hahn;
mod hyp2f1_recur;
mod pfq;
pub(crate) mod quadpack;
pub(crate) mod tridiag;

pub use hahn::hahn_q;
pub use pfq::{hyp1f2, hyp2f2, hyp3f2, hyppfq};

pub(crate) use hyp2f1_recur::*;
