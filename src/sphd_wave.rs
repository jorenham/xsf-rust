use crate::bindings;

/// Characteristic value of prolate spheroidal function
pub fn prolate_segv(m: u64, n: u64, c: f64) -> f64 {
    unsafe { bindings::prolate_segv(m as f64, n as f64, c) }
}

/// Characteristic value of oblate spheroidal function
pub fn oblate_segv(m: u64, n: u64, c: f64) -> f64 {
    unsafe { bindings::oblate_segv(m as f64, n as f64, c) }
}
