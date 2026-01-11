/// Parabolic cylinder function $W$
///
/// The function is a particular solution to the differential equation:
///
/// $$
/// y\'\' + \left( {1 \over 4} x^2 - a \right) y = 0
/// $$
///
/// Corresponds to [`scipy.special.pbwa`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbwa.html
///
/// # Arguments
/// - `a` - Real parameter
/// - `x` - Real argument
///
/// # Returns
/// - *w*: Value of the function
/// - *wp*: Value of the derivative in `x`
///
/// # See also
/// - [`pbdv`]
/// - [`pbvv`]
#[must_use]
#[inline]
pub fn pbwa(a: f64, x: f64) -> (f64, f64) {
    let (mut w, mut wd) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::pbwa(a, x, &raw mut w, &raw mut wd);
    }
    (w, wd)
}

/// Parabolic cylinder function $D$
///
/// Returns `(d, dp)` the parabolic cylinder function $Dv(x)$ in `d` and the derivative,
/// $Dv\'(x)$ in `dp`.
///
/// Corresponds to [`scipy.special.pbdv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbdv.html
///
/// # Arguments
/// - `v` - Real parameter
/// - `x` - Real argument
///
/// # Returns
/// - *d*: Value of the function
/// - *dp*: Value of the derivative in `x`
///
/// # See also
/// - [`pbwa`]
/// - [`pbvv`]
#[doc(alias = "pbdv_seq")]
#[must_use]
#[inline]
pub fn pbdv(v: f64, x: f64) -> (f64, f64) {
    let (mut d, mut dp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::pbdv(v, x, &raw mut d, &raw mut dp);
    }
    (d, dp)
}

/// Parabolic cylinder function $V$
///
/// Returns `(v, vp)` the parabolic cylinder function $Vv(x)$ in `v` and the derivative,
/// $Vv\'(x)$ in `vp`.
///
/// Corresponds to [`scipy.special.pbvv`][scipy].
///
/// [scipy]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.pbvv.html
///
/// # Arguments
/// - `v` - Real parameter
/// - `x` - Real argument
///
/// # Returns
/// - *v*: Value of the function
/// - *vp*: Value of the derivative in `x`
///
/// # See also
/// - [`pbwa`]
/// - [`pbdv`]
#[doc(alias = "pbvv_seq")]
#[must_use]
#[inline]
pub fn pbvv(v: f64, x: f64) -> (f64, f64) {
    let (mut vv, mut vp) = (f64::NAN, f64::NAN);
    unsafe {
        crate::ffi::xsf::pbvv(v, x, &raw mut vv, &raw mut vp);
    }
    (vv, vp)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_pbwa() {
        xsref::test("pbwa", "d_d-d_d", |x| crate::pbwa(x[0], x[1]));
    }

    #[test]
    fn test_pbdv() {
        xsref::test("pbdv", "d_d-d_d", |x| crate::pbdv(x[0], x[1]));
    }

    #[test]
    fn test_pbvv() {
        xsref::test("pbvv", "d_d-d_d", |x| crate::pbvv(x[0], x[1]));
    }
}
