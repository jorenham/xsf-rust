use crate::bindings;

/// Parabolic cylinder function *W*
///
/// The function is a particular solution to the differential equation:
///
/// y″ + (x²/4 - a)y = 0
///
/// # Arguments
///
/// - `a` - Real parameter
/// - `x` - Real argument
///
/// # Returns
///
/// - *w*: Value of the function
/// - *wp*: Value of the derivative in `x`
pub fn pbwa(a: f64, x: f64) -> (f64, f64) {
    let mut w = f64::NAN;
    let mut wd = f64::NAN;
    unsafe {
        bindings::pbwa(a, x, &mut w, &mut wd);
    }
    (w, wd)
}

/// Parabolic cylinder function *D*
///
/// Returns *(d, dp)* the parabolic cylinder function *Dv(x)* in *d* and the derivative,
/// *Dv'(x)* in *dp*.
///
/// # Arguments
///
/// - `v` - Real parameter
/// - `x` - Real argument
///
/// # Returns
///
/// - *d*: Value of the function
/// - *dp*: Value of the derivative in `x`
pub fn pbdv(v: f64, x: f64) -> (f64, f64) {
    let mut d = f64::NAN;
    let mut dp = f64::NAN;
    unsafe {
        bindings::pbdv(v, x, &mut d, &mut dp);
    }
    (d, dp)
}

/// Parabolic cylinder function *V*
///
/// Returns *(v, vp)* the parabolic cylinder function *Vv(x)* in *v* and the derivative,
/// *Vv'(x)* in *vp*.
///
/// # Arguments
///
/// - `v` - Real parameter
/// - `x` - Real argument
///
/// # Returns
///
/// - *v*: Value of the function
/// - *vp*: Value of the derivative in `x`
pub fn pbvv(v: f64, x: f64) -> (f64, f64) {
    let mut vv = f64::NAN;
    let mut vp = f64::NAN;
    unsafe {
        bindings::pbvv(v, x, &mut vv, &mut vp);
    }
    (vv, vp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_pbwa() {
        xsref::test::<(f64, f64), _>("pbwa", "d_d-d_d", |x: &[f64]| pbwa(x[0], x[1]));
    }

    #[test]
    fn test_pbdv() {
        xsref::test::<(f64, f64), _>("pbdv", "d_d-d_d", |x: &[f64]| pbdv(x[0], x[1]));
    }

    #[test]
    fn test_pbvv() {
        xsref::test::<(f64, f64), _>("pbvv", "d_d-d_d", |x: &[f64]| pbvv(x[0], x[1]));
    }
}
