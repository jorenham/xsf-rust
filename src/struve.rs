use crate::bindings::xsf_impl;

xsf_impl!(itstruve0, (x: f64), "Integral of the Struve function of order 0");
xsf_impl!(it2struve0, (x: f64), "Integral related to the Struve function of order 0");
xsf_impl!(itmodstruve0, (x: f64), "Integral of the modified Struve function of order 0");

xsf_impl!(struve_h, (v: f64, x: f64), "Struve `H` function");
xsf_impl!(struve_l, (v: f64, x: f64), "Struve `L` function");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xsref;

    #[test]
    fn test_itstruve0_f64() {
        xsref::test::<f64, _>("itstruve0", "d-d", |x: &[f64]| itstruve0(x[0]));
    }

    #[test]
    fn test_it2struve0_f64() {
        xsref::test::<f64, _>("it2struve0", "d-d", |x: &[f64]| it2struve0(x[0]));
    }

    #[test]
    fn test_itmodstruve0_f64() {
        xsref::test::<f64, _>("itmodstruve0", "d-d", |x: &[f64]| itmodstruve0(x[0]));
    }

    #[test]
    fn test_struve_h_f64() {
        xsref::test::<f64, _>("struve_h", "d_d-d", |x: &[f64]| struve_h(x[0], x[1]));
    }

    #[test]
    fn test_struve_l_f64() {
        xsref::test::<f64, _>("struve_l", "d_d-d", |x: &[f64]| struve_l(x[0], x[1]));
    }
}
