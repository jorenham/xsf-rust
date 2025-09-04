use crate::bindings::xsf_impl;

xsf_impl!(expit, (x: f64), "Expit function, `1/(1 + exp(-x))`");
xsf_impl!(exprel, (x: f64), "Relative error exponential, `(exp(x) - 1)/x`");

xsf_impl!(logit, (x: f64), "Logit function, `log(x / (1 - x))`");
xsf_impl!(log_expit, (x: f64), "Log of the expit function, `log(expit(x))`");
xsf_impl!(log1mexp, (x: f64), "Compute `log(1 - exp(x))`");
