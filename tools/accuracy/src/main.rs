//! Evaluates xsf-rust or Rust std functions for the accuracy table (see `accuracy.py`).
//!
//! Usage: `xsf-accuracy {xsf|std} < INPUT > OUTPUT`
//!
//! Each input line is `<function> <arg>...`; each output line is the f64 result, printed in
//! the shortest form that round-trips. The function names are those of `specs.py`, and must
//! match the C++ harness `harness.cpp`.

use std::io::{self, BufRead, BufWriter, Write};

use num_complex::Complex;

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn eval_xsf(f: &str, a: &[f64]) -> Option<f64> {
    Some(match (f, a) {
        // gamma and related
        ("gamma", &[x]) => xsf::gamma(x),
        ("gammaln", &[x]) => xsf::gammaln(x),
        ("rgamma", &[x]) => xsf::rgamma(x),
        ("digamma", &[x]) => xsf::digamma(x),
        ("digammainv", &[y]) => xsf::digammainv(y),
        ("beta", &[a, b]) => xsf::beta(a, b),
        ("betaln", &[a, b]) => xsf::betaln(a, b),
        ("binom", &[n, k]) => xsf::binom(n, k),
        ("pow_rising", &[a, m]) => xsf::pow_rising(a, m),
        // error function and normal distribution
        ("erf", &[x]) => xsf::erf(x),
        ("erfc", &[x]) => xsf::erfc(x),
        ("erfcx", &[x]) => xsf::erfcx(x),
        ("erfinv", &[y]) => xsf::erfinv(y),
        ("erfcinv", &[y]) => xsf::erfcinv(y),
        ("ndtr", &[x]) => xsf::ndtr(x),
        ("ndtri", &[p]) => xsf::ndtri(p),
        ("log_ndtr", &[x]) => xsf::log_ndtr(x),
        // elementary functions
        ("expm1", &[x]) => xsf::expm1(x),
        ("log1p", &[x]) => xsf::log1p(x),
        ("cbrt", &[x]) => xsf::cbrt(x),
        ("exp2", &[x]) => xsf::exp2(x),
        ("exp10", &[x]) => xsf::exp10(x),
        ("sinpi", &[x]) => xsf::sinpi(x),
        ("cospi", &[x]) => xsf::cospi(x),
        ("cosine_cdf", &[x]) => xsf::cosine_cdf(x),
        ("cosine_invcdf", &[p]) => xsf::cosine_invcdf(p),
        // Bessel and Struve functions
        ("bessel_j0", &[x]) => xsf::bessel_j0(x),
        ("bessel_j1", &[x]) => xsf::bessel_j1(x),
        ("bessel_y0", &[x]) => xsf::bessel_y0(x),
        ("bessel_y1", &[x]) => xsf::bessel_y1(x),
        ("bessel_i0", &[x]) => xsf::bessel_i0(x),
        ("bessel_i1", &[x]) => xsf::bessel_i1(x),
        ("bessel_k0", &[x]) => xsf::bessel_k0(x),
        ("bessel_k1", &[x]) => xsf::bessel_k1(x),
        ("bessel_j", &[v, x]) => xsf::bessel_j(v, x),
        ("bessel_y", &[v, x]) => xsf::bessel_y(v, x),
        ("bessel_i", &[v, x]) => xsf::bessel_i(v, x),
        ("bessel_k", &[v, x]) => xsf::bessel_k(v, x),
        ("struve_h", &[v, x]) => xsf::struve_h(v, x),
        ("struve_l", &[v, x]) => xsf::struve_l(v, x),
        // zeta functions
        ("zeta", &[x, q]) => xsf::zeta(x, q),
        ("riemann_zeta", &[x]) => xsf::riemann_zeta(x),
        ("zetac", &[x]) => xsf::zetac(x),
        // exponential integrals
        ("expi", &[x]) => xsf::expi(x),
        ("expn", &[n, x]) => xsf::expn(n as u32, x),
        // incomplete gamma and beta functions
        ("gammainc", &[a, x]) => xsf::gammainc(a, x),
        ("gammaincc", &[a, x]) => xsf::gammaincc(a, x),
        ("betainc", &[a, b, x]) => xsf::betainc(a, b, x),
        ("fdtr", &[a, b, x]) => xsf::fdtr(a, b, x),
        ("fdtrc", &[a, b, x]) => xsf::fdtrc(a, b, x),
        // elliptic integrals
        ("ellipk", &[m]) => xsf::ellipk(m),
        ("ellipe", &[m]) => xsf::ellipe(m),
        ("ellipkinc", &[phi, m]) => xsf::ellipkinc(phi, m),
        ("ellipeinc", &[phi, m]) => xsf::ellipeinc(phi, m),
        // other
        ("hyp0f1", &[b, z]) => xsf::hyp0f1(b, z),
        ("lambertw", &[x]) => xsf::lambertw(Complex::new(x, 0.0), 0, 1e-8).re,
        ("spence", &[x]) => xsf::spence(x),
        _ => return None,
    })
}

fn eval_std(f: &str, a: &[f64]) -> Option<f64> {
    Some(match (f, a) {
        ("expm1", &[x]) => x.exp_m1(),
        ("log1p", &[x]) => x.ln_1p(),
        ("cbrt", &[x]) => x.cbrt(),
        ("exp2", &[x]) => x.exp2(),
        _ => return None,
    })
}

fn main() {
    let eval = match std::env::args().nth(1).as_deref() {
        Some("xsf") => eval_xsf,
        Some("std") => eval_std,
        _ => panic!("usage: xsf-accuracy {{xsf|std}} < INPUT > OUTPUT"),
    };
    let mut out = BufWriter::new(io::stdout().lock());
    let mut args = Vec::new();
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let f = it.next().unwrap();
        args.clear();
        args.extend(it.map(|s| s.parse::<f64>().unwrap()));
        let y = eval(f, &args).unwrap_or_else(|| panic!("unknown function: {line}"));
        writeln!(out, "{y:?}").unwrap();
    }
}
