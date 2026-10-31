//! Evaluates xsf-rust or Rust std functions for the accuracy table (see `accuracy.py`).
//!
//! Usage: `xsf-accuracy {xsf|std} < INPUT > OUTPUT`
//!
//! Each input line is `<function> <arg>...`; each output line is the f64 result (or the real
//! and imaginary parts of a complex result), printed in the shortest form that round-trips. The function names are those of `specs.py`, and must
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
        ("lambertw", &[x]) => xsf::lambertw(Complex::new(x, 0.0), 0, 1e-8).re,
        ("spence", &[x]) => xsf::spence(x),
        ("wrightomega", &[x]) => xsf::wrightomega(x),
        ("stirling2", &[n, k]) => xsf::stirling2::<f64>(n as u32, k as u32),
        _ => return eval_xsf_2(f, a),
    })
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]
fn eval_xsf_2(f: &str, a: &[f64]) -> Option<f64> {
    Some(match (f, a) {
        // gamma and related
        ("loggamma", &[x]) => xsf::loggamma(x),
        ("gammasgn", &[x]) => xsf::gammasgn(x),
        ("lgam1p", &[x]) => xsf::lgam1p(x),
        ("polygamma", &[n, x]) => xsf::polygamma(n as u32, x),
        ("multigammaln", &[a, d]) => xsf::multigammaln(a, d as u32),
        ("pow_falling", &[x, m]) => xsf::pow_falling(x, m),
        ("perm", &[n, k]) => xsf::perm(n, k as u32),
        // error function and normal distribution
        ("erfi", &[x]) => xsf::erfi(x),
        ("dawsn", &[x]) => xsf::dawsn(x),
        ("voigt_profile", &[x, s, g]) => xsf::voigt_profile(x, s, g),
        ("ndtri_exp", &[y]) => xsf::ndtri_exp(y),
        ("owens_t", &[h, a]) => xsf::owens_t(h, a),
        ("bivariate_normal_cdf", &[h, k, r]) => xsf::bivariate_normal_cdf(h, k, r),
        ("trivariate_normal_cdf", &[h1, h2, h3, r12, r13, r23]) => {
            xsf::trivariate_normal_cdf(h1, h2, h3, r12, r13, r23, 1e-14)
        },
        ("nrdtrimn", &[p, s, x]) => xsf::nrdtrimn(p, s, x),
        ("nrdtrisd", &[m, p, x]) => xsf::nrdtrisd(m, p, x),
        // elementary functions
        ("cosm1", &[x]) => xsf::cosm1(x),
        ("sindg", &[x]) => xsf::sindg(x),
        ("cosdg", &[x]) => xsf::cosdg(x),
        ("tandg", &[x]) => xsf::tandg(x),
        ("cotdg", &[x]) => xsf::cotdg(x),
        ("radian", &[d, m, s]) => xsf::radian(d, m, s),
        ("sinc", &[x]) => xsf::sinc(x),
        ("diric", &[n, x]) => xsf::diric(n as u32, x),
        ("agm", &[a, b]) => xsf::agm(a, b),
        // exponential and logarithmic functions
        ("expit", &[x]) => xsf::expit(x),
        ("logit", &[p]) => xsf::logit(p),
        ("log_expit", &[x]) => xsf::log_expit(x),
        ("exprel", &[x]) => xsf::exprel(x),
        ("log1mexp", &[x]) => xsf::log1mexp(x),
        ("log1pmx", &[x]) => xsf::log1pmx(x),
        ("xlogy", &[x, y]) => xsf::xlogy(x, y),
        ("xlog1py", &[x, y]) => xsf::xlog1py(x, y),
        ("softplus", &[x]) => xsf::softplus(x),
        ("logaddexp", &[x, y]) => xsf::logaddexp(x, y),
        ("logaddexp2", &[x, y]) => xsf::logaddexp2(x, y),
        // Box-Cox transformations and convex analysis
        ("boxcox", &[x, l]) => xsf::boxcox(x, l),
        ("boxcox1p", &[x, l]) => xsf::boxcox1p(x, l),
        ("inv_boxcox", &[y, l]) => xsf::inv_boxcox(y, l),
        ("inv_boxcox1p", &[y, l]) => xsf::inv_boxcox1p(y, l),
        ("entr", &[x]) => xsf::entr(x),
        ("rel_entr", &[x, y]) => xsf::rel_entr(x, y),
        ("kl_div", &[x, y]) => xsf::kl_div(x, y),
        ("huber", &[d, r]) => xsf::huber(d, r),
        ("pseudo_huber", &[d, r]) => xsf::pseudo_huber(d, r),
        // scaled Bessel functions, derivatives, and ratios
        ("bessel_i0e", &[x]) => xsf::bessel_i0e(x),
        ("bessel_i1e", &[x]) => xsf::bessel_i1e(x),
        ("bessel_k0e", &[x]) => xsf::bessel_k0e(x),
        ("bessel_k1e", &[x]) => xsf::bessel_k1e(x),
        ("bessel_je", &[v, x]) => xsf::bessel_je(v, x),
        ("bessel_ye", &[v, x]) => xsf::bessel_ye(v, x),
        ("bessel_ie", &[v, x]) => xsf::bessel_ie(v, x),
        ("bessel_ke", &[v, x]) => xsf::bessel_ke(v, x),
        ("bessel_j_prime", &[v, x, n]) => xsf::bessel_j_prime(v, x, n as u32),
        ("bessel_y_prime", &[v, x, n]) => xsf::bessel_y_prime(v, x, n as u32),
        ("bessel_i_prime", &[v, x, n]) => xsf::bessel_i_prime(v, x, n as u32),
        ("bessel_k_prime", &[v, x, n]) => xsf::bessel_k_prime(v, x, n as u32),
        ("iv_ratio", &[v, x]) => xsf::iv_ratio(v, x),
        ("iv_ratio_c", &[v, x]) => xsf::iv_ratio_c(v, x),
        ("iv_ratioinv", &[v, r]) => xsf::iv_ratioinv(v, r),
        // spherical Bessel functions
        ("sph_bessel_j", &[n, x]) => xsf::sph_bessel_j(n as i64, x),
        ("sph_bessel_y", &[n, x]) => xsf::sph_bessel_y(n as i64, x),
        ("sph_bessel_i", &[n, x]) => xsf::sph_bessel_i(n as i64, x),
        ("sph_bessel_k", &[n, x]) => xsf::sph_bessel_k(n as i64, x),
        ("sph_bessel_j_prime", &[n, x]) => xsf::sph_bessel_j_prime(n as i64, x),
        ("sph_bessel_y_prime", &[n, x]) => xsf::sph_bessel_y_prime(n as i64, x),
        ("sph_bessel_i_prime", &[n, x]) => xsf::sph_bessel_i_prime(n as i64, x),
        ("sph_bessel_k_prime", &[n, x]) => xsf::sph_bessel_k_prime(n as i64, x),
        // integrals of Bessel and Struve functions
        ("it1j0y0.0", &[x]) => xsf::it1j0y0(x).0,
        ("it1j0y0.1", &[x]) => xsf::it1j0y0(x).1,
        ("it2j0y0.0", &[x]) => xsf::it2j0y0(x).0,
        ("it2j0y0.1", &[x]) => xsf::it2j0y0(x).1,
        ("it1i0k0.0", &[x]) => xsf::it1i0k0(x).0,
        ("it1i0k0.1", &[x]) => xsf::it1i0k0(x).1,
        ("it2i0k0.0", &[x]) => xsf::it2i0k0(x).0,
        ("it2i0k0.1", &[x]) => xsf::it2i0k0(x).1,
        ("besselpoly", &[a, l, v]) => xsf::besselpoly(a, l, v),
        ("itstruve0", &[x]) => xsf::itstruve0(x),
        ("it2struve0", &[x]) => xsf::it2struve0(x),
        ("itmodstruve0", &[x]) => xsf::itmodstruve0(x),
        // Kelvin functions
        ("ber", &[x]) => xsf::ber(x),
        ("bei", &[x]) => xsf::bei(x),
        ("ker", &[x]) => xsf::ker(x),
        ("kei", &[x]) => xsf::kei(x),
        ("berp", &[x]) => xsf::berp(x),
        ("beip", &[x]) => xsf::beip(x),
        ("kerp", &[x]) => xsf::kerp(x),
        ("keip", &[x]) => xsf::keip(x),
        // Airy functions
        ("airy[0]", &[x]) => xsf::airy(x)[0],
        ("airy[1]", &[x]) => xsf::airy(x)[1],
        ("airy[2]", &[x]) => xsf::airy(x)[2],
        ("airy[3]", &[x]) => xsf::airy(x)[3],
        ("airy_scaled[0]", &[x]) => xsf::airy_scaled(x)[0],
        ("airy_scaled[1]", &[x]) => xsf::airy_scaled(x)[1],
        ("airy_scaled[2]", &[x]) => xsf::airy_scaled(x)[2],
        ("airy_scaled[3]", &[x]) => xsf::airy_scaled(x)[3],
        ("airy_integrals[0]", &[x]) => xsf::airy_integrals(x)[0],
        ("airy_integrals[1]", &[x]) => xsf::airy_integrals(x)[1],
        ("airy_integrals[2]", &[x]) => xsf::airy_integrals(x)[2],
        ("airy_integrals[3]", &[x]) => xsf::airy_integrals(x)[3],
        // parabolic cylinder functions
        ("pbdv.0", &[v, x]) => xsf::pbdv(v, x).0,
        ("pbdv.1", &[v, x]) => xsf::pbdv(v, x).1,
        ("pbvv.0", &[v, x]) => xsf::pbvv(v, x).0,
        ("pbvv.1", &[v, x]) => xsf::pbvv(v, x).1,
        ("pbwa.0", &[a, x]) => xsf::pbwa(a, x).0,
        ("pbwa.1", &[a, x]) => xsf::pbwa(a, x).1,
        // zeta functions
        ("gen_harmonic", &[n, a]) => xsf::gen_harmonic(n as u64, a),
        ("normalized_gen_harmonic", &[j, k, n, a]) => {
            xsf::normalized_gen_harmonic(j as u64, k as u64, n as u64, a)
        },
        // exponential, sine, cosine, and Fresnel integrals
        ("exp1", &[x]) => xsf::exp1(x),
        ("scaled_exp1", &[x]) => xsf::scaled_exp1(x),
        ("sici.0", &[x]) => xsf::sici(x).0,
        ("sici.1", &[x]) => xsf::sici(x).1,
        ("shichi.0", &[x]) => xsf::shichi(x).0,
        ("shichi.1", &[x]) => xsf::shichi(x).1,
        ("fresnel.0", &[x]) => xsf::fresnel(x).0,
        ("fresnel.1", &[x]) => xsf::fresnel(x).1,
        // inverse incomplete gamma and beta functions
        ("gammaincinv", &[a, p]) => xsf::gammaincinv(a, p),
        ("gammainccinv", &[a, q]) => xsf::gammainccinv(a, q),
        ("betaincinv", &[a, b, y]) => xsf::betaincinv(a, b, y),
        // elliptic functions and integrals
        ("ellipkm1", &[p]) => xsf::ellipkm1(p),
        ("elliprc", &[x, y]) => xsf::elliprc(x, y),
        ("elliprd", &[x, y, z]) => xsf::elliprd(x, y, z),
        ("elliprf", &[x, y, z]) => xsf::elliprf(x, y, z),
        ("elliprg", &[x, y, z]) => xsf::elliprg(x, y, z),
        ("elliprj", &[x, y, z, p]) => xsf::elliprj(x, y, z, p),
        ("ellipj.0", &[u, m]) => xsf::ellipj(u, m).0,
        ("ellipj.1", &[u, m]) => xsf::ellipj(u, m).1,
        ("ellipj.2", &[u, m]) => xsf::ellipj(u, m).2,
        ("ellipj.3", &[u, m]) => xsf::ellipj(u, m).3,
        // hypergeometric functions
        ("hyp0f0", &[z]) => xsf::hyp0f0(z),
        ("hyp1f0", &[a, z]) => xsf::hyp1f0(a, z),
        ("hyp0f1", &[b, z]) => xsf::hyp0f1(b, z),
        ("hyp1f1", &[a, b, x]) => xsf::hyp1f1(a, b, x),
        ("hypu", &[a, b, x]) => xsf::hypu(a, b, x),
        ("hyp2f1", &[a, b, c, x]) => xsf::hyp2f1(a, b, c, x),
        _ => return eval_xsf_3(f, a),
    })
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]
fn eval_xsf_3(f: &str, a: &[f64]) -> Option<f64> {
    Some(match (f, a) {
        // orthogonal polynomials; the `_l` variants have an integer degree
        ("eval_jacobi", &[n, a, b, x]) => xsf::eval_jacobi(n, a, b, x),
        ("eval_jacobi_l", &[n, a, b, x]) => xsf::eval_jacobi(n as i32, a, b, x),
        ("eval_sh_jacobi", &[n, a, b, x]) => xsf::eval_sh_jacobi(n, a, b, x),
        ("eval_sh_jacobi_l", &[n, a, b, x]) => xsf::eval_sh_jacobi(n as i32, a, b, x),
        ("eval_gengegenbauer", &[n, p, q, x]) => xsf::eval_gengegenbauer(n, p, q, x),
        ("eval_gengegenbauer_l", &[n, p, q, x]) => xsf::eval_gengegenbauer(n as i32, p, q, x),
        ("eval_gegenbauer", &[n, a, x]) => xsf::eval_gegenbauer(n, a, x),
        ("eval_gegenbauer_l", &[n, a, x]) => xsf::eval_gegenbauer(n as i32, a, x),
        ("eval_chebyshev_t", &[n, x]) => xsf::eval_chebyshev_t(n, x),
        ("eval_chebyshev_t_l", &[n, x]) => xsf::eval_chebyshev_t(n as i32, x),
        ("eval_chebyshev_u", &[n, x]) => xsf::eval_chebyshev_u(n, x),
        ("eval_chebyshev_u_l", &[n, x]) => xsf::eval_chebyshev_u(n as i32, x),
        ("eval_chebyshev_c", &[n, x]) => xsf::eval_chebyshev_c(n, x),
        ("eval_chebyshev_c_l", &[n, x]) => xsf::eval_chebyshev_c(n as i32, x),
        ("eval_chebyshev_s", &[n, x]) => xsf::eval_chebyshev_s(n, x),
        ("eval_chebyshev_s_l", &[n, x]) => xsf::eval_chebyshev_s(n as i32, x),
        ("eval_sh_chebyshev_t", &[n, x]) => xsf::eval_sh_chebyshev_t(n, x),
        ("eval_sh_chebyshev_t_l", &[n, x]) => xsf::eval_sh_chebyshev_t(n as i32, x),
        ("eval_sh_chebyshev_u", &[n, x]) => xsf::eval_sh_chebyshev_u(n, x),
        ("eval_sh_chebyshev_u_l", &[n, x]) => xsf::eval_sh_chebyshev_u(n as i32, x),
        ("eval_legendre", &[n, x]) => xsf::eval_legendre(n, x),
        ("eval_legendre_l", &[n, x]) => xsf::eval_legendre(n as i32, x),
        ("eval_sh_legendre", &[n, x]) => xsf::eval_sh_legendre(n, x),
        ("eval_sh_legendre_l", &[n, x]) => xsf::eval_sh_legendre(n as i32, x),
        ("eval_genlaguerre", &[n, a, x]) => xsf::eval_genlaguerre(n, a, x),
        ("eval_genlaguerre_l", &[n, a, x]) => xsf::eval_genlaguerre(n as i32, a, x),
        ("eval_laguerre", &[n, x]) => xsf::eval_laguerre(n, x),
        ("eval_laguerre_l", &[n, x]) => xsf::eval_laguerre(n as i32, x),
        ("eval_hermite_h", &[n, x]) => xsf::eval_hermite_h(n, x),
        ("eval_hermite_h_l", &[n, x]) => xsf::eval_hermite_h(n as u32, x),
        ("eval_hermite_he", &[n, x]) => xsf::eval_hermite_he(n, x),
        ("eval_hermite_he_l", &[n, x]) => xsf::eval_hermite_he(n as u32, x),
        // Legendre functions
        ("legendre_p", &[n, x]) => xsf::legendre_p(n as i32, x),
        ("assoc_legendre_p", &[n, m, x]) => xsf::assoc_legendre_p(n as i32, m as i32, x),
        ("assoc_legendre_p_norm", &[n, m, x]) => xsf::assoc_legendre_p_norm(n as i32, m as i32, x),
        ("sph_legendre_p", &[n, m, t]) => xsf::sph_legendre_p(n as i32, m as i32, t),
        ("pmv", &[m, v, x]) => xsf::pmv(m as i32, v, x),
        // discrete distributions
        ("bdtr", &[k, n, p]) => xsf::bdtr(k, n as i32, p),
        ("bdtrc", &[k, n, p]) => xsf::bdtrc(k, n as i32, p),
        ("bdtri", &[k, n, y]) => xsf::bdtri(k, n as i32, y),
        ("nbdtr", &[k, n, p]) => xsf::nbdtr(k as i32, n as i32, p),
        ("nbdtrc", &[k, n, p]) => xsf::nbdtrc(k as i32, n as i32, p),
        ("nbdtri", &[k, n, y]) => xsf::nbdtri(k as i32, n as i32, y),
        ("pdtr", &[k, m]) => xsf::pdtr(k, m),
        ("pdtrc", &[k, m]) => xsf::pdtrc(k, m),
        ("pdtri", &[k, y]) => xsf::pdtri(k as i32, y),
        // continuous distributions
        ("gdtr", &[a, b, x]) => xsf::gdtr(a, b, x),
        ("gdtrc", &[a, b, x]) => xsf::gdtrc(a, b, x),
        ("gdtria", &[p, b, x]) => xsf::gdtria(p, b, x),
        ("gdtrib", &[a, p, x]) => xsf::gdtrib(a, p, x),
        ("gdtrix", &[a, b, p]) => xsf::gdtrix(a, b, p),
        ("chdtr", &[v, x]) => xsf::chdtr(v, x),
        ("chdtrc", &[v, x]) => xsf::chdtrc(v, x),
        ("chdtri", &[v, y]) => xsf::chdtri(v, y),
        ("fdtri", &[a, b, y]) => xsf::fdtri(a, b, y),
        ("stdtr", &[v, t]) => xsf::stdtr(v, t),
        ("stdtri", &[v, p]) => xsf::stdtri(v, p),
        ("tukeylambdacdf", &[x, l]) => xsf::tukeylambdacdf(x, l),
        ("von_mises_cdf", &[k, x]) => xsf::von_mises_cdf(k, x),
        ("geninvgauss_pdf", &[x, p, b]) => xsf::geninvgauss_pdf(x, p, b),
        ("geninvgauss_logpdf", &[x, p, b]) => xsf::geninvgauss_logpdf(x, p, b),
        ("genhyperbolic_pdf", &[x, p, a, b]) => xsf::genhyperbolic_pdf(x, p, a, b),
        ("genhyperbolic_logpdf", &[x, p, a, b]) => xsf::genhyperbolic_logpdf(x, p, a, b),
        ("cdf_cvm_inf", &[x]) => xsf::cdf_cvm_inf(x),
        // Kolmogorov-Smirnov distributions
        ("smirnov", &[n, x]) => xsf::smirnov(n as i32, x),
        ("smirnovc", &[n, x]) => xsf::smirnovc(n as i32, x),
        ("smirnovp", &[n, x]) => xsf::smirnovp(n as i32, x),
        ("smirnovi", &[n, q]) => xsf::smirnovi(n as i32, q),
        ("smirnovci", &[n, p]) => xsf::smirnovci(n as i32, p),
        ("kolmogorov", &[x]) => xsf::kolmogorov(x),
        ("kolmogc", &[x]) => xsf::kolmogc(x),
        ("kolmogp", &[x]) => xsf::kolmogp(x),
        ("kolmogi", &[q]) => xsf::kolmogi(q),
        ("kolmogci", &[p]) => xsf::kolmogci(p),
        // Wright's generalized Bessel function
        ("wright_bessel", &[a, b, x]) => xsf::wright_bessel(a, b, x),
        ("log_wright_bessel", &[a, b, x]) => xsf::log_wright_bessel(a, b, x),
        _ => return None,
    })
}

/// Complex-valued functions; the complex arguments are passed as (real, imaginary) pairs.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn eval_xsf_complex(f: &str, a: &[f64]) -> Option<Complex<f64>> {
    let c = Complex::new;
    Some(match (f, a) {
        ("gamma", &[x, y]) => xsf::gamma(c(x, y)),
        ("loggamma", &[x, y]) => xsf::loggamma(c(x, y)),
        ("rgamma", &[x, y]) => xsf::rgamma(c(x, y)),
        ("digamma", &[x, y]) => xsf::digamma(c(x, y)),
        ("erf", &[x, y]) => xsf::erf(c(x, y)),
        ("erfc", &[x, y]) => xsf::erfc(c(x, y)),
        ("wofz", &[x, y]) => xsf::wofz(c(x, y)),
        ("exp1", &[x, y]) => xsf::exp1(c(x, y)),
        ("expi", &[x, y]) => xsf::expi(c(x, y)),
        ("lambertw", &[x, y]) => xsf::lambertw(c(x, y), 0, 1e-8),
        ("hyp1f1", &[a, b, x, y]) => xsf::hyp1f1(a, b, c(x, y)),
        ("hyp2f1", &[a, b, cc, x, y]) => xsf::hyp2f1(a, b, cc, c(x, y)),
        ("bessel_j", &[v, x, y]) => xsf::bessel_j(v, c(x, y)),
        ("bessel_y", &[v, x, y]) => xsf::bessel_y(v, c(x, y)),
        ("bessel_i", &[v, x, y]) => xsf::bessel_i(v, c(x, y)),
        ("bessel_k", &[v, x, y]) => xsf::bessel_k(v, c(x, y)),
        ("hankel_1", &[v, x, y]) => xsf::hankel_1(v, c(x, y)),
        ("hankel_2", &[v, x, y]) => xsf::hankel_2(v, c(x, y)),
        ("hankel_1e", &[v, x, y]) => xsf::hankel_1e(v, c(x, y)),
        ("hankel_2e", &[v, x, y]) => xsf::hankel_2e(v, c(x, y)),
        ("hankel_1_prime", &[v, x, y, n]) => xsf::hankel_1_prime(v, c(x, y), n as u32),
        ("hankel_2_prime", &[v, x, y, n]) => xsf::hankel_2_prime(v, c(x, y), n as u32),
        ("modified_fresnel_plus.0", &[x]) => xsf::modified_fresnel_plus(x).0,
        ("modified_fresnel_plus.1", &[x]) => xsf::modified_fresnel_plus(x).1,
        ("modified_fresnel_minus.0", &[x]) => xsf::modified_fresnel_minus(x).0,
        ("modified_fresnel_minus.1", &[x]) => xsf::modified_fresnel_minus(x).1,
        ("sph_harm_y", &[n, m, t, p]) => xsf::sph_harm_y(n as usize, m as isize, t, p),
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

type Eval = fn(&str, &[f64]) -> Option<f64>;
type EvalComplex = fn(&str, &[f64]) -> Option<Complex<f64>>;

fn main() {
    let (eval, complex): (Eval, Option<EvalComplex>) = match std::env::args().nth(1).as_deref() {
        Some("xsf") => (eval_xsf, Some(eval_xsf_complex)),
        Some("std") => (eval_std, None),
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
        if let Some(y) = eval(f, &args) {
            writeln!(out, "{y:?}").unwrap();
        } else if let Some(z) = complex.and_then(|g| g(f, &args)) {
            writeln!(out, "{:?} {:?}", z.re, z.im).unwrap();
        } else {
            panic!("unknown function: {line}");
        }
    }
}
