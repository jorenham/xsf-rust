use std::env;
use std::path::PathBuf;

const XSF_HEADERS: &[&str] = &[
    // "xsf/airy.h",
    "xsf/alg.h",
    // "xsf/bessel.h",
    "xsf/beta.h",
    "xsf/binom.h",
    "xsf/digamma.h",
    // "xsf/ellip.h",
    "xsf/erf.h",
    "xsf/exp.h",
    "xsf/expint.h",
    "xsf/fp_error_metrics.h",
    // "xsf/fresnel.h",
    "xsf/gamma.h",
    "xsf/hyp2f1.h",
    "xsf/iv_ratio.h",
    "xsf/kelvin.h",
    // "xsf/lambertw.h",
    "xsf/legendre.h",
    "xsf/log_exp.h",
    "xsf/log.h",
    "xsf/loggamma.h",
    "xsf/mathieu.h",
    // "xsf/par_cyl.h",
    // "xsf/recur.h",
    // "xsf/sici.h",
    "xsf/specfun.h",
    // "xsf/sph_bessel.h",
    // "xsf/sph_harm.h",
    "xsf/sphd_wave.h",
    "xsf/stats.h",
    "xsf/struve.h",
    "xsf/trig.h",
    "xsf/wright_bessel.h",
    "xsf/zeta.h",
];

const WRAPPER_FUNCTIONS: &[(&str, &str)] = &[
    // alg.h
    ("cbrt", "double x"),
    // beta.h
    ("beta", "double a, double b"),
    ("betaln", "double a, double b"),
    // binom.h
    ("binom", "double n, double k"),
    // digamma.h
    ("digamma", "double z"),
    // erf.h
    ("erf", "double x"),
    ("erfc", "double x"),
    ("erfcx", "double x"),
    ("erfi", "double x"),
    ("voigt_profile", "double x, double sigma, double gamma"),
    ("dawsn", "double x"),
    // exp.h
    ("expm1", "double x"),
    ("exp2", "double x"),
    ("exp10", "double x"),
    // expint.h
    ("exp1", "double x"),
    ("expi", "double x"),
    ("scaled_exp1", "double x"),
    // fp_error_metrics.h
    ("extended_absolute_error", "double actual, double desired"),
    ("extended_relative_error", "double actual, double desired"),
    // gamma.h
    ("gamma", "double x"),
    ("gammaln", "double x"),
    ("gammasgn", "double x"),
    ("gammainc", "double a, double x"),
    ("gammaincinv", "double a, double p"),
    ("gammaincc", "double a, double x"),
    ("gammainccinv", "double a, double p"),
    ("gamma_ratio", "double a, double b"),
    // hyp2f1.h
    ("hyp2f1", "double a, double b, double c, double x"),
    // iv_ratio.h
    ("iv_ratio", "double v, double x"),
    ("iv_ratio_c", "double v, double x"),
    // kelvin.h (TODO: `kelvin`, `klvnzo`)
    ("ber", "double x"),
    ("bei", "double x"),
    ("ker", "double x"),
    ("kei", "double x"),
    ("berp", "double x"),
    ("beip", "double x"),
    ("kerp", "double x"),
    ("keip", "double x"),
    // legendre.h  (TODO: `assoc_legendre_p`, `lqn`, `lqmn`)
    ("legendre_p", "int n, double z"),
    ("sph_legendre_p", "int n, int m, double theta"),
    // log_exp.h
    ("expit", "double x"),
    ("exprel", "double x"),
    ("logit", "double x"),
    ("log_expit", "double x"),
    ("log1mexp", "double x"),
    // log.h
    ("log1p", "double x"),
    ("log1pmx", "double x"),
    ("xlogy", "double x, double y"),
    ("xlog1py", "double x, double y"),
    // loggamma.h
    ("loggamma", "double x"),
    ("rgamma", "double z"),
    // mathieu.h  (TODO: `cen`, `sem`, `mcm1`, `msm1`, `mcm2`, `msm2`)
    ("cem_cva", "double m, double q"),
    ("sem_cva", "double m, double q"),
    // specfun.h  (TODO: `chyp2f1`, `cerf`)
    ("hypu", "double a, double b, double x"),
    ("hyp1f1", "double a, double b, double x"),
    ("pmv", "double m, double v, double x"),
    // sphd_wave.h
    ("prolate_segv", "double m, double n, double c"),
    ("oblate_segv", "double m, double n, double c"),
    // stats.h
    ("bdtr", "double k, int n, double p"),
    ("bdtrc", "double k, int n, double p"),
    ("bdtri", "double k, int n, double y"),
    ("chdtr", "double df, double x"),
    ("chdtrc", "double df, double x"),
    ("chdtri", "double df, double y"),
    ("fdtr", "double a, double b, double x"),
    ("fdtrc", "double a, double b, double x"),
    ("fdtri", "double a, double b, double y"),
    ("gdtr", "double a, double b, double x"),
    ("gdtrc", "double a, double b, double x"),
    ("kolmogorov", "double x"),
    ("kolmogc", "double x"),
    ("kolmogi", "double x"),
    ("kolmogp", "double x"),
    ("ndtr", "double x"),
    ("ndtri", "double x"),
    ("log_ndtr", "double x"),
    ("nbdtr", "int k, int n, double p"),
    ("nbdtrc", "int k, int n, double p"),
    ("nbdtri", "int k, int n, double p"),
    ("owens_t", "double h, double a"),
    ("pdtr", "double k, double m"),
    ("pdtrc", "double k, double m"),
    ("pdtri", "int k, double y"),
    ("smirnov", "int n, double x"),
    ("smirnovc", "int n, double x"),
    ("smirnovi", "int n, double x"),
    ("smirnovp", "int n, double x"),
    ("tukeylambdacdf", "double x, double lmbda"),
    // struve.h
    ("itstruve0", "double x"),
    ("it2struve0", "double x"),
    ("itmodstruve0", "double x"),
    ("struve_h", "double v, double z"),
    ("struve_l", "double v, double z"),
    // trig.h
    ("sinpi", "double x"),
    ("cospi", "double x"),
    ("sindg", "double x"),
    ("cosdg", "double x"),
    ("tandg", "double x"),
    ("cotdg", "double x"),
    ("radian", "double d, double m, double s"),
    ("cosm1", "double x"),
    // wright_bessel.h
    // ("wright_bessel_t", "double a, double b, double x"),
    ("wright_bessel", "double a, double b, double x"),
    ("log_wright_bessel", "double a, double b, double x"),
    // zeta.h
    ("riemann_zeta", "double x"),
    ("zeta", "double x, double q"),
    ("zetac", "double x"),
];

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let xsf_path = format!("{manifest_dir}/xsf");
    let out_dir = env::var("OUT_DIR").unwrap();

    setup_build_dependencies(&xsf_path);

    let wrapper_header = generate_wrapper_header(&out_dir);
    let wrapper_cpp = generate_wrapper_cpp(&wrapper_header, &out_dir);
    build_cpp_library(&wrapper_cpp, &xsf_path);
    generate_bindings(&wrapper_header);
}

fn setup_build_dependencies(xsf_path: &str) {
    println!("cargo:rerun-if-changed={xsf_path}/include");
}

fn extract_param_names(params: &str) -> String {
    params
        .split(',')
        .map(|param| param.split_whitespace().last().unwrap_or(""))
        .collect::<Vec<_>>()
        .join(", ")
}

fn generate_wrapper_header(out_dir: &str) -> String {
    let wrapper_header = format!("{out_dir}/xsf_wrapper.h");

    let mut header_content =
        String::from("#pragma once\n\n#ifdef __cplusplus\nextern \"C\" {\n#endif\n\n");

    for (func_name, params) in WRAPPER_FUNCTIONS {
        header_content.push_str(&format!("double xsf_{func_name}({params});\n"));
    }

    header_content.push_str("\n#ifdef __cplusplus\n}\n#endif\n");

    std::fs::write(&wrapper_header, header_content).expect("Failed to write wrapper header");

    wrapper_header
}

fn generate_wrapper_cpp(wrapper_header: &str, out_dir: &str) -> String {
    let wrapper_cpp = format!("{out_dir}/xsf_wrapper_impl.cpp");

    let mut cpp_content = format!("#include \"{wrapper_header}\"\n");

    // Include all specified XSF headers
    for header in XSF_HEADERS {
        cpp_content.push_str(&format!("#include \"{header}\"\n"));
    }

    cpp_content.push_str("\nextern \"C\" {\n\n");

    for (func_name, params) in WRAPPER_FUNCTIONS {
        let args = extract_param_names(params);
        cpp_content.push_str(&format!(
            "double xsf_{func_name}({params}) {{ return xsf::{func_name}({args}); }}\n"
        ));
    }

    cpp_content.push_str("\n}\n");

    std::fs::write(&wrapper_cpp, cpp_content).expect("Failed to write wrapper implementation");

    wrapper_cpp
}

fn build_cpp_library(wrapper_cpp: &str, xsf_path: &str) {
    cc::Build::new()
        .cpp(true)
        .std("c++17")
        .file(wrapper_cpp)
        .include(format!("{xsf_path}/include"))
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-logical-op-parentheses")
        .compile("xsf_wrapper_impl");
}

fn generate_bindings(wrapper_header: &str) {
    let bindings = bindgen::Builder::default()
        .header(wrapper_header)
        .allowlist_function("xsf_.*")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
