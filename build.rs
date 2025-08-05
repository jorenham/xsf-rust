use std::env;
use std::path::PathBuf;

const WRAPPER_FUNCTIONS: &[(&str, &str, &str)] = &[
    // gamma.h
    ("gamma", "double x", "x"),
    ("gammaln", "double x", "x"),
    ("gammasgn", "double x", "x"),
    ("gammainc", "double a, double x", "a, x"),
    ("gammaincinv", "double a, double p", "a, p"),
    ("gammaincc", "double a, double x", "a, x"),
    ("gammainccinv", "double a, double p", "a, p"),
    ("gamma_ratio", "double a, double b", "a, b"),
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

fn generate_wrapper_header(out_dir: &str) -> String {
    let wrapper_header = format!("{out_dir}/xsf_wrapper.h");

    let mut header_content = String::from(
        "#pragma once\n\n#ifdef __cplusplus\nextern \"C\" {\n#endif\n\n// Generated C function declarations\n",
    );

    for (func_name, params, _) in WRAPPER_FUNCTIONS {
        header_content.push_str(&format!("double xsf_{func_name}({params});\n"));
    }

    header_content.push_str("\n#ifdef __cplusplus\n}\n#endif\n");

    std::fs::write(&wrapper_header, header_content).expect("Failed to write wrapper header");

    wrapper_header
}

fn generate_wrapper_cpp(wrapper_header: &str, out_dir: &str) -> String {
    let wrapper_cpp = format!("{out_dir}/xsf_wrapper_impl.cpp");

    let mut cpp_content =
        format!("#include \"{wrapper_header}\"\n#include \"xsf/gamma.h\"\n\nextern \"C\" {{\n\n");

    for (func_name, params, args) in WRAPPER_FUNCTIONS {
        cpp_content.push_str(&format!(
            "double xsf_{func_name}({params}) {{ return xsf::{func_name}({args}); }}\n"
        ));
    }

    cpp_content.push_str("\n} // extern \"C\"\n");

    std::fs::write(&wrapper_cpp, cpp_content).expect("Failed to write wrapper implementation");

    wrapper_cpp
}

fn build_cpp_library(wrapper_cpp: &str, xsf_path: &str) {
    cc::Build::new()
        .cpp(true)
        .std("c++17")
        .file(wrapper_cpp)
        .include(format!("{xsf_path}/include"))
        .flag("-Wno-unused-parameter")
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
