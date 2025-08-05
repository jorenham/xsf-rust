use std::env;
use std::path::PathBuf;

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let xsf_path = format!("{manifest_dir}/xsf");
    let wrapper_path = format!("{manifest_dir}/wrapper.cpp");

    // Build the C++ wrapper
    cc::Build::new()
        .cpp(true)
        .std("c++17")
        .file(&wrapper_path)
        .include(format!("{xsf_path}/include"))
        .flag("-Wno-unused-parameter")
        .compile("xsf_wrapper");

    // Tell cargo to invalidate the built crate if the wrapper changes
    println!("cargo:rerun-if-changed={wrapper_path}");
    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-changed={xsf_path}/include");

    // Generate bindings
    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
