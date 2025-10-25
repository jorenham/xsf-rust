# Copilot Instructions

## Project Overview

- xsf-rust provides Rust bindings to the https://github.com/scipy/xsf/ C++ library, which
  SciPy uses for (a subset of) the `scipy.special` Python functions
- the structure of the Rust bindings closely follows that of scipy.special, rather than the
  underlying xsf C++ library
- in some cases, xsf-rust deviates from the `scipy.special` API for a better fit with Rust
  conventions (e.g., `const` generic type parameters for return array sizes, or `Result`
  return types for functions that can fail)
- most function names are the same as in `scipy.special`, but some have been renamed to
  for the sake of consistency (e.g., `bessel_jn` instead of `jn`, and `assoc_legendre_q_all`
  instead of `lqmn`)
- functions that support both real and complex arguments in `scipy.special` are implemented using
  sealed traits in xsf-rust, emulating function overloading for e.g. `f64` and
  `num_complex::Complex<f64>` types
- the project is currently in beta stage, and backwards-incompatible changes may still occur

## Documentation

- 100% of the public API is documented using Rustdoc comments
- the rustdoc comments should be consistent with the scipy.special documentation, adapted
  to Rust conventions where necessary, and include a link to the relevant scipy.special
  documentation page
- the `xsf` crate-level documentation mimics the structure of the scipy.special
  documentation, with sections for each group of related functions (e.g., Bessel functions,
  elliptic integrals, etc.), and a brief overview of each function, i.e.
  <https://docs.scipy.org/doc/scipy/reference/special.html>
- include relevant "See also" links to related functions in the rustdoc comments, as
  appropriate, following the scipy.special documentation style
- rustdoc does not support Latex math rendering, so use simple html tags (e.g. `<sub>` and
  `<sup>`), Unicode characters (like π, ∫, and ∛), and italic/bold formatting to represent
  mathematical expressions as closely as possible
- in case of multiple return values, document each one individually in the `# Returns` section
  using bullet points in a single line, including the mathematical expression for each return value
  where applicable

## Testing

- xsf itself uses the https://github.com/scipy/xsref test suite, containing "tables" of
  precomputed values for many special functions
- most xsf-rust functions are tested against these precomputed values, mimicking the testing
  approach used in xsf
- for functions that have no xsref tables available, rust translations of the relevant
  scipy.special Python tests are used instead
- 100% coverage is required
- `src/xsref.rs` contains the testing infrastructure for xsref table tests
- `src/macros.rs` contains internal helper macros `np_assert_array_eq!` and
  `np_assert_equal!`, which mimic the behavior of the similarly named NumPy testing
  functions from `numpy.testing`

## Bindings

- The C++ wrapper functions are generated in `build.rs`
- Bindgen is used to generate the bindings from these generated C++ wrapper functions at build time
- the generated bindings are located in `src/ffi.rs`

## Module structure

- the `numpy` module contains pure-rust translations from `numpy` functions that are wrapped or
  reexported in `scipy.special`
- the `scipy_special` module contains pure-rust implementations of some functions from
  `scipy.special` that are not available in xsf
- the `xsf` module contains the bindings to (or sometimes pure-rust implementations of) the xsf
  C++ functions
- these modules are not publicly exposed, but their contents are reexported at the crate root
  for easier access, i.e. the public API of the `xsf` crate is completely "flat"

## Code style

- Aim for idiomatic Rust code, following the Rust API Guidelines
  (https://rust-lang.github.io/api-guidelines/)
- Prefer rust code over C++ code where possible
- Direct translations of C++ or Python code should include the original code as inline comments
  for reference, and a permalink to the relevant source file
- Everywhere else comments should be avoided, and self-explanatory code should be preferred
  wherever possible
- prefer `core` over `std` where possible
- use `num-complex` crate for complex number support
- prefer arrays and `const` generics over `Vec`s for fixed-size collections where possible
