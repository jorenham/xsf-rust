//! Tests for xsf special functions
//!
//! Each function gets individual basic and reference data tests that appear
//! as separate tests in `cargo test` output.

//! Comprehensive parquet table validation and testing utilities for xsf functions
//!
//! Provides extensive validation of parquet table integrity including checksums,
//! metadata consistency, type validation, and reference data testing.
//! Closely matches the functionality of xsref/tests/test_tables.py.

use std::env;
use std::fs::File;
use std::path::PathBuf;

use arrow::array::{Array, Float32Array, Float64Array, Int32Array, Int64Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

#[derive(Debug)]
#[allow(dead_code)]
pub enum TestError {
    Io(std::io::Error),
    DataFormat,
}

impl From<std::io::Error> for TestError {
    fn from(err: std::io::Error) -> Self {
        TestError::Io(err)
    }
}

impl From<arrow::error::ArrowError> for TestError {
    fn from(_: arrow::error::ArrowError) -> Self {
        TestError::DataFormat
    }
}

impl From<parquet::errors::ParquetError> for TestError {
    fn from(_: parquet::errors::ParquetError) -> Self {
        TestError::DataFormat
    }
}

/// Test case with inputs, expected output, and tolerance
#[derive(Debug, Clone)]
pub struct TestCase {
    pub inputs: Vec<f64>,
    pub expected: f64,
    pub tolerance: f64,
}

/// Get path to SciPy reference test data
pub fn reference_data_path() -> PathBuf {
    env::var("XSREF_TABLES_PATH")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/xsref/tables"))
}

/// Load test cases for a function
pub fn load_test_cases(name: &str, signature: &str) -> Result<Vec<TestCase>, TestError> {
    let base_path = reference_data_path().join("scipy_special_tests").join(name);

    let input_file = File::open(base_path.join(format!("In_{}.parquet", signature)))?;
    let output_file = File::open(base_path.join(format!("Out_{}.parquet", signature)))?;

    // Try platform-specific tolerance file, fall back to generic
    let platform = if cfg!(all(target_arch = "x86_64", target_os = "linux")) {
        "gcc-linux-x86_64"
    } else if cfg!(all(target_arch = "aarch64", target_os = "macos")) {
        "clang-darwin-aarch64"
    } else {
        "other"
    };

    let tolerance_file =
        File::open(base_path.join(format!("Err_{}_{}.parquet", signature, platform)))
            .or_else(|_| File::open(base_path.join(format!("Err_{}_other.parquet", signature))))?;

    let inputs = read_parquet_rows(input_file)?;
    let outputs = read_parquet_column(output_file)?;
    let tolerances = read_parquet_column(tolerance_file)?;

    let test_cases = inputs
        .into_iter()
        .zip(outputs)
        .zip(tolerances)
        .map(|((inputs, expected), tolerance)| TestCase {
            inputs,
            expected,
            tolerance,
        })
        .collect();

    Ok(test_cases)
}

/// Read all rows from a parquet file (for multi-column input data)
fn read_parquet_rows(file: File) -> Result<Vec<Vec<f64>>, TestError> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let mut rows = Vec::new();
    for batch_result in reader {
        let batch = batch_result?;
        for row_idx in 0..batch.num_rows() {
            let mut row = Vec::new();
            for col_idx in 0..batch.num_columns() {
                let column = batch.column(col_idx);
                // Handle different data types
                if let Some(f64_array) = column.as_any().downcast_ref::<Float64Array>() {
                    row.push(f64_array.value(row_idx));
                } else if let Some(f32_array) = column.as_any().downcast_ref::<Float32Array>() {
                    row.push(f32_array.value(row_idx) as f64);
                } else if let Some(i32_array) = column.as_any().downcast_ref::<Int32Array>() {
                    row.push(i32_array.value(row_idx) as f64);
                } else if let Some(i64_array) = column.as_any().downcast_ref::<Int64Array>() {
                    row.push(i64_array.value(row_idx) as f64);
                }
            }
            rows.push(row);
        }
    }
    Ok(rows)
}

/// Read first column from a parquet file (for single-column output/tolerance data)
fn read_parquet_column(file: File) -> Result<Vec<f64>, TestError> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let mut values = Vec::new();
    for batch_result in reader {
        let batch = batch_result?;
        if let Some(column) = batch.column(0).as_any().downcast_ref::<Float64Array>() {
            for i in 0..column.len() {
                values.push(column.value(i));
            }
        } else if let Some(column) = batch.column(0).as_any().downcast_ref::<Float32Array>() {
            for i in 0..column.len() {
                values.push(column.value(i) as f64);
            }
        }
    }
    Ok(values)
}

/// based on `xsref.float_tools.extended_relative_error`
fn extended_relative_error(actual: f64, desired: f64) -> f64 {
    let abs_error = extended_absolute_error(actual, desired);
    let abs_desired = if desired == 0.0 {
        // If the desired result is 0.0, normalize by smallest subnormal instead
        // of zero. Some answers are still better than others and we want to guard
        f64::MIN_POSITIVE
    } else if desired.is_infinite() {
        f64::MAX
    } else if desired.is_nan() {
        // extended_relative_error(nan, nan) = 0, otherwise
        // extended_relative_error(x0, x1) with one of x0 or x1 NaN is infinity.
        1.0
    } else {
        desired.abs()
    };

    (abs_error / abs_desired).min(abs_error)
}

/// based on `xsref.float_tools.extended_absolute_error`
fn extended_absolute_error(actual: f64, desired: f64) -> f64 {
    if actual == desired || (actual.is_nan() && desired.is_nan()) {
        return 0.0;
    }
    if desired.is_nan() || actual.is_nan() {
        return f64::INFINITY;
    }

    let mantissa_bits = 52; // f64 has 52 mantissa bits
    let max_float = f64::MAX;
    if actual.is_infinite() {
        let sgn = actual.signum();
        let ulp = max_float * 2.0_f64.powi(-(mantissa_bits + 1));
        return ((sgn * max_float - desired) + sgn * ulp).abs();
    }
    if desired.is_infinite() {
        let sgn = desired.signum();
        let ulp = max_float * 2.0_f64.powi(-(mantissa_bits + 1));
        return ((sgn * max_float - actual) + sgn * ulp).abs();
    }
    (actual - desired).abs()
}

pub fn test_function<F>(name: &str, signature: &str, test_fn: F) -> Result<(), TestError>
where
    F: Fn(&[f64]) -> f64,
{
    let mut failed = 0;

    let test_cases = load_test_cases(name, signature)?;
    for (i, case) in test_cases.iter().enumerate() {
        let actual = test_fn(&case.inputs);

        let is_huge = case.expected.abs() > 1e16;
        let max_error_factor = if is_huge { 4096.0 } else { 16.0 };
        let max_error = case.tolerance.max(f64::EPSILON) * max_error_factor;

        let error = extended_relative_error(actual, case.expected);
        if error > max_error && (!is_huge || error > 1e-13) {
            failed += 1;
            eprintln!(
                "{}: test {}, expected {:.6e}, got {:.6e}, error {:.2e}, tolerance {:.2e}",
                name, i, case.expected, actual, error, max_error
            );
        }
    }

    let tested = test_cases.len();
    assert!(failed == 0, "{}: {}/{} tests failed", name, failed, tested,);
    println!("{}: {}/{} tests passed", name, tested - failed, tested,);

    Ok(())
}

/// Helper macro to generate the test body
macro_rules! test_body {
    ($name:ident, $sig:literal, $test_fn:expr) => {
        paste::paste! {
            #[test]
            fn [<test_ $name>]() {
                let test_fn = $test_fn;
                test_function(stringify!($name), $sig, test_fn).unwrap();
            }
        }
    };
}

/// Generate a test function for xsf functions
macro_rules! generate_tests {
    ($name:ident, "d-d") => {
        test_body!($name, "d-d", |xs: &[f64]| xsf::$name(xs[0]));
    };
    ($name:ident, "d_d-d") => {
        test_body!($name, "d_d-d", |xs: &[f64]| xsf::$name(xs[0], xs[1]));
    };
    ($name:ident, "d_d_d-d") => {
        test_body!($name, "d_d_d-d", |xs: &[f64]| xsf::$name(
            xs[0], xs[1], xs[2]
        ));
    };
    ($name:ident, "d_d_d_d-d") => {
        test_body!($name, "d_d_d_d-d", |xs: &[f64]| xsf::$name(
            xs[0], xs[1], xs[2], xs[3]
        ));
    };
    ($name:ident, "p_d-d") => {
        test_body!($name, "p_d-d", |xs: &[f64]| xsf::$name(xs[0] as i32, xs[1]));
    };
    ($name:ident, "p_p_d-d") => {
        test_body!($name, "p_p_d-d", |xs: &[f64]| xsf::$name(
            xs[0] as i32,
            xs[1] as i32,
            xs[2]
        ));
    };
    ($name:ident, "d_p_d-d") => {
        test_body!($name, "d_p_d-d", |xs: &[f64]| xsf::$name(
            xs[0],
            xs[1] as i32,
            xs[2]
        ));
    };
}

// Basic mathematical functions
generate_tests!(gamma, "d-d");
generate_tests!(erf, "d-d");
generate_tests!(erfc, "d-d");
generate_tests!(cbrt, "d-d");
generate_tests!(expm1, "d-d");
generate_tests!(exp2, "d-d");
generate_tests!(exp10, "d-d");

// Gamma-related functions
generate_tests!(digamma, "d-d");
generate_tests!(gammaln, "d-d");
generate_tests!(gammasgn, "d-d");
generate_tests!(gammainc, "d_d-d");
generate_tests!(gammaincinv, "d_d-d");
generate_tests!(gammaincc, "d_d-d");
generate_tests!(gammainccinv, "d_d-d");
generate_tests!(loggamma, "d-d");
generate_tests!(rgamma, "d-d");

// Bessel functions (single parameter)
generate_tests!(cyl_bessel_j0, "d-d");
generate_tests!(cyl_bessel_j1, "d-d");
generate_tests!(cyl_bessel_y0, "d-d");
generate_tests!(cyl_bessel_y1, "d-d");
generate_tests!(cyl_bessel_i0, "d-d");
generate_tests!(cyl_bessel_i0e, "d-d");
generate_tests!(cyl_bessel_i1, "d-d");
generate_tests!(cyl_bessel_i1e, "d-d");
generate_tests!(cyl_bessel_k0, "d-d");
generate_tests!(cyl_bessel_k0e, "d-d");
generate_tests!(cyl_bessel_k1, "d-d");
generate_tests!(cyl_bessel_k1e, "d-d");

// Bessel functions (two parameters)
generate_tests!(cyl_bessel_j, "d_d-d");
generate_tests!(cyl_bessel_je, "d_d-d");
generate_tests!(cyl_bessel_y, "d_d-d");
generate_tests!(cyl_bessel_ye, "d_d-d");
generate_tests!(cyl_bessel_i, "d_d-d");
generate_tests!(cyl_bessel_ie, "d_d-d");
generate_tests!(cyl_bessel_k, "d_d-d");
generate_tests!(cyl_bessel_ke, "d_d-d");
generate_tests!(iv_ratio, "d_d-d");
generate_tests!(iv_ratio_c, "d_d-d");

// Bessel functions (three parameters)
generate_tests!(besselpoly, "d_d_d-d");

// Beta functions
generate_tests!(beta, "d_d-d");
generate_tests!(betaln, "d_d-d");

// Binomial functions
generate_tests!(binom, "d_d-d");

// Error functions
generate_tests!(erfcx, "d-d");
generate_tests!(erfi, "d-d");
generate_tests!(dawsn, "d-d");
generate_tests!(voigt_profile, "d_d_d-d");

// Exponential integral functions
generate_tests!(exp1, "d-d");
generate_tests!(expi, "d-d");
generate_tests!(scaled_exp1, "d-d");

// Hypergeometric functions
generate_tests!(hyp2f1, "d_d_d_d-d");
// generate_tests!(hyp1f1, "d_d_d-d");  // xsref table only exists for complex

// Kelvin functions
generate_tests!(ber, "d-d");
generate_tests!(bei, "d-d");
generate_tests!(ker, "d-d");
generate_tests!(kei, "d-d");
generate_tests!(berp, "d-d");
generate_tests!(beip, "d-d");
generate_tests!(kerp, "d-d");
generate_tests!(keip, "d-d");

// Legendre functions
// generate_tests!(legendre_p, "p_d-d");  // no xsref table
// generate_tests!(sph_legendre_p, "p_p_d-d");  // no xsref table
generate_tests!(pmv, "d_d_d-d");

// Log and exponential functions
generate_tests!(expit, "d-d");
generate_tests!(exprel, "d-d");
generate_tests!(logit, "d-d");
generate_tests!(log_expit, "d-d");
// generate_tests!(log1mexp, "d-d");  // no xsref table
generate_tests!(log1pmx, "d-d");
generate_tests!(xlogy, "d_d-d");
generate_tests!(xlog1py, "d_d-d");

// Mathieu functions
generate_tests!(cem_cva, "d_d-d");
generate_tests!(sem_cva, "d_d-d");

// Spheroidal wave functions
generate_tests!(prolate_segv, "d_d_d-d");
// generate_tests!(oblate_segv, "d_d_d-d"); // no xsref table??

// Statistical functions
generate_tests!(bdtr, "d_p_d-d");
generate_tests!(bdtrc, "d_p_d-d");
generate_tests!(bdtri, "d_p_d-d");
generate_tests!(chdtr, "d_d-d");
generate_tests!(chdtrc, "d_d-d");
generate_tests!(chdtri, "d_d-d");
generate_tests!(fdtr, "d_d_d-d");
generate_tests!(fdtrc, "d_d_d-d");
generate_tests!(fdtri, "d_d_d-d");
generate_tests!(gdtr, "d_d_d-d");
generate_tests!(gdtrc, "d_d_d-d");
generate_tests!(kolmogorov, "d-d");
generate_tests!(kolmogc, "d-d");
generate_tests!(kolmogi, "d-d");
generate_tests!(kolmogp, "d-d");
generate_tests!(ndtr, "d-d");
generate_tests!(ndtri, "d-d");
// generate_tests!(log_ndtr, "d-d");  // no xsref table
generate_tests!(nbdtr, "p_p_d-d");
generate_tests!(nbdtrc, "p_p_d-d");
// generate_tests!(nbdtri, "p_p_d-d");  // no xsref table
generate_tests!(owens_t, "d_d-d");
generate_tests!(pdtr, "d_d-d");
generate_tests!(pdtrc, "d_d-d");
generate_tests!(pdtri, "p_d-d");
generate_tests!(smirnov, "p_d-d");
generate_tests!(smirnovc, "p_d-d");
generate_tests!(smirnovi, "p_d-d");
generate_tests!(smirnovp, "p_d-d");
// generate_tests!(tukeylambdacdf, "d_d-d");  // no xsref table

// Struve functions
generate_tests!(itstruve0, "d-d");
generate_tests!(it2struve0, "d-d");
generate_tests!(itmodstruve0, "d-d");
generate_tests!(struve_h, "d_d-d");
generate_tests!(struve_l, "d_d-d");

// Trigonometric functions
generate_tests!(sinpi, "d-d");
generate_tests!(cospi, "d-d");
generate_tests!(sindg, "d-d");
generate_tests!(cosdg, "d-d");
generate_tests!(tandg, "d-d");
generate_tests!(cotdg, "d-d");
generate_tests!(cosm1, "d-d");
generate_tests!(radian, "d_d_d-d");

// Wright Bessel functions
generate_tests!(wright_bessel, "d_d_d-d");
generate_tests!(log_wright_bessel, "d_d_d-d");

// Zeta functions
generate_tests!(riemann_zeta, "d-d");
generate_tests!(zeta, "d_d-d");
generate_tests!(zetac, "d-d");
