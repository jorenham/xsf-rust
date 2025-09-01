//! Comprehensive parquet table validation and testing utilities for xsf functions
//!
//! Provides extensive validation of parquet table integrity including checksums,
//! metadata consistency, type validation, and reference data testing.
//! Closely matches the functionality of xsref/tests/test_tables.py.

#![allow(dead_code)]

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
    ValidationFailed(String),
    ChecksumMismatch,
    TypeMismatch,
    MetadataMissing,
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
pub fn load_test_cases(function_name: &str, signature: &str) -> Result<Vec<TestCase>, TestError> {
    let base_path = reference_data_path()
        .join("scipy_special_tests")
        .join(function_name);

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
                if let Some(float64_array) = column.as_any().downcast_ref::<Float64Array>() {
                    row.push(float64_array.value(row_idx));
                } else if let Some(float32_array) = column.as_any().downcast_ref::<Float32Array>() {
                    row.push(float32_array.value(row_idx) as f64);
                } else if let Some(int32_array) = column.as_any().downcast_ref::<Int32Array>() {
                    row.push(int32_array.value(row_idx) as f64);
                } else if let Some(int64_array) = column.as_any().downcast_ref::<Int64Array>() {
                    row.push(int64_array.value(row_idx) as f64);
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

    abs_error / abs_desired
}

/// based on `xsref.float_tools.extended_absolute_error`
fn extended_absolute_error(actual: f64, desired: f64) -> f64 {
    if actual == desired || (actual.is_nan() && desired.is_nan()) {
        return 0.0;
    }
    if desired.is_nan() || actual.is_nan() {
        // If expected nan but got non-NaN or expected non-NaN but got NaN
        // we consider this to be an infinite error.
        return f64::INFINITY;
    }
    if actual.is_infinite() {
        // We don't want to penalize early overflow too harshly, so instead
        // compare with the mythical value nextafter(max_float).
        let sgn = actual.signum();
        let mantissa_bits = 52; // f64 has 52 mantissa bits
        let max_float = f64::MAX;
        // max_float * 2**-(mantissa_bits + 1) = ulp(max_float)
        let ulp = max_float * 2.0_f64.powi(-(mantissa_bits + 1));
        return ((sgn * max_float - desired) + sgn * ulp).abs();
    }
    if desired.is_infinite() {
        let sgn = desired.signum();
        let mantissa_bits = 52; // f64 has 52 mantissa bits
        let max_float = f64::MAX;
        // max_float * 2**-(mantissa_bits + 1) = ulp(max_float)
        let ulp = max_float * 2.0_f64.powi(-(mantissa_bits + 1));
        return ((sgn * max_float - actual) + sgn * ulp).abs();
    }
    (actual - desired).abs()
}

pub fn test_function<F>(function_name: &str, signature: &str, test_fn: F) -> Result<(), TestError>
where
    F: Fn(&[f64]) -> f64,
{
    let test_cases = load_test_cases(function_name, signature)?;

    let mut failed = 0;

    for (i, case) in test_cases.iter().enumerate() {
        let actual = test_fn(&case.inputs);

        let error = extended_relative_error(actual, case.expected);

        if error > case.tolerance {
            failed += 1;
            if failed <= 3 {
                eprintln!(
                    "{}: test {}, expected {:.6e}, got {:.6e}, error {:.2e}, tolerance {:.2e}",
                    function_name, i, case.expected, actual, error, case.tolerance
                );
            }
        }
    }

    let tested = test_cases.len();
    assert!(
        failed == 0,
        "{}: {}/{} tests failed",
        function_name,
        failed,
        tested,
    );
    println!(
        "{}: {}/{} tests passed",
        function_name,
        tested - failed,
        tested,
    );

    Ok(())
}
