//! Simple parquet test utilities for xsf functions
//!
//! Provides clean, minimal utilities for loading and testing against
//! SciPy reference data when available.

use std::env;
use std::fs::File;
use std::path::PathBuf;

use arrow::array::{Array, Float64Array};
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

/// Get path to SciPy reference test data
pub fn reference_data_path() -> PathBuf {
    env::var("XSREF_TABLES_PATH")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/xsref/tables"))
}

/// Test case with inputs, expected output, and tolerance
#[derive(Debug, Clone)]
pub struct TestCase {
    pub inputs: Vec<f64>,
    pub expected: f64,
    pub tolerance: f64,
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
                if let Some(float_array) = column.as_any().downcast_ref::<Float64Array>() {
                    row.push(float_array.value(row_idx));
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
        }
    }
    Ok(values)
}

/// Test a function against reference data
pub fn test_function<F>(function_name: &str, signature: &str, test_fn: F) -> Result<(), TestError>
where
    F: Fn(&[f64]) -> f64,
{
    let test_cases = load_test_cases(function_name, signature)?;

    let mut failed = 0;

    for (i, case) in test_cases.iter().enumerate() {
        let actual = test_fn(&case.inputs);
        let tolerance = (case.tolerance * 4.0).max(f64::EPSILON);

        if !case.expected.is_finite() {
            continue;
        }

        let error = if case.expected == 0.0 {
            actual.abs()
        } else {
            ((actual - case.expected) / case.expected).abs()
        };

        if error > tolerance {
            failed += 1;
            if failed <= 3 {
                eprintln!(
                    "{}: test {}, expected {:.6e}, got {:.6e}, error {:.2e}",
                    function_name, i, case.expected, actual, error
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
