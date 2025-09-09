use core::fmt::Display;
use std::fs::File;
use std::io::Error as IOError;
use std::path::PathBuf;

use arrow::array::{Array, Float64Array, Int32Array, Int64Array};
use arrow::error::ArrowError;
use num_complex::{Complex, c64};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;

use crate::fp_error_metrics::{ExtendedErrorArg, extended_relative_error};

pub(crate) trait TestOutputValue: Copy + ExtendedErrorArg + Display {
    fn magnitude(self) -> f64;
}

impl TestOutputValue for f64 {
    fn magnitude(self) -> f64 {
        self.abs()
    }
}

impl TestOutputValue for Complex<f64> {
    fn magnitude(self) -> f64 {
        self.norm()
    }
}

pub trait TestOutput: Copy {
    type Value: TestOutputValue;

    fn values(&self) -> Vec<Self::Value>;

    fn from_parquet_row(row: Vec<f64>) -> Self;

    fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
        rows.into_iter().map(Self::from_parquet_row).collect()
    }

    fn magnitude(self) -> f64 {
        let values = self.values();
        values.iter().map(|x| x.magnitude()).sum::<f64>() / values.len() as f64
    }

    fn is_huge(self) -> bool {
        self.magnitude() > f64::EPSILON
    }

    fn error(self, expected: Self) -> f64 {
        // max adjusted relative error
        self.values()
            .iter()
            .zip(expected.values().iter())
            .map(|(&a, &e)| extended_relative_error(a, e))
            .fold(0.0, |acc, x| acc.max(x))
    }

    fn format(self) -> String {
        let values = self.values();
        if values.len() == 1 {
            format!("{}", values[0])
        } else {
            format!(
                "({})",
                values
                    .iter()
                    .map(|x| format!("{x}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        }
    }
}

impl TestOutput for f64 {
    type Value = f64;

    fn values(&self) -> Vec<Self::Value> {
        vec![*self]
    }

    fn from_parquet_row(row: Vec<f64>) -> Self {
        row[0]
    }
}

impl TestOutput for Complex<f64> {
    type Value = Complex<f64>;

    fn values(&self) -> Vec<Self::Value> {
        vec![*self]
    }

    fn from_parquet_row(row: Vec<f64>) -> Self {
        c64(row[0], row[1])
    }
}

impl TestOutput for (f64, f64) {
    type Value = f64;

    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1]
    }

    fn from_parquet_row(row: Vec<f64>) -> Self {
        (row[0], row[1])
    }
}

impl TestOutput for (Complex<f64>, Complex<f64>) {
    type Value = Complex<f64>;

    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1]
    }

    fn from_parquet_row(row: Vec<f64>) -> Self {
        (c64(row[0], row[1]), c64(row[2], row[3]))
    }
}

impl TestOutput for (f64, f64, f64, f64) {
    type Value = f64;

    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1, self.2, self.3]
    }

    fn from_parquet_row(row: Vec<f64>) -> Self {
        (row[0], row[1], row[2], row[3])
    }
}

impl TestOutput for (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>) {
    type Value = Complex<f64>;

    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1, self.2, self.3]
    }

    fn from_parquet_row(row: Vec<f64>) -> Self {
        (
            c64(row[0], row[1]),
            c64(row[2], row[3]),
            c64(row[4], row[5]),
            c64(row[6], row[7]),
        )
    }
}

#[derive(Debug)]
struct TestCase<T: TestOutput> {
    pub inputs: Vec<f64>,
    pub expected: T,
    pub tolerance: f64,
}

#[derive(Debug)]
#[allow(dead_code)]
pub enum TestError {
    Io(IOError),
    DataFormat,
}

impl From<IOError> for TestError {
    fn from(err: IOError) -> Self {
        TestError::Io(err)
    }
}

impl From<ArrowError> for TestError {
    fn from(_: arrow::error::ArrowError) -> Self {
        TestError::DataFormat
    }
}

impl From<ParquetError> for TestError {
    fn from(_: parquet::errors::ParquetError) -> Self {
        TestError::DataFormat
    }
}

fn read_parquet_rows(file: File) -> Vec<Vec<f64>> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    let reader = builder.build().unwrap();

    let mut rows = Vec::new();
    for batch_result in reader {
        let batch = batch_result.unwrap();
        for row_idx in 0..batch.num_rows() {
            let mut row = Vec::new();
            for col_idx in 0..batch.num_columns() {
                let column = batch.column(col_idx).as_any();
                if let Some(f64_array) = column.downcast_ref::<Float64Array>() {
                    row.push(f64_array.value(row_idx));
                } else if let Some(i32_array) = column.downcast_ref::<Int32Array>() {
                    row.push(i32_array.value(row_idx) as f64);
                } else if let Some(i64_array) = column.downcast_ref::<Int64Array>() {
                    row.push(i64_array.value(row_idx) as f64);
                }
            }
            rows.push(row);
        }
    }
    rows
}

fn read_parquet_column(file: File) -> Vec<f64> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    let reader = builder.build().unwrap();

    let mut values = Vec::new();
    for batch in reader {
        if let Some(column) = batch
            .unwrap()
            .column(0)
            .as_any()
            .downcast_ref::<Float64Array>()
        {
            values.extend(column.iter().map(|v| v.unwrap()));
        }
    }
    values
}

fn read_parquet_output<T: TestOutput>(file: File) -> Vec<T> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    let reader = builder.build().unwrap();

    let mut rows = Vec::new();
    for batch_result in reader {
        let batch = batch_result.unwrap();
        for row_idx in 0..batch.num_rows() {
            let mut row = Vec::new();
            for col_idx in 0..batch.num_columns() {
                if let Some(f64_array) = batch
                    .column(col_idx)
                    .as_any()
                    .downcast_ref::<Float64Array>()
                {
                    row.push(f64_array.value(row_idx));
                }
            }
            rows.push(row);
        }
    }
    T::from_parquet_rows(rows)
}

fn load_testcases<T: TestOutput>(
    name: &str,
    signature: &str,
) -> Result<Vec<TestCase<T>>, TestError> {
    let tables_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("xsref/tables/scipy_special_tests")
        .join(name);

    let get_table = |name: &str| File::open(tables_dir.join(format!("{name}.parquet")));

    let table_in = read_parquet_rows(get_table(&format!("In_{signature}"))?);
    let table_out = read_parquet_output::<T>(get_table(&format!("Out_{signature}"))?);

    let platform = if cfg!(all(target_arch = "x86_64", target_os = "linux")) {
        "gcc-linux-x86_64"
    } else if cfg!(all(target_arch = "aarch64", target_os = "macos")) {
        "clang-darwin-aarch64"
    } else {
        "other"
    };

    let table_err_file = get_table(&format!("Err_{signature}_{platform}"))
        .or_else(|_| get_table(&format!("Err_{signature}_other")))?;
    let table_err = read_parquet_column(table_err_file);

    let test_cases = table_in
        .into_iter()
        .zip(table_out)
        .zip(table_err)
        .map(|((inputs, expected), tolerance)| TestCase {
            inputs,
            expected,
            tolerance,
        })
        .collect();

    Ok(test_cases)
}

pub fn test<T, F>(name: &str, signature: &str, test_fn: F)
where
    T: TestOutput,
    F: Fn(&[f64]) -> T,
{
    let cases = load_testcases::<T>(name, signature).unwrap();

    let mut failed = 0;
    for (i, case) in cases.iter().enumerate() {
        let actual = test_fn(&case.inputs);
        let desired = case.expected;
        let max_error = case.tolerance.max(f64::EPSILON) * 16.0;

        if max_error.is_huge() || desired.is_huge() && actual.is_huge() {
            continue;
        }

        let error = T::error(actual, desired);
        if error > max_error {
            failed += 1;
            eprintln!(
                "{}: test {}, expected {}, got {}, error {:.2e}, max {:.2e}",
                name,
                i,
                desired.format(),
                actual.format(),
                error,
                max_error
            );
        }
    }

    assert!(
        failed == 0,
        "{}: {}/{} tests failed",
        name,
        failed,
        cases.len()
    );
}
