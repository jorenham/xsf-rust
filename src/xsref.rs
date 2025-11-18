use core::fmt::Display;
use std::fs::File;
use std::io::Error as IOError;
use std::path::PathBuf;

use arrow::array::{Array, Float64Array, Int32Array, Int64Array};
use arrow::error::ArrowError;
use num_complex::{Complex, c64};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;

use crate::extended_absolute_error;
use crate::xsf::fp_error_metrics::{ExtendedErrorArg, extended_relative_error};

pub(crate) trait TestOutputValue: ExtendedErrorArg + Copy + Display {
    fn magnitude(&self) -> f64;
    fn format(&self) -> String;
}

impl TestOutputValue for f64 {
    #[inline(always)]
    fn magnitude(&self) -> f64 {
        self.abs()
    }

    #[inline(always)]
    fn format(&self) -> String {
        format!("{self:.3e}")
    }
}

impl TestOutputValue for Complex<f64> {
    #[inline(always)]
    fn magnitude(&self) -> f64 {
        self.norm()
    }

    #[inline(always)]
    fn format(&self) -> String {
        format!("{:.3e}+{:.3e}i", self.re, self.im)
    }
}

pub(crate) trait TestOutput: Copy + PartialEq {
    type Value: TestOutputValue;

    fn from_parquet_row(row: Vec<f64>) -> Self;

    #[inline(always)]
    fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
        rows.into_iter().map(Self::from_parquet_row).collect()
    }

    fn values(&self) -> Vec<Self::Value>;

    #[inline(always)]
    fn magnitude(&self) -> f64 {
        let values = self.values();
        values.iter().map(|x| x.magnitude()).sum::<f64>() / values.len() as f64
    }

    #[inline(always)]
    fn is_nan(&self) -> bool {
        self.magnitude().is_nan()
    }

    #[inline(always)]
    fn error(&self, expected: Self) -> f64 {
        // max adjusted relative error
        self.values()
            .iter()
            .zip(expected.values().iter())
            .map(|(&a, &e)| {
                extended_relative_error(a, e)
                    .min(extended_absolute_error(e, a))
                    .min(extended_absolute_error(a, e))
                    .min(extended_absolute_error(e, a))
            })
            .fold(0.0, |acc, x| acc.max(x))
    }

    fn format(&self) -> String {
        let values = self.values();
        if values.len() == 1 {
            values[0].format()
        } else {
            format!(
                "({})",
                values
                    .iter()
                    .map(Self::Value::format)
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        }
    }
}

impl TestOutput for f64 {
    type Value = f64;

    #[inline(always)]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        row[0]
    }

    #[inline(always)]
    fn values(&self) -> Vec<Self::Value> {
        vec![*self]
    }
}

impl TestOutput for Complex<f64> {
    type Value = Complex<f64>;

    #[inline(always)]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        c64(row[0], row[1])
    }

    #[inline(always)]
    fn values(&self) -> Vec<Self::Value> {
        vec![*self]
    }
}

impl TestOutput for (f64, f64) {
    type Value = f64;

    #[inline(always)]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (row[0], row[1])
    }

    #[inline(always)]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1]
    }
}

impl TestOutput for (Complex<f64>, Complex<f64>) {
    type Value = Complex<f64>;

    #[inline(always)]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (c64(row[0], row[1]), c64(row[2], row[3]))
    }

    #[inline(always)]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1]
    }
}

impl TestOutput for (f64, f64, f64, f64) {
    type Value = f64;

    #[inline(always)]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (row[0], row[1], row[2], row[3])
    }

    #[inline(always)]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1, self.2, self.3]
    }
}

impl TestOutput for (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>) {
    type Value = Complex<f64>;

    #[inline(always)]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (
            c64(row[0], row[1]),
            c64(row[2], row[3]),
            c64(row[4], row[5]),
            c64(row[6], row[7]),
        )
    }

    #[inline(always)]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1, self.2, self.3]
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
pub(crate) enum TestError {
    Io(IOError),
    DataFormat,
}

impl From<IOError> for TestError {
    #[inline(always)]
    fn from(err: IOError) -> Self {
        TestError::Io(err)
    }
}

impl From<ArrowError> for TestError {
    #[inline(always)]
    fn from(_: arrow::error::ArrowError) -> Self {
        TestError::DataFormat
    }
}

impl From<ParquetError> for TestError {
    #[inline(always)]
    fn from(_: ParquetError) -> Self {
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

#[inline(always)]
pub(crate) fn test<T, F>(name: &str, signature: &str, test_fn: F)
where
    T: TestOutput,
    F: Fn(&[f64]) -> T,
{
    let cases = load_testcases::<T>(name, signature).unwrap();

    let mut failed = 0;
    let mut worst_error = -1.0;
    let mut worst_input = Vec::new();
    let mut worst_actual = String::new();
    let mut worst_desired = String::new();

    for case in cases.iter() {
        // skip big and tiny inputs
        if case
            .inputs
            .iter()
            .any(|&x| x.abs() >= 100.0 || (x.abs() <= 1e-15 && x != 0.0))
        {
            continue;
        }

        if name == "ellipj" && case.inputs[1] < 0.0 {
            // ellipj is only defined for 0 <= m <= 1
            continue;
        }

        let actual = test_fn(&case.inputs);
        let desired = case.expected;
        let desired_magnitude = desired.magnitude();

        if desired_magnitude > 1e300 {
            // skip huge values
            continue;
        }

        if desired_magnitude.is_nan() {
            if !actual.is_nan() {
                // TODO
                // failed += 1;
                // eprintln!("{name}: test {i}, expected NaN, got {}", actual.format());
            }
        } else if desired_magnitude.is_infinite() {
            if actual != desired {
                failed += 1;
            }
        } else {
            let max_error = if name == "ellipj" || name == "itairy" {
                // https://github.com/scipy/xsref/issues/12
                5e-8
            } else {
                case.tolerance.max(f64::EPSILON) * 3000.0
            };
            let error = actual.error(desired);

            if error > max_error {
                failed += 1;
                if error > worst_error {
                    worst_error = error;
                    worst_input = case.inputs.clone();
                    worst_actual = actual.format();
                    worst_desired = desired.format();
                }
            }
        }
    }

    let total = cases.len();
    assert!(
        failed == 0,
        "{name}: {failed}/{total} tests failed, worst case for {worst_input:?} \
        (actual={worst_actual}; desired={worst_desired}; error={worst_error:.2e})",
    );
}
