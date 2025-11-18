use core::any::Any;
use core::fmt::Display;
use std::fs::File;
use std::io::Error as IOError;
use std::path::PathBuf;

use arrow::array::{Array, Float64Array, Int32Array, Int64Array};
use arrow::error::ArrowError;
use num_complex::{Complex, c64};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;

use crate::xsf::fp_error_metrics::{
    ExtendedErrorArg, extended_absolute_error, extended_relative_error,
};

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
        format!("{self:.5e}")
    }
}

impl TestOutputValue for Complex<f64> {
    #[inline(always)]
    fn magnitude(&self) -> f64 {
        self.norm()
    }

    #[inline(always)]
    fn format(&self) -> String {
        format!("{:.5e}{:+.5e}i", self.re, self.im)
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
    fn error(&self, expected: Self) -> f64 {
        // max adjusted relative error
        self.values()
            .iter()
            .zip(expected.values().iter())
            .map(|(&a, &e)| extended_relative_error(a, e).min(extended_absolute_error(a, e)))
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
    pub r#in: Vec<f64>,
    pub out: T,
    pub err: f64,
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

fn map_parquet_batches<T, U, F>(file: File, f: F) -> Vec<T>
where
    U: IntoIterator<Item = T>,
    F: FnMut(arrow::array::RecordBatch) -> U,
{
    ParquetRecordBatchReaderBuilder::try_new(file)
        .unwrap()
        .build()
        .unwrap()
        .map(|batch| batch.unwrap())
        .flat_map(f)
        .collect()
}

fn map_parquet_columns<T, F>(file: File, mut f: F) -> Vec<Vec<T>>
where
    F: FnMut(usize, &dyn Any) -> Option<T>,
{
    map_parquet_batches(file, move |batch| {
        let mut rows = Vec::with_capacity(batch.num_rows());
        for row_idx in 0..batch.num_rows() {
            let mut row = Vec::new();
            for col in batch.columns().iter() {
                if let Some(value) = f(row_idx, col.as_any()) {
                    row.push(value);
                }
            }
            rows.push(row);
        }
        rows
    })
}

fn read_parquet_rows(file: File) -> Vec<Vec<f64>> {
    map_parquet_columns(file, |row_idx, col| {
        col.downcast_ref::<Float64Array>()
            .map(|arr| arr.value(row_idx))
            .or_else(|| {
                col.downcast_ref::<Int32Array>()
                    .map(|arr| arr.value(row_idx) as f64)
            })
            .or_else(|| {
                col.downcast_ref::<Int64Array>()
                    .map(|arr| arr.value(row_idx) as f64)
            })
    })
}

fn read_parquet_column(file: File) -> Vec<f64> {
    map_parquet_batches(file, |batch| {
        batch
            .column(0)
            .as_any()
            .downcast_ref::<Float64Array>()
            .map(|arr| arr.values().to_vec())
            .unwrap_or_else(Vec::new)
    })
}

fn read_parquet_output<T: TestOutput>(file: File) -> Vec<T> {
    T::from_parquet_rows(map_parquet_columns(file, |row_idx, col| {
        col.downcast_ref::<Float64Array>()
            .map(|arr| arr.value(row_idx))
    }))
}

fn load_testcases<T: TestOutput>(name: &str, sig: &str) -> Result<Vec<TestCase<T>>, TestError> {
    let tables_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("xsref/tables/scipy_special_tests")
        .join(name);

    let open = |name: &str| File::open(tables_dir.join(format!("{name}.parquet")));

    let table_in = read_parquet_rows(open(&format!("In_{sig}"))?);
    let table_out = read_parquet_output::<T>(open(&format!("Out_{sig}"))?);
    let table_err = read_parquet_column(open(&format!("Err_{sig}_other"))?);

    let test_cases = (0..table_in.len())
        .map(|i| TestCase {
            r#in: table_in[i].clone(),
            out: table_out[i],
            err: table_err[i],
        })
        .collect::<Vec<_>>();

    Ok(test_cases)
}

#[inline(always)]
pub(crate) fn test<T, F>(name: &str, signature: &str, test_fn: F)
where
    T: TestOutput + std::fmt::Debug,
    F: Fn(&[f64]) -> T,
{
    let mut failed = 0;
    let mut total = 0;
    let mut worst_error = -1.0;
    let mut worst_input = Vec::new();
    let mut worst_actual = String::new();
    let mut worst_desired = String::new();

    let cases = load_testcases::<T>(name, signature).unwrap();
    for case in cases.iter() {
        // skip big inputs, as these tend to have inaccurate xsref values
        if case.r#in.iter().any(|&x| x.abs() >= 1e3) {
            continue;
        }

        if name == "ellipj" && case.r#in[1] < 0.0 {
            // ellipj is only defined for 0 <= m <= 1, see
            // https://github.com/scipy/xsref/issues/11#issuecomment-3545242126
            continue;
        }

        let desired = case.out;
        let desired_magnitude = desired.magnitude();

        if desired_magnitude > 1e100 || (desired_magnitude.is_nan() && cases.len() > 1) {
            // skip NaN's and huge values
            continue;
        }

        let actual = test_fn(&case.r#in);

        // special casing for certain inaccurate xsref tables
        let max_error = if name == "itairy" {
            // https://github.com/scipy/xsref/issues/12
            case.err.max(5e-8)
        } else if name == "ellipj" {
            case.err.max(5e-9)
        } else if name == "airy" || name == "it1i0k0" {
            case.err.max(2e-11)
        } else {
            (case.err * 4.0).max(5e-14)
        };
        let error = actual.error(desired);

        if error > max_error {
            failed += 1;
            if error > worst_error {
                worst_error = error;
                worst_input = case.r#in.clone();
                worst_actual = actual.format();
                worst_desired = desired.format();
            }
        }
        total += 1;
    }

    assert!(
        total > 0,
        "{name}: no test cases run (skipped {})\n{cases:?}",
        cases.len()
    );

    assert!(
        failed == 0,
        "{name}: {failed}/{total} tests failed, worst case for {worst_input:?} \
        (actual={worst_actual}; desired={worst_desired}; error={worst_error:.2e})",
    );
}
