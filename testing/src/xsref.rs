use core::fmt::LowerExp;
use std::fs::File;
use std::io::Error as IOError;
use std::path::PathBuf;

use num_complex::{Complex, c64};
use polars::error::PolarsError;
use polars::io::SerReader;
use polars::io::parquet::read::ParquetReader;
use polars::prelude::{AnyValue, Column, DataFrame};

use crate::{ExtendedErrorArg, extended_absolute_error, extended_relative_error};

pub trait TestOutputValue: ExtendedErrorArg + Copy + LowerExp {
    fn magnitude(&self) -> f64;

    fn format(&self) -> String {
        format!("{self:.5e}")
    }
}

impl TestOutputValue for f32 {
    #[inline]
    fn magnitude(&self) -> f64 {
        self.abs() as f64
    }
}

impl TestOutputValue for f64 {
    #[inline]
    fn magnitude(&self) -> f64 {
        self.abs()
    }
}

impl TestOutputValue for Complex<f64> {
    #[inline]
    fn magnitude(&self) -> f64 {
        self.norm()
    }
}

pub trait TestOutput: Copy + PartialEq {
    type Value: TestOutputValue;

    fn from_parquet_row(row: Vec<f64>) -> Self;

    #[inline]
    fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
        rows.into_iter().map(Self::from_parquet_row).collect()
    }

    fn values(&self) -> Vec<Self::Value>;

    #[inline]
    fn magnitude(&self) -> f64 {
        let values = self.values();
        values.iter().map(|x| x.magnitude()).sum::<f64>() / values.len() as f64
    }

    #[inline]
    fn error(&self, expected: Self) -> f64 {
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

impl TestOutput for f32 {
    type Value = f32;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        row[0] as f32
    }

    #[inline]
    fn values(&self) -> Vec<Self::Value> {
        vec![*self]
    }
}

impl TestOutput for f64 {
    type Value = f64;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        row[0]
    }

    #[inline]
    fn values(&self) -> Vec<Self::Value> {
        vec![*self]
    }
}

impl TestOutput for Complex<f64> {
    type Value = Complex<f64>;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        c64(row[0], row[1])
    }

    #[inline]
    fn values(&self) -> Vec<Self::Value> {
        vec![*self]
    }
}

impl TestOutput for (f32, f32) {
    type Value = f32;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (row[0] as f32, row[1] as f32)
    }

    #[inline]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1]
    }
}

impl TestOutput for (f64, f64) {
    type Value = f64;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (row[0], row[1])
    }

    #[inline]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1]
    }
}

impl TestOutput for (Complex<f64>, Complex<f64>) {
    type Value = Complex<f64>;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (c64(row[0], row[1]), c64(row[2], row[3]))
    }

    #[inline]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1]
    }
}

impl TestOutput for (f64, f64, f64, f64) {
    type Value = f64;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (row[0], row[1], row[2], row[3])
    }

    #[inline]
    fn values(&self) -> Vec<Self::Value> {
        vec![self.0, self.1, self.2, self.3]
    }
}

impl TestOutput for (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>) {
    type Value = Complex<f64>;

    #[inline]
    fn from_parquet_row(row: Vec<f64>) -> Self {
        (
            c64(row[0], row[1]),
            c64(row[2], row[3]),
            c64(row[4], row[5]),
            c64(row[6], row[7]),
        )
    }

    #[inline]
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
pub enum TestError {
    Io(IOError),
    DataFormat,
}

impl From<IOError> for TestError {
    #[inline]
    fn from(err: IOError) -> Self {
        TestError::Io(err)
    }
}

impl From<PolarsError> for TestError {
    #[inline]
    fn from(_: PolarsError) -> Self {
        TestError::DataFormat
    }
}

#[inline]
fn read_parquet_df(file: File) -> Result<DataFrame, TestError> {
    Ok(ParquetReader::new(file).finish()?)
}

#[inline]
fn for_each_column_value<F>(column: &Column, mut f: F) -> Result<(), TestError>
where
    F: FnMut(f64) -> Result<(), TestError>,
{
    let series = column.as_materialized_series();

    macro_rules! consume_chunked {
        ($method:ident, $cast:expr) => {
            if let Ok(ca) = series.$method() {
                for value in ca {
                    let value = value.map($cast).unwrap_or(f64::NAN);
                    f(value)?;
                }
                return Ok(());
            }
        };
    }

    consume_chunked!(f64, |v: f64| v);
    consume_chunked!(i64, |v: i64| v as f64);
    consume_chunked!(i32, |v: i32| v as f64);

    for value in series.iter() {
        let _ = f(match value {
            AnyValue::Null => f64::NAN,
            AnyValue::Float64(v) => v,
            AnyValue::Int64(v) => v as f64,
            AnyValue::Int32(v) => v as f64,
            _ => value.try_extract::<f64>().map_err(TestError::from)?,
        });
    }

    Ok(())
}

#[inline]
fn read_parquet_rows(file: File) -> Result<Vec<Vec<f64>>, TestError> {
    let df = read_parquet_df(file)?;
    let height = df.height();
    let width = df.width();
    let mut rows = vec![Vec::with_capacity(width); height];

    for column in df.get_columns() {
        let mut row_idx = 0;
        for_each_column_value(column, |value| {
            if row_idx >= rows.len() {
                Err(TestError::DataFormat)
            } else {
                rows[row_idx].push(value);
                row_idx += 1;
                Ok(())
            }
        })?;

        if row_idx != height {
            return Err(TestError::DataFormat);
        }
    }

    Ok(rows)
}

#[inline]
fn read_parquet_column(file: File) -> Result<Vec<f64>, TestError> {
    let df = read_parquet_df(file)?;
    let column = df.get_columns().first().ok_or(TestError::DataFormat)?;
    let mut values = Vec::with_capacity(df.height());

    for_each_column_value(column, |value| {
        values.push(value);
        Ok(())
    })?;

    Ok(values)
}

#[inline]
fn load_testcases<T: TestOutput>(name: &str, sig: &str) -> Result<Vec<TestCase<T>>, TestError> {
    let tables_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("xsref/tables/scipy_special_tests")
        .join(name);

    let open = |name: &str| File::open(tables_dir.join(format!("{name}.parquet")));

    let table_in = read_parquet_rows(open(&format!("In_{sig}"))?)?;
    let table_out = T::from_parquet_rows(read_parquet_rows(open(&format!("Out_{sig}"))?)?);
    let table_err = read_parquet_column(open(&format!("Err_{sig}_other"))?)?;

    let test_cases = (0..table_in.len())
        .map(|i| TestCase {
            r#in: table_in[i].clone(),
            out: table_out[i],
            err: table_err[i],
        })
        .collect::<Vec<_>>();

    Ok(test_cases)
}

#[inline]
pub fn test<T, F>(name: &str, signature: &str, test_fn: F)
where
    T: TestOutput + std::fmt::Debug,
    F: Fn(&[f64]) -> T,
{
    let mut failed = 0;
    let mut total = 0;
    let mut worst_err = -1.0;
    let mut worst_in = Vec::new();
    let mut worst_actual = String::new();
    let mut worst_expect = String::new();

    let cases = load_testcases::<T>(name, signature).unwrap();
    for case in cases.iter() {
        if case.r#in.iter().any(|&x| x.abs() >= 1e3) {
            continue;
        }

        if name == "ellipj" && case.r#in[1] < 0.0 {
            continue;
        }

        let desired = case.out;
        let desired_magnitude = desired.magnitude();

        if desired_magnitude > 1e100 || (desired_magnitude.is_nan() && cases.len() > 1) {
            continue;
        }

        let actual = test_fn(&case.r#in);

        let max_error = if name == "itairy" {
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
            if error > worst_err {
                worst_err = error;
                worst_in = case.r#in.clone();
                worst_actual = actual.format();
                worst_expect = desired.format();
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
        "{name}: {failed}/{total} tests failed, worst case for {worst_in:?} \
        (actual={worst_actual}; desired={worst_expect}; error={worst_err:.2e})",
    );
}
