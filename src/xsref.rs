use core::fmt::Display;
use std::fs::File;
use std::io::Error as IOError;
use std::path::PathBuf;

use num_complex::{Complex, c64};
use polars::error::PolarsError;
use polars::io::SerReader;
use polars::io::parquet::read::ParquetReader;
use polars::prelude::{AnyValue, DataFrame};

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

impl From<PolarsError> for TestError {
    #[inline(always)]
    fn from(_: PolarsError) -> Self {
        TestError::DataFormat
    }
}
fn read_parquet_df(file: File) -> Result<DataFrame, TestError> {
    Ok(ParquetReader::new(file).finish()?)
}

fn any_value_to_f64(value: AnyValue<'_>) -> Result<f64, TestError> {
    if matches!(value, AnyValue::Null) {
        return Ok(f64::NAN);
    }

    value.try_extract::<f64>().map_err(TestError::from)
}

fn read_parquet_rows(file: File) -> Result<Vec<Vec<f64>>, TestError> {
    let df = read_parquet_df(file)?;
    let height = df.height();
    let width = df.width();
    let mut rows = vec![Vec::with_capacity(width); height];

    for column in df.get_columns() {
        let series = column.as_materialized_series();
        for (row_idx, value) in series.iter().enumerate() {
            rows[row_idx].push(any_value_to_f64(value)?);
        }
    }

    Ok(rows)
}

fn read_parquet_column(file: File) -> Result<Vec<f64>, TestError> {
    let df = read_parquet_df(file)?;
    let column = df.get_columns().first().ok_or(TestError::DataFormat)?;
    let series = column.as_materialized_series();

    series.iter().map(any_value_to_f64).collect()
}

fn read_parquet_output<T: TestOutput>(file: File) -> Result<Vec<T>, TestError> {
    Ok(T::from_parquet_rows(read_parquet_rows(file)?))
}

fn load_testcases<T: TestOutput>(name: &str, sig: &str) -> Result<Vec<TestCase<T>>, TestError> {
    let tables_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("xsref/tables/scipy_special_tests")
        .join(name);

    let open = |name: &str| File::open(tables_dir.join(format!("{name}.parquet")));

    let table_in = read_parquet_rows(open(&format!("In_{sig}"))?)?;
    let table_out = read_parquet_output::<T>(open(&format!("Out_{sig}"))?)?;
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
