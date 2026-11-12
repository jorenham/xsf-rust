use core::fmt::LowerExp;
use std::fs::File;
use std::io::Error as IOError;
use std::path::PathBuf;

use num_complex::{Complex, c64};
use polars::error::PolarsError;
use polars::io::SerReader;
use polars::io::parquet::read::ParquetReader;
use polars::prelude::{Column, DataFrame, DataType};

use crate::{ExtendedErrorArg, extended_relative_error};

pub trait TestOutputValue: ExtendedErrorArg + Copy + LowerExp {
    fn format(&self) -> String {
        format!("{self:.5e}")
    }
}

impl TestOutputValue for f32 {}
impl TestOutputValue for f64 {}
impl TestOutputValue for Complex<f64> {}

pub trait TestOutput: Copy {
    type Value: TestOutputValue;

    fn from_parquet_row(row: Vec<f64>) -> Self;

    #[inline]
    fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
        rows.into_iter().map(Self::from_parquet_row).collect()
    }

    fn values(&self) -> Vec<Self::Value>;

    /// All outputs NaN
    #[inline]
    fn nan() -> Self {
        Self::from_parquet_row(vec![f64::NAN; 8])
    }

    /// The extended relative error of each output
    #[inline]
    fn errors(&self, expected: Self) -> Vec<f64> {
        self.values()
            .iter()
            .zip(expected.values())
            .map(|(&a, e)| extended_relative_error(a, e))
            .collect()
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
    pub err: Vec<f64>,
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
fn column_values(column: &Column) -> Result<Vec<f64>, TestError> {
    let series = column.as_materialized_series().cast(&DataType::Float64)?;
    Ok(series
        .f64()?
        .iter()
        .map(|x| x.unwrap_or(f64::NAN))
        .collect())
}

#[inline]
fn read_parquet_rows(file: File) -> Result<Vec<Vec<f64>>, TestError> {
    let df = read_parquet_df(file)?;
    let mut rows = vec![Vec::with_capacity(df.width()); df.height()];

    // every column of a `DataFrame` has `height` values
    for column in df.columns() {
        for (row, value) in rows.iter_mut().zip(column_values(column)?) {
            row.push(value);
        }
    }

    Ok(rows)
}

#[inline]
fn load_testcases<T: TestOutput>(name: &str, sig: &str) -> Result<Vec<TestCase<T>>, TestError> {
    let tables_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("xsref/src/xsref/tables/scipy_special_tests")
        .join(name);

    let open = |name: &str| File::open(tables_dir.join(format!("{name}.parquet")));

    let table_in = read_parquet_rows(open(&format!("In_{sig}"))?)?;
    let table_out = T::from_parquet_rows(read_parquet_rows(open(&format!("Out_{sig}"))?)?);
    let table_err = read_parquet_rows(open(&format!("Err_{sig}_other"))?)?;

    let test_cases = (0..table_in.len())
        .map(|i| TestCase {
            r#in: table_in[i].clone(),
            out: table_out[i],
            err: table_err[i].clone(),
        })
        .collect::<Vec<_>>();

    Ok(test_cases)
}

/// `f(n)` if `x` is an integer below 2^63 in magnitude that fits in an `N`, and NaN otherwise,
/// like in xsf for the integer arguments that it takes as `double`
#[inline]
pub fn with_int<N: TryFrom<i64>, T: TestOutput>(x: f64, f: impl FnOnce(N) -> T) -> T {
    if x.fract() == 0.0 && x.abs() < 2f64.powi(63) {
        N::try_from(x as i64).map_or_else(|_| T::nan(), f)
    } else {
        T::nan()
    }
}

/// The smallest relative tolerance of output `i` of a test case, for where the reference or the
/// function is less accurate than the tolerance of the table; infinite if the reference is wrong
fn rtol_min(name: &str, x: &[f64], i: usize) -> f64 {
    match name {
        // https://github.com/scipy/xsf/pull/85#issuecomment-3738637926
        "hyp2f1" if x[0] == -1.0 => f64::INFINITY,
        // xsf computes the wrong imaginary part of F+-(x) for x < 0, which is fixed in xsf-rust
        "modified_fresnel_plus" | "modified_fresnel_minus" if x[0] < 0.0 => f64::INFINITY,
        // the reference for x = 1 is xsf's inaccurate value, which is fixed in xsf-rust
        "kolmogp" if x[0] == 1.0 => f64::INFINITY,
        // the Kelvin derivative references at |x| = 10 are xsf's inaccurate values (scipy/xsf#232)
        "berp" | "beip" | "kerp" | "keip" if x[0].abs() == 10.0 => f64::INFINITY,
        "kelvin" if i >= 2 && x[0].abs() == 10.0 => f64::INFINITY,
        // xsf returns NaN for m = 10000
        "cem" | "sem" if x[0] == 1e4 => f64::INFINITY,
        // xsf returns NaN instead of the limit for m = 1 (scipy/xsf#375)
        "ellipj" if x[0].is_infinite() && x[1] == 1.0 => f64::INFINITY,
        // the phase is reconstructed from sn in the test (scipy/xsref#11)
        "ellipj" if i == 3 => 5e-9,
        // xsf's `fcoef` is inaccurate for even m (scipy/xsf#343)
        "sem" if x[0] % 2.0 == 0.0 => 1e-8,
        // the references for even m are xsf's values, which differ by up to 1e-7 (unverified,
        // presumably also from `fcoef`)
        "msm1" if x[0] % 2.0 == 0.0 => 2e-7,
        "msm2" if x[0] % 2.0 == 0.0 => 2e-13,
        // the references of the derivatives are xsf's inaccurate values
        "pbwa" if i == 1 => 1e-11,
        "pbdv" if i == 1 && x[0] <= -5.0 => 1e-12,
        "pbvv" if i == 1 && x[0] > 0.0 && x[1] == 10.0 => 1e-8,
        // the references of V' are xsf's wrong values (scipy/xsf#374)
        "pbvv" if i == 1 && x[0] > 0.0 && x[1] == -10.0 => f64::INFINITY,
        // ill-conditioned for small a
        "gammaincinv" | "gammainccinv" if (0.0..0.01).contains(&x[0]) => 4e-12,
        // U(a, 1, x) is ill-conditioned near its zeros
        "hyperu" if x[1] == 1.0 && HYPERU_ZEROS.contains(&(x[0], x[2])) => 1e-12,
        // xsf overflows early for a complex z = x[1] + x[2] i
        "cyl_bessel_i" if x.len() == 3 && BESSEL_I_OVERFLOWS.contains(&x[0]) => f64::INFINITY,
        // sin(pi v) of the reflection is ill-conditioned for large |v|
        "cyl_bessel_i" if x[0] < -100.0 && x[0].fract() != 0.0 => 5e-15 * x[0].abs(),
        _ => 0.0,
    }
}

/// The `(a, x)` of the `hyperu` test cases that are nearest to a zero of `U(a, 1, x)`
const HYPERU_ZEROS: [(f64, f64); 5] = [
    (-0.321_643_286_573_146_34, 0.047_094_717_434_869_746),
    (-0.351_703_406_813_627_3, 0.065_331_008_016_032_07),
    (-0.371_743_486_973_947_9, 0.080_360_917_835_671_36),
    (-0.381_763_527_054_108_2, 0.085_170_488_977_955_93),
    (-0.393_787_575_150_300_6, 0.094_589_232_464_929_87),
];

/// The orders of the complex `cyl_bessel_i` test cases where xsf overflows and the reference
/// isn't finite
const BESSEL_I_OVERFLOWS: [f64; 4] = [
    -3_335.494_571_461_091,
    -2_788.048_924_847_581_5,
    -2_668.868_419_500_347_5,
    -1_132.109_773_438_097_1,
];

#[inline]
pub fn test<T, F>(name: &str, signature: &str, test_fn: F)
where
    T: TestOutput + std::fmt::Debug,
    F: Fn(&[f64]) -> T,
{
    let mut failed = 0;
    let mut worst_err = -1.0;
    let mut worst_in = Vec::new();
    let mut worst_actual = String::new();
    let mut worst_expect = String::new();

    let cases = load_testcases::<T>(name, signature).unwrap();
    for case in &cases {
        let actual = test_fn(&case.r#in);
        let errors = actual.errors(case.out);
        assert_eq!(
            errors.len(),
            case.err.len(),
            "{name}: a tolerance per output"
        );

        // like in xsf, each output has its own tolerance
        let error = (errors.iter().zip(&case.err).enumerate())
            .filter(|&(i, (&error, &tol))| {
                let rtol = (4.0 * tol.max(f64::EPSILON)).max(rtol_min(name, &case.r#in, i));
                error.is_nan() || error > rtol
            })
            .map(|(_, (&error, _))| error)
            .reduce(f64::max);

        if let Some(error) = error {
            failed += 1;
            if error > worst_err {
                worst_err = error;
                worst_in = case.r#in.clone();
                worst_actual = actual.format();
                worst_expect = case.out.format();
            }
        }
    }
    let total = cases.len();
    assert!(total > 0, "{name}: no test cases");
    assert!(
        failed == 0,
        "{name}: {failed}/{total} tests failed, worst case for {worst_in:?} \
        (actual={worst_actual}; desired={worst_expect}; error={worst_err:.2e})",
    );
}
