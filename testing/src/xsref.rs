use core::fmt::LowerExp;
use std::fs;
use std::io::Error as IOError;
use std::path::PathBuf;

use num_complex::{Complex, c64};

use crate::{ExtendedErrorArg, extended_absolute_error, extended_relative_error};

pub trait TestOutputValue: ExtendedErrorArg + Copy + LowerExp {
    /// The real and imaginary parts
    fn parts(self) -> [f64; 2];

    fn format(&self) -> String {
        format!("{self:.5e}")
    }
}

impl TestOutputValue for f32 {
    #[inline]
    fn parts(self) -> [f64; 2] {
        [f64::from(self), 0.0]
    }
}

impl TestOutputValue for f64 {
    #[inline]
    fn parts(self) -> [f64; 2] {
        [self, 0.0]
    }
}

impl TestOutputValue for Complex<f64> {
    #[inline]
    fn parts(self) -> [f64; 2] {
        [self.re, self.im]
    }
}

pub trait TestOutput: Copy {
    type Value: TestOutputValue;

    fn from_row(row: &[f64]) -> Self;

    fn values(&self) -> Vec<Self::Value>;

    /// All outputs NaN
    #[inline]
    fn nan() -> Self {
        Self::from_row(&[f64::NAN; 8])
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
    fn from_row(row: &[f64]) -> Self {
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
    fn from_row(row: &[f64]) -> Self {
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
    fn from_row(row: &[f64]) -> Self {
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
    fn from_row(row: &[f64]) -> Self {
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
    fn from_row(row: &[f64]) -> Self {
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
    fn from_row(row: &[f64]) -> Self {
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
    fn from_row(row: &[f64]) -> Self {
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

/// How the reference of an output was derived
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Status {
    /// xsref's reference, with xsref's tolerance
    Legacy,
    /// Computed independently of xsf
    Computed,
    /// Fixed by a convention, e.g. NaN outside of the domain
    Convention,
}

impl TryFrom<u8> for Status {
    type Error = TestError;

    #[inline]
    fn try_from(status: u8) -> Result<Self, TestError> {
        match status {
            0 => Ok(Self::Legacy),
            1 => Ok(Self::Computed),
            2 => Ok(Self::Convention),
            _ => Err(TestError::DataFormat),
        }
    }
}

/// How an output is compared with its reference
#[derive(Clone, Copy, Debug)]
struct Check {
    tol: f64,
    /// The scale of the absolute error, or 0 for none
    scale: f64,
    status: Status,
}

impl Check {
    /// A legacy output has no scale; any other one has a finite tolerance, and a scale that is
    /// 0 or positive and finite
    fn is_valid(&self) -> bool {
        let finite = 0.0..f64::INFINITY;
        match self.status {
            Status::Legacy => self.scale == 0.0,
            _ => finite.contains(&self.tol) && finite.contains(&self.scale),
        }
    }

    /// The error of an output if it fails, for the smallest relative tolerance `rtol_min`
    fn failure<T: TestOutputValue>(&self, actual: T, expected: T, rtol_min: f64) -> Option<f64> {
        // the error can't tell +inf from `f64::MAX`, or -0 from +0
        let zeros = self.status == Status::Convention;
        if self.status != Status::Legacy && !same_class(actual, expected, zeros) {
            return Some(f64::INFINITY);
        }
        let error = table_error(actual, expected, self.scale);
        let rtol = (4.0 * self.tol.max(f64::EPSILON)).max(rtol_min);
        (error.is_nan() || error > rtol).then_some(error)
    }
}

/// The table metric: the extended relative error, or the extended absolute error over the
/// `scale` (if it isn't 0) where that is smaller
fn table_error<T: ExtendedErrorArg>(actual: T, expected: T, scale: f64) -> f64 {
    let error = extended_relative_error(actual, expected);
    if scale == 0.0 {
        error
    } else {
        error.min(extended_absolute_error(actual, expected) / scale)
    }
}

/// Whether each part of `actual` is in the class of that of `expected`: NaN, an infinity of the
/// same sign, or finite. With `zeros`, a zero of the same sign is a class of its own.
fn same_class<T: TestOutputValue>(actual: T, expected: T, zeros: bool) -> bool {
    (actual.parts().into_iter().zip(expected.parts())).all(|(a, e)| {
        if e.is_nan() {
            a.is_nan()
        } else if e.is_infinite() || zeros && e == 0.0 {
            a.to_bits() == e.to_bits()
        } else {
            a.is_finite() && !(zeros && a == 0.0)
        }
    })
}

#[derive(Debug)]
struct TestCase<T: TestOutput> {
    pub r#in: Vec<f64>,
    pub out: T,
    /// The check of each output
    pub checks: Vec<Check>,
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

/// A table of `tests/data/xsref`, whose format is described in the README.md there
#[derive(Debug)]
struct Table {
    n_in: usize,
    n_out: usize,
    n_tol: usize,
    /// The number of `f64`s per row
    width: usize,
    /// The inputs, outputs, tolerances, and scales (if any) of each row (row-major)
    values: Vec<f64>,
    /// The check of each output of each row (row-major)
    checks: Vec<Check>,
}

impl Table {
    const MAGIC: &[u8] = b"XSFT";
    const VERSION: u32 = 1;
    const HEADER_LEN: usize = 20;
    /// The header flag for a status byte per logical output, after the row flags
    const STATUS: u16 = 1;
    /// The header flag for an error scale per logical output, after the tolerances
    const SCALES: u16 = 2;
    /// The row flag for a reference that is xsf's own value (xsref's `fallback`)
    const FALLBACK: u8 = 1;

    fn parse(bytes: &[u8]) -> Result<Self, TestError> {
        let (header, data) = bytes
            .split_at_checked(Self::HEADER_LEN)
            .ok_or(TestError::DataFormat)?;
        let u16_at = |i: usize| u16::from_le_bytes([header[i], header[i + 1]]);
        let u32_at = |i: usize| u32::from_le_bytes(header[i..i + 4].try_into().unwrap());

        let flags = u16_at(18);
        if &header[..4] != Self::MAGIC
            || u32_at(4) != Self::VERSION
            || flags & !(Self::STATUS | Self::SCALES) != 0
        {
            return Err(TestError::DataFormat);
        }
        let n_rows = u32_at(8) as usize;
        let [n_in, n_out, n_tol] = [12, 14, 16].map(|i| usize::from(u16_at(i)));
        let scaled = flags & Self::SCALES != 0;
        let n_status = if flags & Self::STATUS == 0 { 0 } else { n_tol };
        let width = n_in + n_out + n_tol * (1 + usize::from(scaled));

        // each row has its `f64` values, a flags byte, and the status bytes (if any)
        let len = n_rows
            .checked_mul(width * size_of::<f64>() + 1 + n_status)
            .and_then(|n| n.checked_add(Self::HEADER_LEN));
        if n_rows == 0 || n_tol == 0 || len != Some(bytes.len()) {
            return Err(TestError::DataFormat);
        }
        let (values, rest) = data.split_at(n_rows * width * size_of::<f64>());
        let (row_flags, status) = rest.split_at(n_rows);
        if row_flags.iter().any(|&flags| flags & !Self::FALLBACK != 0) {
            return Err(TestError::DataFormat);
        }

        let values = values
            .chunks_exact(size_of::<f64>())
            .map(|x| f64::from_le_bytes(x.try_into().unwrap()))
            .collect::<Vec<_>>();
        let mut checks = Vec::with_capacity(n_rows * n_tol);
        for (r, row) in values.chunks_exact(width).enumerate() {
            let (tols, scales) = row[n_in + n_out..].split_at(n_tol);
            for (i, &tol) in tols.iter().enumerate() {
                let check = Check {
                    tol,
                    scale: scales.get(i).copied().unwrap_or(0.0),
                    // without the status bytes, every reference is xsref's
                    status: status.get(r * n_tol + i).copied().unwrap_or(0).try_into()?,
                };
                if !check.is_valid() {
                    return Err(TestError::DataFormat);
                }
                checks.push(check);
            }
        }
        Ok(Self {
            n_in,
            n_out,
            n_tol,
            width,
            values,
            checks,
        })
    }

    /// The inputs, outputs, and checks of each row
    fn rows(&self) -> impl Iterator<Item = (&[f64], &[f64], &[Check])> {
        (self.values.chunks_exact(self.width))
            .zip(self.checks.chunks_exact(self.n_tol))
            .map(|(row, checks)| {
                let (inputs, rest) = row.split_at(self.n_in);
                (inputs, &rest[..self.n_out], checks)
            })
    }

    fn test_cases<T: TestOutput>(&self) -> Vec<TestCase<T>> {
        self.rows()
            .map(|(inputs, outputs, checks)| TestCase {
                r#in: inputs.to_vec(),
                out: T::from_row(outputs),
                checks: checks.to_vec(),
            })
            .collect()
    }
}

#[inline]
fn load_testcases<T: TestOutput>(name: &str, sig: &str) -> Result<Vec<TestCase<T>>, TestError> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../tests/data/xsref")
        .join(name)
        .join(format!("{sig}.bin"));
    Ok(Table::parse(&fs::read(path)?)?.test_cases())
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

/// The smallest relative tolerance of the legacy output `i` of a test case, for where xsref's
/// reference is less accurate than its tolerance; infinite if the reference is wrong
fn rtol_min(name: &str, x: &[f64], i: usize) -> f64 {
    match name {
        // https://github.com/scipy/xsf/pull/85#issuecomment-3738637926
        "hyp2f1" if x[0] == -1.0 => f64::INFINITY,
        // xsf computes the wrong imaginary part of F+-(x) for x < 0, which is fixed in xsf-rust
        "modified_fresnel_plus" | "modified_fresnel_minus" if x[0] < 0.0 => f64::INFINITY,
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
        // xsf overflows early for a complex z = x[1] + x[2] i
        "cyl_bessel_i" if x.len() == 3 && BESSEL_I_OVERFLOWS.contains(&x[0]) => f64::INFINITY,
        _ => 0.0,
    }
}

/// The smallest relative tolerance of a test case where the function is ill-conditioned, for
/// any reference
fn rtol_ill_conditioned(name: &str, x: &[f64]) -> f64 {
    match name {
        // for small a
        "gammaincinv" | "gammainccinv" if (0.0..0.01).contains(&x[0]) => 4e-12,
        // U(a, 1, x) near its zeros
        "hyperu" if x[1] == 1.0 && HYPERU_ZEROS.contains(&(x[0], x[2])) => 1e-12,
        // sin(pi v) of the reflection for large |v|
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
    let cases = load_testcases::<T>(name, signature).unwrap();
    if let Err(message) = check(name, &cases, test_fn) {
        panic!("{message}");
    }
}

/// Whether `test_fn` passes the test cases of the table `name`, or else the worst failure
fn check<T, F>(name: &str, cases: &[TestCase<T>], test_fn: F) -> Result<(), String>
where
    T: TestOutput,
    F: Fn(&[f64]) -> T,
{
    let mut failed = 0;
    let mut worst_err = -1.0;
    let mut worst_in = Vec::new();
    let mut worst_actual = String::new();
    let mut worst_expect = String::new();

    for case in cases {
        let actual = test_fn(&case.r#in);
        let values = actual.values();
        assert_eq!(
            values.len(),
            case.checks.len(),
            "{name}: a tolerance per output"
        );

        // like in xsf, each output has its own tolerance
        let error = (values
            .iter()
            .zip(case.out.values())
            .zip(&case.checks)
            .enumerate())
        .filter_map(|(i, ((&actual, expected), check))| {
            let mut rtol = rtol_ill_conditioned(name, &case.r#in);
            if check.status == Status::Legacy {
                rtol = rtol.max(rtol_min(name, &case.r#in, i));
            }
            check.failure(actual, expected, rtol)
        })
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
    if failed == 0 {
        Ok(())
    } else {
        Err(format!(
            "{name}: {failed}/{total} tests failed, worst case for {worst_in:?} \
            (actual={worst_actual}; desired={worst_expect}; error={worst_err:.2e})",
        ))
    }
}

#[cfg(test)]
mod tests {
    use num_complex::{Complex, c64};

    use super::{Check, Status, Table, TestError, TestOutput, check, table_error};

    /// A file with the dimensions `n_in, n_out, n_tol`, the header `flags`, and `rows` of
    /// `(values, row flags, statuses)`
    fn file(dims: [u16; 3], flags: u16, rows: &[(&[f64], u8, &[u8])]) -> Vec<u8> {
        let mut bytes = Table::MAGIC.to_vec();
        bytes.extend(Table::VERSION.to_le_bytes());
        bytes.extend(u32::try_from(rows.len()).unwrap().to_le_bytes());
        bytes.extend(dims.into_iter().chain([flags]).flat_map(u16::to_le_bytes));
        for (row, _, _) in rows {
            bytes.extend(row.iter().flat_map(|x| x.to_le_bytes()));
        }
        bytes.extend(rows.iter().map(|&(_, flags, _)| flags));
        bytes.extend(rows.iter().flat_map(|&(_, _, status)| status));
        bytes
    }

    /// Whether `test_fn` passes a one-row table with an input, `expected`, and its check
    fn passes<T: TestOutput>(name: &str, expected: &[f64], check: Check, test_fn: T) -> bool {
        let row = [&[1.0], expected, &[check.tol, check.scale]].concat();
        let n_out = u16::try_from(expected.len()).unwrap();
        let flags = Table::STATUS | Table::SCALES;
        let table = file([1, n_out, 1], flags, &[(&row, 0, &[check.status as u8])]);
        let cases = Table::parse(&table).unwrap().test_cases::<T>();
        super::check(name, &cases, |_| test_fn).is_ok()
    }

    #[test]
    fn test_table_singleton_bits() {
        // one row: an input, two complex outputs, and their tolerances
        let nan = f64::from_bits(0xfff8_0000_dead_beef);
        let row = [-0.0, 1.5, nan, -2.0, f64::INFINITY, 1e-16, 2.5e-13];
        let table = Table::parse(&file([1, 4, 2], 0, &[(&row, Table::FALLBACK, &[])])).unwrap();

        let rows = table.rows().collect::<Vec<_>>();
        assert_eq!(rows.len(), 1);
        let (inputs, outputs, checks) = rows[0];
        let bits = |x: &[f64]| x.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(inputs), bits(&row[..1]));
        assert_eq!(bits(outputs), bits(&row[1..5]));
        assert_eq!(
            bits(&checks.iter().map(|c| c.tol).collect::<Vec<_>>()),
            bits(&row[5..])
        );
        // without the optional sections, the outputs are xsref's and have no scale
        assert!(
            checks
                .iter()
                .all(|c| c.status == Status::Legacy && c.scale == 0.0)
        );
    }

    #[test]
    fn test_table_sections() {
        // an input, two outputs, their tolerances, and their scales
        let rows = [
            [1.0, 2.0, 3.0, 1e-16, f64::INFINITY, 0.0, 0.0],
            [4.0, f64::NAN, 5.0, 0.0, 2e-15, 0.0, 8.0],
        ];
        let status = [
            Status::Computed,
            Status::Legacy,
            Status::Convention,
            Status::Computed,
        ];
        for flags in 0..4 {
            let (has_status, has_scales) = (flags & Table::STATUS != 0, flags == 3);
            let n = if flags & Table::SCALES == 0 { 5 } else { 7 };
            let rows = rows.map(|mut row| {
                // without a status, the outputs are legacy ones, which have no scale
                row[6] = if has_scales { row[6] } else { 0.0 };
                row
            });
            let bytes = status.map(|s| s as u8);
            let (status0, status1) = if has_status {
                bytes.split_at(2)
            } else {
                Default::default()
            };
            let bytes = file(
                [1, 2, 2],
                flags,
                &[(&rows[0][..n], 0, status0), (&rows[1][..n], 1, status1)],
            );
            assert_eq!(bytes.len(), 20 + 2 * (8 * n + 1 + status0.len()));

            let cases = Table::parse(&bytes).unwrap().test_cases::<(f64, f64)>();
            assert_eq!((&cases[0].r#in, cases[0].out), (&vec![1.0], (2.0, 3.0)));
            let checks = [
                cases[0].checks[0],
                cases[0].checks[1],
                cases[1].checks[0],
                cases[1].checks[1],
            ];
            assert_eq!(checks.map(|c| c.tol), [1e-16, f64::INFINITY, 0.0, 2e-15]);
            assert_eq!(
                checks.map(|c| c.scale),
                [0.0, 0.0, 0.0, if has_scales { 8.0 } else { 0.0 }]
            );
            assert_eq!(
                checks.map(|c| c.status),
                if has_status {
                    status
                } else {
                    [Status::Legacy; 4]
                }
            );
        }
    }

    #[test]
    fn test_table_invalid() {
        let valid = file(
            [1, 1, 1],
            0,
            &[(&[1.0, 2.0, 0.0], 0, &[]), (&[3.0, 4.0, 0.0], 0, &[])],
        );
        assert!(Table::parse(&valid).is_ok());

        let invalid = |edit: fn(&mut Vec<u8>)| {
            let mut bytes = valid.clone();
            edit(&mut bytes);
            matches!(Table::parse(&bytes), Err(TestError::DataFormat))
        };
        // truncated, wrong magic, unknown version, unknown header flags, unknown row flags
        assert!(invalid(|b| _ = b.pop()));
        assert!(invalid(|b| b.truncate(Table::HEADER_LEN - 1)));
        assert!(invalid(|b| b[0] = b'Y'));
        assert!(invalid(|b| b[4] = 2));
        assert!(invalid(|b| b[18] = 4));
        assert!(invalid(|b| *b.last_mut().unwrap() = 2));
        // a section that the header flags announce is missing
        assert!(invalid(|b| b[18] = 1));
        assert!(invalid(|b| b[18] = 2));
        assert!(invalid(|b| b[18] = 3));
    }

    #[test]
    fn test_table_invalid_sections() {
        // an input, an output, its tolerance, and its scale; a computed and a legacy output
        let flags = Table::STATUS | Table::SCALES;
        let valid = |row0: [f64; 4], status0: u8| {
            let rows: [(&[f64], u8, &[u8]); 2] = [
                (&row0, 0, &[status0]),
                (&[3.0, 4.0, f64::INFINITY, 0.0], 0, &[0]),
            ];
            Table::parse(&file([1, 1, 1], flags, &rows)).is_ok()
        };
        assert!(valid([1.0, 2.0, 1e-15, 0.0], 1));
        assert!(valid([1.0, 2.0, 0.0, 5.0], 2));
        // an unknown status
        assert!(!valid([1.0, 2.0, 1e-15, 0.0], 3));
        // a scale that is negative, infinite, or NaN, or that a legacy output has
        assert!(!valid([1.0, 2.0, 1e-15, -1.0], 1));
        assert!(!valid([1.0, 2.0, 1e-15, f64::INFINITY], 1));
        assert!(!valid([1.0, 2.0, 1e-15, f64::NAN], 1));
        assert!(!valid([1.0, 2.0, 1e-15, 5.0], 0));
        // a tolerance of a computed output or a convention that is negative, infinite, or NaN
        assert!(!valid([1.0, 2.0, -1e-15, 0.0], 1));
        assert!(!valid([1.0, 2.0, f64::INFINITY, 0.0], 2));
        assert!(!valid([1.0, 2.0, f64::NAN, 0.0], 1));

        // truncated in the status bytes, or with one too many
        let mut bytes = file([1, 1, 1], flags, &[(&[1.0, 2.0, 0.0, 0.0], 0, &[1])]);
        assert!(Table::parse(&bytes).is_ok());
        bytes.push(1);
        assert!(Table::parse(&bytes).is_err());
        bytes.truncate(bytes.len() - 2);
        assert!(Table::parse(&bytes).is_err());
    }

    #[test]
    fn test_metric_vectors() {
        // shared with `tools/tables/generate.py`
        let number = |s: &str| s.parse::<f64>().unwrap();
        for line in include_str!("../metric_vectors.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let [actual, expected] =
                [0, 1].map(|i| fields[i].split(',').map(number).collect::<Vec<_>>());
            let scale = number(fields[2]);
            let error = if actual.len() == 2 {
                table_error(
                    c64(actual[0], actual[1]),
                    c64(expected[0], expected[1]),
                    scale,
                )
            } else {
                table_error(actual[0], expected[0], scale)
            };
            assert_eq!(error, number(fields[3]), "{line}");
        }
    }

    #[test]
    fn test_class() {
        let check = |status| Check {
            tol: 0.0,
            scale: 0.0,
            status,
        };
        let (inf, max) = (f64::INFINITY, f64::MAX);
        for status in [Status::Computed, Status::Convention] {
            assert!(passes("f", &[inf], check(status), inf));
            assert!(passes("f", &[0.0], check(status), 0.0));
            assert!(passes("f", &[f64::NAN], check(status), -f64::NAN));
            assert!(!passes("f", &[inf], check(status), -inf));
            assert!(!passes("f", &[f64::NAN], check(status), 1.0));
            // the error can't tell `f64::MAX` from +inf
            assert!(!passes("f", &[inf], check(status), max));
            assert!(!passes("f", &[1.0, -inf], check(status), c64(1.0, -max)));
            // nor -0 from +0, which only a convention checks
            let convention = status == Status::Convention;
            assert_eq!(passes("f", &[0.0], check(status), -0.0), !convention);
            assert_eq!(passes("f", &[-0.0], check(status), 0.0), !convention);
            assert_eq!(
                passes("f", &[inf, 0.0], check(status), c64(inf, -0.0)),
                !convention
            );
        }
        // a legacy output only has its error
        assert!(passes("f", &[inf], check(Status::Legacy), max));
        // a finite convention has its tolerance
        let convention = Check {
            tol: 1e-15,
            scale: 0.0,
            status: Status::Convention,
        };
        assert!(passes("f", &[1.0], convention, 1.0 + 2e-15));
        assert!(!passes("f", &[1.0], convention, 1.0 + 8e-15));
        assert!(!passes("f", &[1.0], convention, 0.0));
    }

    #[test]
    fn test_scale() {
        // near a zero, where the error is relative to the scale
        let check = |scale| Check {
            tol: 1e-15,
            scale,
            status: Status::Computed,
        };
        let (expected, actual) = (1e-3, 1e-3 + 1e-15);
        assert!(!passes("f", &[expected], check(0.0), actual));
        assert!(passes("f", &[expected], check(1.0), actual));
        assert!(!passes("f", &[expected], check(1.0), actual + 1e-14));
        // a scale below the reference changes nothing
        assert!(!passes("f", &[expected], check(1e-6), actual));
        assert!(!passes("f", &[expected], check(1.0), f64::NAN));
    }

    #[test]
    fn test_rtol_min() {
        // the input `x`, and two outputs with tolerances of 0 and the statuses `status`
        let cases = |x: f64, status: [u8; 2]| {
            let row = [x, 1.0, 1.0, 0.0, 0.0];
            let table = file([1, 2, 2], Table::STATUS, &[(&row, 0, &status)]);
            Table::parse(&table).unwrap().test_cases::<(f64, f64)>()
        };
        // the rule for xsref's second output of `pbwa` only applies to that legacy output
        let actual = |_: &[f64]| (1.0, 1.0 + 5e-12);
        assert!(check("pbwa", &cases(1.0, [0, 0]), actual).is_ok());
        assert!(check("pbwa", &cases(1.0, [1, 0]), actual).is_ok());
        assert!(check("pbwa", &cases(1.0, [0, 1]), actual).is_err());
        assert!(check("pbwa", &cases(1.0, [0, 0]), |_| (1.0 + 5e-12, 1.0)).is_err());
        assert!(check("f", &cases(1.0, [0, 0]), actual).is_err());
        // the rule for the ill-conditioned `gammaincinv` at small a applies to any output
        let actual = |_: &[f64]| (1.0 + 2e-12, 1.0 + 2e-12);
        assert!(check("gammaincinv", &cases(0.0, [0, 0]), actual).is_ok());
        assert!(check("gammaincinv", &cases(0.0, [1, 2]), actual).is_ok());
        assert!(check("gammaincinv", &cases(1.0, [1, 2]), actual).is_err());
    }

    #[test]
    fn test_complex_parts() {
        // a complex output is one logical output
        let check = Check {
            tol: 0.0,
            scale: 0.0,
            status: Status::Convention,
        };
        let nan = Complex::new(f64::NAN, f64::NAN);
        assert!(passes("f", &[f64::NAN, f64::NAN], check, nan));
        assert!(!passes(
            "f",
            &[f64::NAN, f64::NAN],
            check,
            c64(f64::NAN, 0.0)
        ));
    }
}
