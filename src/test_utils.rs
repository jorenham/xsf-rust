#[cfg(test)]
pub mod xsref {
    use arrow::array::{Array, Float64Array, Int32Array, Int64Array};
    use num_complex::{Complex, c64};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use std::env;
    use std::fs::File;
    use std::path::PathBuf;

    pub trait TestOutputValue: Copy {
        fn error(self, expected: Self) -> f64;
        fn magnitude(self) -> f64;
        fn format(self) -> String;
    }

    impl TestOutputValue for f64 {
        fn error(self, expected: Self) -> f64 {
            relative_error(self, expected)
        }

        fn magnitude(self) -> f64 {
            self.abs()
        }

        fn format(self) -> String {
            format!("{:.6e}", self)
        }
    }

    impl TestOutputValue for Complex<f64> {
        fn error(self, expected: Self) -> f64 {
            complex_relative_error(self, expected)
        }

        fn magnitude(self) -> f64 {
            self.norm()
        }

        fn format(self) -> String {
            format!("{:.6e}+{:.6e}i", self.re, self.im)
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
                .map(|(&a, &e)| a.error(e))
                .fold(0.0, |acc, x| acc.max(x))
        }

        fn format(self) -> String {
            let values = self.values();
            if values.len() == 1 {
                values[0].format()
            } else {
                format!(
                    "({})",
                    values
                        .iter()
                        .map(|x| x.format())
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
        let tables_dir = env::var("XSREF_TABLES_PATH")
            .map(PathBuf::from)
            .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("xsref/tables"))
            .join("scipy_special_tests")
            .join(name);

        let get_table = |name: &str| File::open(tables_dir.join(format!("{}.parquet", name)));

        let table_in = read_parquet_rows(get_table(&format!("In_{}", signature))?);
        let table_out = read_parquet_output::<T>(get_table(&format!("Out_{}", signature))?);

        let platform = if cfg!(all(target_arch = "x86_64", target_os = "linux")) {
            "gcc-linux-x86_64"
        } else if cfg!(all(target_arch = "aarch64", target_os = "macos")) {
            "clang-darwin-aarch64"
        } else {
            "other"
        };

        let table_err_file = get_table(&format!("Err_{}_{}", signature, platform))
            .or_else(|_| get_table(&format!("Err_{}_other", signature)))?;
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

    /// based on `xsref.float_tools.extended_absolute_error`
    fn absolute_error(actual: f64, desired: f64) -> f64 {
        if actual == desired || (actual.is_nan() && desired.is_nan()) {
            return 0.0;
        } else if desired.is_nan() || actual.is_nan() {
            return f64::INFINITY;
        }

        // f64::MAX + ULP
        let max1p = f64::MAX * (1.0 + 2.0_f64.powi(-(f64::MANTISSA_DIGITS as i32)));
        if actual.is_infinite() {
            (actual.signum() * max1p - desired).abs()
        } else if desired.is_infinite() {
            (actual - desired.signum() * max1p).abs()
        } else {
            (actual - desired).abs()
        }
    }

    /// based on `xsref.float_tools.extended_relative_error`
    fn relative_error(actual: f64, desired: f64) -> f64 {
        let abserr0 = if desired == 0.0 {
            f64::MIN_POSITIVE
        } else if desired.is_infinite() {
            f64::MAX
        } else if desired.is_nan() {
            1.0
        } else {
            desired.abs()
        };

        let abserr = absolute_error(actual, desired);
        (abserr / abserr0).min(abserr)
    }

    /// based on `xsref.float_tools.extended_absolute_error_complex`
    fn complex_absolute_error(actual: Complex<f64>, expected: Complex<f64>) -> f64 {
        let real_error = (actual.re - expected.re).abs();
        let imag_error = (actual.im - expected.im).abs();

        let max_error = real_error.max(imag_error);
        if max_error == 0.0 {
            0.0
        } else if max_error.is_infinite() {
            f64::INFINITY
        } else {
            real_error.hypot(imag_error)
        }
    }

    /// based on `xsref.float_tools.extended_relative_error_complex`
    fn complex_relative_error(actual: Complex<f64>, expected: Complex<f64>) -> f64 {
        if actual.re.is_nan() || actual.im.is_nan() {
            if expected.re.is_nan() || expected.im.is_nan() {
                0.0
            } else {
                f64::INFINITY
            }
        } else if expected.re.is_infinite() || expected.im.is_infinite() {
            // If any component of expected is infinite, use component-wise relative error
            fn rel_error(x: f64, y: f64) -> f64 {
                if y.is_infinite() {
                    if x.is_infinite() && x.signum() == y.signum() {
                        0.0
                    } else {
                        f64::INFINITY
                    }
                } else if y == 0.0 {
                    x.abs()
                } else {
                    (x / y - 1.0).abs()
                }
            }

            let real_rel_err = rel_error(actual.re, expected.re);
            let imag_rel_err = rel_error(actual.im, expected.im);
            real_rel_err.max(imag_rel_err)
        } else {
            // For finite expected values, use the magnitude-based relative error
            let abs_error = complex_absolute_error(actual, expected);
            let expected_magnitude = expected.re.hypot(expected.im);
            if expected_magnitude == 0.0 {
                abs_error
            } else {
                abs_error / expected_magnitude
            }
        }
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
}
