use num_complex::{Complex, c64};

mod xsref {
    use super::*;
    use arrow::array::{
        Array, Float32Array as ArrayF32, Float64Array as ArrayF64, Int32Array as ArrayI32,
        Int64Array as ArrayI64,
    };
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use std::env;
    use std::fs::File;
    use std::path::PathBuf;

    pub trait TestOutput: Copy {
        fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self>;
        fn is_huge(self) -> bool {
            self.magnitude() > f64::EPSILON
        }
        fn error(actual: Self, expected: Self) -> f64;
        fn magnitude(self) -> f64;
        fn format_value(self) -> String;
    }

    impl TestOutput for f64 {
        fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
            rows.into_iter().map(|row| row[0]).collect()
        }

        fn error(actual: Self, expected: Self) -> f64 {
            relative_error(actual, expected)
        }

        fn magnitude(self) -> f64 {
            self.abs()
        }

        fn format_value(self) -> String {
            format!("{:.6e}", self)
        }
    }

    impl TestOutput for Complex<f64> {
        fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
            rows.into_iter().map(|row| c64(row[0], row[1])).collect()
        }

        fn error(actual: Self, expected: Self) -> f64 {
            complex_relative_error(actual, expected)
        }

        fn magnitude(self) -> f64 {
            self.norm()
        }

        fn format_value(self) -> String {
            format!("{:.6e}+{:.6e}i", self.re, self.im)
        }
    }

    impl TestOutput for (f64, f64) {
        fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
            rows.into_iter().map(|row| (row[0], row[1])).collect()
        }

        fn error(actual: Self, expected: Self) -> f64 {
            let errors = [
                relative_error(actual.0, expected.0),
                relative_error(actual.1, expected.1),
            ];
            errors.iter().fold(0.0, |acc, &x| acc.max(x))
        }

        fn magnitude(self) -> f64 {
            (self.0.abs() + self.1.abs()) / 2.0
        }

        fn format_value(self) -> String {
            format!("({:.6e}, {:.6e})", self.0, self.1)
        }
    }

    impl TestOutput for (f64, f64, f64, f64) {
        fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
            rows.into_iter()
                .map(|row| (row[0], row[1], row[2], row[3]))
                .collect()
        }

        fn error(actual: Self, expected: Self) -> f64 {
            let errors = [
                relative_error(actual.0, expected.0),
                relative_error(actual.1, expected.1),
                relative_error(actual.2, expected.2),
                relative_error(actual.3, expected.3),
            ];
            errors.iter().fold(0.0, |acc, &x| acc.max(x))
        }

        fn magnitude(self) -> f64 {
            (self.0.abs() + self.1.abs() + self.2.abs() + self.3.abs()) / 4.0
        }

        fn format_value(self) -> String {
            format!(
                "({:.6e}, {:.6e}, {:.6e}, {:.6e})",
                self.0, self.1, self.2, self.3
            )
        }
    }

    impl TestOutput for (Complex<f64>, Complex<f64>) {
        fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
            rows.into_iter()
                .map(|row| (c64(row[0], row[1]), c64(row[2], row[3])))
                .collect()
        }

        fn error(actual: Self, expected: Self) -> f64 {
            let errors = [
                complex_relative_error(actual.0, expected.0),
                complex_relative_error(actual.1, expected.1),
            ];
            errors.iter().fold(0.0, |acc, &x| acc.max(x))
        }

        fn magnitude(self) -> f64 {
            (self.0.norm() + self.1.norm()) / 2.0
        }

        fn format_value(self) -> String {
            format!(
                "({:.6e}+{:.6e}i, {:.6e}+{:.6e}i)",
                self.0.re, self.0.im, self.1.re, self.1.im,
            )
        }
    }

    impl TestOutput for (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>) {
        fn from_parquet_rows(rows: Vec<Vec<f64>>) -> Vec<Self> {
            rows.into_iter()
                .map(|row| {
                    (
                        c64(row[0], row[1]),
                        c64(row[2], row[3]),
                        c64(row[4], row[5]),
                        c64(row[6], row[7]),
                    )
                })
                .collect()
        }

        fn error(actual: Self, expected: Self) -> f64 {
            let errors = [
                complex_relative_error(actual.0, expected.0),
                complex_relative_error(actual.1, expected.1),
                complex_relative_error(actual.2, expected.2),
                complex_relative_error(actual.3, expected.3),
            ];
            errors.iter().fold(0.0, |acc, &x| acc.max(x))
        }

        fn magnitude(self) -> f64 {
            (self.0.norm() + self.1.norm() + self.2.norm() + self.3.norm()) / 4.0
        }

        fn format_value(self) -> String {
            format!(
                "({:.6e}+{:.6e}i, {:.6e}+{:.6e}i, {:.6e}+{:.6e}i, {:.6e}+{:.6e}i)",
                self.0.re,
                self.0.im,
                self.1.re,
                self.1.im,
                self.2.re,
                self.2.im,
                self.3.re,
                self.3.im
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
                    if let Some(f64_array) = column.downcast_ref::<ArrayF64>() {
                        row.push(f64_array.value(row_idx));
                    } else if let Some(f32_array) = column.downcast_ref::<ArrayF32>() {
                        row.push(f32_array.value(row_idx) as f64);
                    } else if let Some(i32_array) = column.downcast_ref::<ArrayI32>() {
                        row.push(i32_array.value(row_idx) as f64);
                    } else if let Some(i64_array) = column.downcast_ref::<ArrayI64>() {
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
        for batch_result in reader {
            let batch = batch_result.unwrap();
            let col0 = batch.column(0).as_any();

            if let Some(column) = col0.downcast_ref::<ArrayF64>() {
                for i in 0..column.len() {
                    values.push(column.value(i));
                }
            } else if let Some(column) = col0.downcast_ref::<ArrayF32>() {
                for i in 0..column.len() {
                    values.push(column.value(i) as f64);
                }
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
                    let column = batch.column(col_idx).as_any();
                    if let Some(f64_array) = column.downcast_ref::<ArrayF64>() {
                        row.push(f64_array.value(row_idx));
                    } else if let Some(f32_array) = column.downcast_ref::<ArrayF32>() {
                        row.push(f32_array.value(row_idx) as f64);
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
            .unwrap_or_else(|_| {
                PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/xsref/tables")
            })
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
                    desired.format_value(),
                    actual.format_value(),
                    error,
                    max_error
                );
            }
        }

        let tested = cases.len();
        assert!(failed == 0, "{}: {}/{} tests failed", name, failed, tested);
        println!("{}: {}/{} tests passed", name, tested - failed, tested);
    }
}

#[test]
fn test_airy_f64() {
    xsref::test::<(f64, f64, f64, f64), _>("airy", "d-d_d_d_d", |x: &[f64]| xsf::airy(x[0]));
}

#[test]
fn test_airy_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>), _>(
        "airy",
        "cd-cd_cd_cd_cd",
        |x: &[f64]| xsf::airy(c64(x[0], x[1])),
    );
}

#[test]
fn test_airye_f64() {
    xsref::test::<(f64, f64, f64, f64), _>("airye", "d-d_d_d_d", |x: &[f64]| xsf::airye(x[0]));
}

#[test]
fn test_airye_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>), _>(
        "airye",
        "cd-cd_cd_cd_cd",
        |x: &[f64]| xsf::airye(c64(x[0], x[1])),
    );
}

#[test]
fn test_itairy_f64() {
    xsref::test::<(f64, f64, f64, f64), _>("itairy", "d-d_d_d_d", |x: &[f64]| xsf::itairy(x[0]));
}

#[test]
fn test_cbrt_f64() {
    xsref::test::<f64, _>("cbrt", "d-d", |x: &[f64]| xsf::cbrt(x[0]));
}

#[test]
fn test_cyl_bessel_j0_f64() {
    xsref::test::<f64, _>("cyl_bessel_j0", "d-d", |x: &[f64]| xsf::cyl_bessel_j0(x[0]));
}

#[test]
fn test_cyl_bessel_j1_f64() {
    xsref::test::<f64, _>("cyl_bessel_j1", "d-d", |x: &[f64]| xsf::cyl_bessel_j1(x[0]));
}

#[test]
fn test_cyl_bessel_y0_f64() {
    xsref::test::<f64, _>("cyl_bessel_y0", "d-d", |x: &[f64]| xsf::cyl_bessel_y0(x[0]));
}

#[test]
fn test_cyl_bessel_y1_f64() {
    xsref::test::<f64, _>("cyl_bessel_y1", "d-d", |x: &[f64]| xsf::cyl_bessel_y1(x[0]));
}

#[test]
fn test_cyl_bessel_i0_f64() {
    xsref::test::<f64, _>("cyl_bessel_i0", "d-d", |x: &[f64]| xsf::cyl_bessel_i0(x[0]));
}

#[test]
fn test_cyl_bessel_i0e_f64() {
    xsref::test::<f64, _>("cyl_bessel_i0e", "d-d", |x: &[f64]| {
        xsf::cyl_bessel_i0e(x[0])
    });
}

#[test]
fn test_cyl_bessel_i1_f64() {
    xsref::test::<f64, _>("cyl_bessel_i1", "d-d", |x: &[f64]| xsf::cyl_bessel_i1(x[0]));
}

#[test]
fn test_cyl_bessel_i1e_f64() {
    xsref::test::<f64, _>("cyl_bessel_i1e", "d-d", |x: &[f64]| {
        xsf::cyl_bessel_i1e(x[0])
    });
}

#[test]
fn test_cyl_bessel_k0_f64() {
    xsref::test::<f64, _>("cyl_bessel_k0", "d-d", |x: &[f64]| xsf::cyl_bessel_k0(x[0]));
}

#[test]
fn test_cyl_bessel_k0e_f64() {
    xsref::test::<f64, _>("cyl_bessel_k0e", "d-d", |x: &[f64]| {
        xsf::cyl_bessel_k0e(x[0])
    });
}

#[test]
fn test_cyl_bessel_k1_f64() {
    xsref::test::<f64, _>("cyl_bessel_k1", "d-d", |x: &[f64]| xsf::cyl_bessel_k1(x[0]));
}

#[test]
fn test_cyl_bessel_k1e_f64() {
    xsref::test::<f64, _>("cyl_bessel_k1e", "d-d", |x: &[f64]| {
        xsf::cyl_bessel_k1e(x[0])
    });
}
#[test]
fn test_cyl_bessel_j_f64() {
    xsref::test::<f64, _>("cyl_bessel_j", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_j(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_j_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_j", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_j(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_bessel_je_f64() {
    xsref::test::<f64, _>("cyl_bessel_je", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_je(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_je_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_je", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_je(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_bessel_y_f64() {
    xsref::test::<f64, _>("cyl_bessel_y", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_y(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_y_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_y", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_y(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_bessel_ye_f64() {
    xsref::test::<f64, _>("cyl_bessel_ye", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_ye(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_ye_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_ye", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_ye(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_bessel_i_f64() {
    xsref::test::<f64, _>("cyl_bessel_i", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_i(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_i_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_i", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_i(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_bessel_ie_f64() {
    xsref::test::<f64, _>("cyl_bessel_ie", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_ie(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_ie_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_ie", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_ie(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_bessel_k_f64() {
    xsref::test::<f64, _>("cyl_bessel_k", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_k(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_k_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_k", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_k(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_bessel_ke_f64() {
    xsref::test::<f64, _>("cyl_bessel_ke", "d_d-d", |x: &[f64]| {
        xsf::cyl_bessel_ke(x[0], x[1])
    });
}

#[test]
fn test_cyl_bessel_ke_c64() {
    xsref::test::<Complex<f64>, _>("cyl_bessel_ke", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_bessel_ke(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_hankel_1_c64() {
    xsref::test::<Complex<f64>, _>("cyl_hankel_1", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_hankel_1(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_hankel_1e_c64() {
    xsref::test::<Complex<f64>, _>("cyl_hankel_1e", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_hankel_1e(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_hankel_2_c64() {
    xsref::test::<Complex<f64>, _>("cyl_hankel_2", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_hankel_2(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_cyl_hankel_2e_c64() {
    xsref::test::<Complex<f64>, _>("cyl_hankel_2e", "d_cd-cd", |x: &[f64]| {
        xsf::cyl_hankel_2e(x[0], c64(x[1], x[2]))
    });
}

#[test]
fn test_besselpoly_f64() {
    xsref::test::<f64, _>("besselpoly", "d_d_d-d", |x: &[f64]| {
        xsf::besselpoly(x[0], x[1], x[2])
    });
}

#[test]
fn test_it1j0y0_f64() {
    xsref::test::<(f64, f64), _>("it1j0y0", "d-d_d", |x: &[f64]| xsf::it1j0y0(x[0]));
}

#[test]
fn test_it2j0y0_f64() {
    xsref::test::<(f64, f64), _>("it2j0y0", "d-d_d", |x: &[f64]| xsf::it2j0y0(x[0]));
}

#[test]
fn test_it1i0k0_f64() {
    xsref::test::<(f64, f64), _>("it1i0k0", "d-d_d", |x: &[f64]| xsf::it1i0k0(x[0]));
}

#[test]
fn test_it2i0k0_f64() {
    xsref::test::<(f64, f64), _>("it2i0k0", "d-d_d", |x: &[f64]| xsf::it2i0k0(x[0]));
}

#[test]
fn test_beta_f64() {
    xsref::test::<f64, _>("beta", "d_d-d", |x: &[f64]| xsf::beta(x[0], x[1]));
}

#[test]
fn test_betaln_f64() {
    xsref::test::<f64, _>("betaln", "d_d-d", |x: &[f64]| xsf::betaln(x[0], x[1]));
}

#[test]
fn test_binom_f64() {
    xsref::test::<f64, _>("binom", "d_d-d", |x: &[f64]| xsf::binom(x[0], x[1]));
}

#[test]
fn test_gdtrib_f64() {
    xsref::test::<f64, _>("gdtrib", "d_d_d-d", |x: &[f64]| {
        xsf::gdtrib(x[0], x[1], x[2])
    });
}

#[test]
fn test_digamma_f64() {
    xsref::test::<f64, _>("digamma", "d-d", |x: &[f64]| xsf::digamma(x[0]));
}

#[test]
fn test_digamma_c64() {
    xsref::test::<Complex<f64>, _>("digamma", "cd-cd", |x: &[f64]| {
        xsf::digamma(c64(x[0], x[1]))
    });
}

#[test]
fn test_ellipk_f64() {
    xsref::test::<f64, _>("ellipk", "d-d", |x: &[f64]| xsf::ellipk(x[0]));
}

#[test]
fn test_ellipkm1_f64() {
    xsref::test::<f64, _>("ellipkm1", "d-d", |x: &[f64]| xsf::ellipkm1(x[0]));
}

#[test]
fn test_ellipkinc_f64() {
    xsref::test::<f64, _>("ellipkinc", "d_d-d", |x: &[f64]| xsf::ellipkinc(x[0], x[1]));
}

#[test]
fn test_ellipe_f64() {
    xsref::test::<f64, _>("ellipe", "d-d", |x: &[f64]| xsf::ellipe(x[0]));
}

#[test]
fn test_ellipeinc_f64() {
    xsref::test::<f64, _>("ellipeinc", "d_d-d", |x: &[f64]| xsf::ellipeinc(x[0], x[1]));
}

#[test]
fn test_ellipj_f64() {
    xsref::test::<(f64, f64, f64, f64), _>("ellipj", "d_d-d_d_d_d", |x: &[f64]| {
        xsf::ellipj(x[0], x[1])
    });
}

#[test]
fn test_erf_f64() {
    xsref::test::<f64, _>("erf", "d-d", |x: &[f64]| xsf::erf(x[0]));
}

#[test]
fn test_erf_c64() {
    xsref::test::<Complex<f64>, _>("erf", "cd-cd", |x: &[f64]| xsf::erf(c64(x[0], x[1])));
}

#[test]
fn test_erfc_f64() {
    xsref::test::<f64, _>("erfc", "d-d", |x: &[f64]| xsf::erfc(x[0]));
}

#[test]
fn test_erfc_c64() {
    xsref::test::<Complex<f64>, _>("erfc", "cd-cd", |x: &[f64]| xsf::erfc(c64(x[0], x[1])));
}

#[test]
fn test_erfcx_f64() {
    xsref::test::<f64, _>("erfcx", "d-d", |x: &[f64]| xsf::erfcx(x[0]));
}

#[test]
fn test_erfcx_c64() {
    xsref::test::<Complex<f64>, _>("erfcx", "cd-cd", |x: &[f64]| xsf::erfcx(c64(x[0], x[1])));
}

#[test]
fn test_erfi_f64() {
    xsref::test::<f64, _>("erfi", "d-d", |x: &[f64]| xsf::erfi(x[0]));
}

#[test]
fn test_erfi_c64() {
    xsref::test::<Complex<f64>, _>("erfi", "cd-cd", |x: &[f64]| xsf::erfi(c64(x[0], x[1])));
}

#[test]
fn test_dawsn_f64() {
    xsref::test::<f64, _>("dawsn", "d-d", |x: &[f64]| xsf::dawsn(x[0]));
}

#[test]
fn test_dawsn_c64() {
    xsref::test::<Complex<f64>, _>("dawsn", "cd-cd", |x: &[f64]| xsf::dawsn(c64(x[0], x[1])));
}

#[test]
fn test_wofz_c64() {
    xsref::test::<Complex<f64>, _>("wofz", "cd-cd", |x: &[f64]| xsf::wofz(c64(x[0], x[1])));
}

#[test]
fn test_voigt_profile_f64() {
    xsref::test::<f64, _>("voigt_profile", "d_d_d-d", |x: &[f64]| {
        xsf::voigt_profile(x[0], x[1], x[2])
    });
}

#[test]
fn test_expm1_f64() {
    xsref::test::<f64, _>("expm1", "d-d", |x: &[f64]| xsf::expm1(x[0]));
}

#[test]
fn test_expm1_c64() {
    xsref::test::<Complex<f64>, _>("expm1", "cd-cd", |x: &[f64]| xsf::expm1(c64(x[0], x[1])));
}

#[test]
fn test_exp2_f64() {
    xsref::test::<f64, _>("exp2", "d-d", |x: &[f64]| xsf::exp2(x[0]));
}

#[test]
fn test_exp10_f64() {
    xsref::test::<f64, _>("exp10", "d-d", |x: &[f64]| xsf::exp10(x[0]));
}

#[test]
fn test_expi_f64() {
    xsref::test::<f64, _>("expi", "d-d", |x: &[f64]| xsf::expi(x[0]));
}

#[test]
fn test_expi_c64() {
    xsref::test::<Complex<f64>, _>("expi", "cd-cd", |x: &[f64]| xsf::expi(c64(x[0], x[1])));
}

#[test]
fn test_exp1_f64() {
    xsref::test::<f64, _>("exp1", "d-d", |x: &[f64]| xsf::exp1(x[0]));
}

#[test]
fn test_exp1_c64() {
    xsref::test::<Complex<f64>, _>("exp1", "cd-cd", |x: &[f64]| xsf::exp1(c64(x[0], x[1])));
}

#[test]
fn test_scaled_exp1_f64() {
    xsref::test::<f64, _>("scaled_exp1", "d-d", |x: &[f64]| xsf::scaled_exp1(x[0]));
}

#[test]
fn test_fresnel_f64() {
    xsref::test::<(f64, f64), _>("fresnel", "d-d_d", |x: &[f64]| xsf::fresnel(x[0]));
}

#[test]
fn test_fresnel_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>), _>("fresnel", "cd-cd_cd", |x: &[f64]| {
        xsf::fresnel(c64(x[0], x[1]))
    });
}

#[test]
fn test_modified_fresnel_plus_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>), _>(
        "modified_fresnel_plus",
        "d-cd_cd",
        |x: &[f64]| xsf::modified_fresnel_plus(x[0]),
    );
}

#[test]
fn test_modified_fresnel_minus_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>), _>(
        "modified_fresnel_minus",
        "d-cd_cd",
        |x: &[f64]| xsf::modified_fresnel_minus(x[0]),
    );
}

#[test]
fn test_gamma_f64() {
    xsref::test::<f64, _>("gamma", "d-d", |x: &[f64]| xsf::gamma(x[0]));
}

#[test]
fn test_gamma_c64() {
    xsref::test::<Complex<f64>, _>("gamma", "cd-cd", |x: &[f64]| xsf::gamma(c64(x[0], x[1])));
}

#[test]
fn test_gammainc_f64() {
    xsref::test::<f64, _>("gammainc", "d_d-d", |x: &[f64]| xsf::gammainc(x[0], x[1]));
}

#[test]
fn test_gammaincc_f64() {
    xsref::test::<f64, _>("gammaincc", "d_d-d", |x: &[f64]| xsf::gammaincc(x[0], x[1]));
}

#[test]
fn test_gammaincinv_f64() {
    xsref::test::<f64, _>("gammaincinv", "d_d-d", |x: &[f64]| {
        xsf::gammaincinv(x[0], x[1])
    });
}

#[test]
fn test_gammainccinv_f64() {
    xsref::test::<f64, _>("gammainccinv", "d_d-d", |x: &[f64]| {
        xsf::gammainccinv(x[0], x[1])
    });
}

#[test]
fn test_gammaln_f64() {
    xsref::test::<f64, _>("gammaln", "d-d", |x: &[f64]| xsf::gammaln(x[0]));
}

#[test]
fn test_gammasgn_f64() {
    xsref::test::<f64, _>("gammasgn", "d-d", |x: &[f64]| xsf::gammasgn(x[0]));
}

#[test]
fn test_hyp2f1_f64() {
    xsref::test::<f64, _>("hyp2f1", "d_d_d_d-d", |x: &[f64]| {
        xsf::hyp2f1(x[0], x[1], x[2], x[3])
    });
}

#[test]
fn test_hyp2f1_c64() {
    xsref::test::<Complex<f64>, _>("hyp2f1", "d_d_d_cd-cd", |x: &[f64]| {
        xsf::hyp2f1(x[0], x[1], x[2], c64(x[3], x[4]))
    });
}

#[test]
fn test_iv_ratio_f64() {
    xsref::test::<f64, _>("iv_ratio", "d_d-d", |x: &[f64]| xsf::iv_ratio(x[0], x[1]));
}

#[test]
fn test_iv_ratio_c_f64() {
    xsref::test::<f64, _>("iv_ratio_c", "d_d-d", |x: &[f64]| {
        xsf::iv_ratio_c(x[0], x[1])
    });
}

#[test]
fn test_ber_f64() {
    xsref::test::<f64, _>("ber", "d-d", |x: &[f64]| xsf::ber(x[0]));
}

#[test]
fn test_bei_f64() {
    xsref::test::<f64, _>("bei", "d-d", |x: &[f64]| xsf::bei(x[0]));
}

#[test]
fn test_ker_f64() {
    xsref::test::<f64, _>("ker", "d-d", |x: &[f64]| xsf::ker(x[0]));
}

#[test]
fn test_kei_f64() {
    xsref::test::<f64, _>("kei", "d-d", |x: &[f64]| xsf::kei(x[0]));
}

#[test]
fn test_berp_f64() {
    xsref::test::<f64, _>("berp", "d-d", |x: &[f64]| xsf::berp(x[0]));
}

#[test]
fn test_beip_f64() {
    xsref::test::<f64, _>("beip", "d-d", |x: &[f64]| xsf::beip(x[0]));
}

#[test]
fn test_kerp_f64() {
    xsref::test::<f64, _>("kerp", "d-d", |x: &[f64]| xsf::kerp(x[0]));
}

#[test]
fn test_keip_f64() {
    xsref::test::<f64, _>("keip", "d-d", |x: &[f64]| xsf::keip(x[0]));
}

#[test]
fn test_kelvin_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>), _>(
        "kelvin",
        "d-cd_cd_cd_cd",
        |x: &[f64]| xsf::kelvin(x[0]),
    );
}

#[test]
fn test_lambertw_c64() {
    xsref::test::<Complex<f64>, _>("lambertw", "cd_p_d-cd", |x: &[f64]| {
        xsf::lambertw(c64(x[0], x[1]), x[2] as std::os::raw::c_long, x[3])
    });
}

#[test]
fn test_expit_f64() {
    xsref::test::<f64, _>("expit", "d-d", |x: &[f64]| xsf::expit(x[0]));
}

#[test]
fn test_exprel_f64() {
    xsref::test::<f64, _>("exprel", "d-d", |x: &[f64]| xsf::exprel(x[0]));
}

#[test]
fn test_logit_f64() {
    xsref::test::<f64, _>("logit", "d-d", |x: &[f64]| xsf::logit(x[0]));
}

#[test]
fn test_log_expit_f64() {
    xsref::test::<f64, _>("log_expit", "d-d", |x: &[f64]| xsf::log_expit(x[0]));
}

#[test]
fn test_log1p_f64() {
    xsref::test::<f64, _>("log1p", "d-d", |x: &[f64]| xsf::log1p(x[0]));
}

#[test]
fn test_log1p_c64() {
    xsref::test::<Complex<f64>, _>("log1p", "cd-cd", |x: &[f64]| xsf::log1p(c64(x[0], x[1])));
}

#[test]
fn test_log1pmx_f64() {
    xsref::test::<f64, _>("log1pmx", "d-d", |x: &[f64]| xsf::log1pmx(x[0]));
}

#[test]
fn test_xlogy_f64() {
    xsref::test::<f64, _>("xlogy", "d_d-d", |x: &[f64]| xsf::xlogy(x[0], x[1]));
}

#[test]
fn test_xlogy_c64() {
    xsref::test::<Complex<f64>, _>("xlogy", "cd_cd-cd", |x: &[f64]| {
        xsf::xlogy(c64(x[0], x[1]), c64(x[2], x[3]))
    });
}

#[test]
fn test_xlog1py_f64() {
    xsref::test::<f64, _>("xlog1py", "d_d-d", |x: &[f64]| xsf::xlog1py(x[0], x[1]));
}

#[test]
fn test_xlog1py_c64() {
    xsref::test::<Complex<f64>, _>("xlog1py", "cd_cd-cd", |x: &[f64]| {
        xsf::xlog1py(c64(x[0], x[1]), c64(x[2], x[3]))
    });
}

#[test]
fn test_loggamma_f64() {
    xsref::test::<f64, _>("loggamma", "d-d", |x: &[f64]| xsf::loggamma(x[0]));
}

#[test]
fn test_loggamma_c64() {
    xsref::test::<Complex<f64>, _>("loggamma", "cd-cd", |x: &[f64]| {
        xsf::loggamma(c64(x[0], x[1]))
    });
}

#[test]
fn test_rgamma_f64() {
    xsref::test::<f64, _>("rgamma", "d-d", |x: &[f64]| xsf::rgamma(x[0]));
}

#[test]
fn test_rgamma_c64() {
    xsref::test::<Complex<f64>, _>("rgamma", "cd-cd", |x: &[f64]| xsf::rgamma(c64(x[0], x[1])));
}

#[test]
fn test_cem_cva_f64() {
    xsref::test::<f64, _>("cem_cva", "d_d-d", |x: &[f64]| xsf::cem_cva(x[0], x[1]));
}

#[test]
fn test_sem_cva_f64() {
    xsref::test::<f64, _>("sem_cva", "d_d-d", |x: &[f64]| xsf::sem_cva(x[0], x[1]));
}

#[test]
fn test_cem_f64() {
    xsref::test::<(f64, f64), _>("cem", "d_d_d-d_d", |x: &[f64]| xsf::cem(x[0], x[1], x[2]));
}

#[test]
fn test_sem_f64() {
    xsref::test::<(f64, f64), _>("sem", "d_d_d-d_d", |x: &[f64]| xsf::sem(x[0], x[1], x[2]));
}

#[test]
fn test_mcm1_f64() {
    xsref::test::<(f64, f64), _>("mcm1", "d_d_d-d_d", |x: &[f64]| xsf::mcm1(x[0], x[1], x[2]));
}

#[test]
fn test_msm1_f64() {
    xsref::test::<(f64, f64), _>("msm1", "d_d_d-d_d", |x: &[f64]| xsf::msm1(x[0], x[1], x[2]));
}

#[test]
fn test_mcm2_f64() {
    xsref::test::<(f64, f64), _>("mcm2", "d_d_d-d_d", |x: &[f64]| xsf::mcm2(x[0], x[1], x[2]));
}

#[test]
fn test_msm2_f64() {
    xsref::test::<(f64, f64), _>("msm2", "d_d_d-d_d", |x: &[f64]| xsf::msm2(x[0], x[1], x[2]));
}

#[test]
fn test_pbwa_f64() {
    xsref::test::<(f64, f64), _>("pbwa", "d_d-d_d", |x: &[f64]| xsf::pbwa(x[0], x[1]));
}

#[test]
fn test_pbdv_f64() {
    xsref::test::<(f64, f64), _>("pbdv", "d_d-d_d", |x: &[f64]| xsf::pbdv(x[0], x[1]));
}

#[test]
fn test_pbvv_f64() {
    xsref::test::<(f64, f64), _>("pbvv", "d_d-d_d", |x: &[f64]| xsf::pbvv(x[0], x[1]));
}

#[test]
fn test_sici_f64() {
    xsref::test::<(f64, f64), _>("sici", "d-d_d", |x: &[f64]| xsf::sici(x[0]));
}

#[test]
fn test_sici_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>), _>("sici", "cd-cd_cd", |x: &[f64]| {
        xsf::sici(c64(x[0], x[1]))
    });
}

#[test]
fn test_shichi_f64() {
    xsref::test::<(f64, f64), _>("shichi", "d-d_d", |x: &[f64]| xsf::shichi(x[0]));
}

#[test]
fn test_shichi_c64() {
    xsref::test::<(Complex<f64>, Complex<f64>), _>("shichi", "cd-cd_cd", |x: &[f64]| {
        xsf::shichi(c64(x[0], x[1]))
    });
}

#[test]
fn test_hyperu_f64() {
    xsref::test::<f64, _>("hyperu", "d_d_d-d", |x: &[f64]| {
        xsf::hyperu(x[0], x[1], x[2])
    });
}

#[test]
fn test_hyp1f1_c64() {
    xsref::test::<Complex<f64>, _>("hyp1f1", "d_d_cd-cd", |x: &[f64]| {
        xsf::hyp1f1(x[0], x[1], c64(x[2], x[3]))
    });
}

#[test]
fn test_pmv_f64() {
    xsref::test::<f64, _>("pmv", "d_d_d-d", |x: &[f64]| {
        xsf::pmv(x[0] as i64, x[1], x[2])
    });
}

// sph_bessel.h
// no xsref tables?
// xsref_test!(sph_bessel_j, "dd->d", "dD->D");
// xsref_test!(sph_bessel_y, "dd->d", "dD->D");
// xsref_test!(sph_bessel_i, "dd->d", "dD->D");
// xsref_test!(sph_bessel_k, "dd->d", "dD->D");
// xsref_test!(sph_bessel_j_jac, "dd->d", "dD->D");
// xsref_test!(sph_bessel_y_jac, "dd->d", "dD->D");
// xsref_test!(sph_bessel_i_jac, "dd->d", "dD->D");
// xsref_test!(sph_bessel_k_jac, "dd->d", "dD->D");

// sph_harm.h
// xsref_test!(sph_harm_y, "iidd->D");  // no xsref table?

#[test]
fn test_prolate_segv_f64() {
    xsref::test::<f64, _>("prolate_segv", "d_d_d-d", |x: &[f64]| {
        xsf::prolate_segv(x[0] as u64, x[1] as u64, x[2])
    });
}

#[test]
fn test_prolate_aswfa_nocv_f64() {
    xsref::test::<(f64, f64), _>("prolate_aswfa_nocv", "d_d_d_d-d_d", |x: &[f64]| {
        xsf::prolate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_oblate_aswfa_nocv_f64() {
    xsref::test::<(f64, f64), _>("oblate_aswfa_nocv", "d_d_d_d-d_d", |x: &[f64]| {
        xsf::oblate_aswfa_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_prolate_radial1_nocv_f64() {
    xsref::test::<(f64, f64), _>("prolate_radial1_nocv", "d_d_d_d-d_d", |x: &[f64]| {
        xsf::prolate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_oblate_radial1_nocv_f64() {
    xsref::test::<(f64, f64), _>("oblate_radial1_nocv", "d_d_d_d-d_d", |x: &[f64]| {
        xsf::oblate_radial1_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_prolate_radial2_nocv_f64() {
    xsref::test::<(f64, f64), _>("prolate_radial2_nocv", "d_d_d_d-d_d", |x: &[f64]| {
        xsf::prolate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_oblate_radial2_nocv_f64() {
    xsref::test::<(f64, f64), _>("oblate_radial2_nocv", "d_d_d_d-d_d", |x: &[f64]| {
        xsf::oblate_radial2_nocv(x[0] as u64, x[1] as u64, x[2], x[3])
    });
}

#[test]
fn test_prolate_aswfa_f64() {
    xsref::test::<(f64, f64), _>("prolate_aswfa", "d_d_d_d_d-d_d", |x: &[f64]| {
        xsf::prolate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_oblate_aswfa_f64() {
    xsref::test::<(f64, f64), _>("oblate_aswfa", "d_d_d_d_d-d_d", |x: &[f64]| {
        xsf::oblate_aswfa(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_prolate_radial1_f64() {
    xsref::test::<(f64, f64), _>("prolate_radial1", "d_d_d_d_d-d_d", |x: &[f64]| {
        xsf::prolate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_oblate_radial1_f64() {
    xsref::test::<(f64, f64), _>("oblate_radial1", "d_d_d_d_d-d_d", |x: &[f64]| {
        xsf::oblate_radial1(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_prolate_radial2_f64() {
    xsref::test::<(f64, f64), _>("prolate_radial2", "d_d_d_d_d-d_d", |x: &[f64]| {
        xsf::prolate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_oblate_radial2_f64() {
    xsref::test::<(f64, f64), _>("oblate_radial2", "d_d_d_d_d-d_d", |x: &[f64]| {
        xsf::oblate_radial2(x[0] as u64, x[1] as u64, x[2], x[3], x[4])
    });
}

#[test]
fn test_ndtr_f64() {
    xsref::test::<f64, _>("ndtr", "d-d", |x: &[f64]| xsf::ndtr(x[0]));
}

#[test]
fn test_ndtr_c64() {
    xsref::test::<Complex<f64>, _>("ndtr", "cd-cd", |x: &[f64]| xsf::ndtr(c64(x[0], x[1])));
}

#[test]
fn test_ndtri_f64() {
    xsref::test::<f64, _>("ndtri", "d-d", |x: &[f64]| xsf::ndtri(x[0]));
}

#[test]
fn test_kolmogorov_f64() {
    xsref::test::<f64, _>("kolmogorov", "d-d", |x: &[f64]| xsf::kolmogorov(x[0]));
}

#[test]
fn test_kolmogc_f64() {
    xsref::test::<f64, _>("kolmogc", "d-d", |x: &[f64]| xsf::kolmogc(x[0]));
}

#[test]
fn test_kolmogi_f64() {
    xsref::test::<f64, _>("kolmogi", "d-d", |x: &[f64]| xsf::kolmogi(x[0]));
}

#[test]
fn test_kolmogci_f64() {
    xsref::test::<f64, _>("kolmogci", "d-d", |x: &[f64]| xsf::kolmogci(x[0]));
}

#[test]
fn test_kolmogp_f64() {
    xsref::test::<f64, _>("kolmogp", "d-d", |x: &[f64]| xsf::kolmogp(x[0]));
}

#[test]
fn test_smirnov_f64() {
    xsref::test::<f64, _>("smirnov", "p_d-d", |x: &[f64]| {
        xsf::smirnov(x[0] as i32, x[1])
    });
}

#[test]
fn test_smirnovc_f64() {
    xsref::test::<f64, _>("smirnovc", "p_d-d", |x: &[f64]| {
        xsf::smirnovc(x[0] as i32, x[1])
    });
}

#[test]
fn test_smirnovi_f64() {
    xsref::test::<f64, _>("smirnovi", "p_d-d", |x: &[f64]| {
        xsf::smirnovi(x[0] as i32, x[1])
    });
}

#[test]
fn test_smirnovci_f64() {
    xsref::test::<f64, _>("smirnovci", "p_d-d", |x: &[f64]| {
        xsf::smirnovci(x[0] as i32, x[1])
    });
}

#[test]
fn test_smirnovp_f64() {
    xsref::test::<f64, _>("smirnovp", "p_d-d", |x: &[f64]| {
        xsf::smirnovp(x[0] as i32, x[1])
    });
}

#[test]
fn test_owens_t_f64() {
    xsref::test::<f64, _>("owens_t", "d_d-d", |x: &[f64]| xsf::owens_t(x[0], x[1]));
}

#[test]
fn test_chdtr_f64() {
    xsref::test::<f64, _>("chdtr", "d_d-d", |x: &[f64]| xsf::chdtr(x[0], x[1]));
}

#[test]
fn test_chdtrc_f64() {
    xsref::test::<f64, _>("chdtrc", "d_d-d", |x: &[f64]| xsf::chdtrc(x[0], x[1]));
}

#[test]
fn test_chdtri_f64() {
    xsref::test::<f64, _>("chdtri", "d_d-d", |x: &[f64]| xsf::chdtri(x[0], x[1]));
}

#[test]
fn test_fdtr_f64() {
    xsref::test::<f64, _>("fdtr", "d_d_d-d", |x: &[f64]| xsf::fdtr(x[0], x[1], x[2]));
}

#[test]
fn test_fdtrc_f64() {
    xsref::test::<f64, _>("fdtrc", "d_d_d-d", |x: &[f64]| xsf::fdtrc(x[0], x[1], x[2]));
}

#[test]
fn test_fdtri_f64() {
    xsref::test::<f64, _>("fdtri", "d_d_d-d", |x: &[f64]| xsf::fdtri(x[0], x[1], x[2]));
}

#[test]
fn test_gdtr_f64() {
    xsref::test::<f64, _>("gdtr", "d_d_d-d", |x: &[f64]| xsf::gdtr(x[0], x[1], x[2]));
}

#[test]
fn test_gdtrc_f64() {
    xsref::test::<f64, _>("gdtrc", "d_d_d-d", |x: &[f64]| xsf::gdtrc(x[0], x[1], x[2]));
}

#[test]
fn test_pdtr_f64() {
    xsref::test::<f64, _>("pdtr", "d_d-d", |x: &[f64]| xsf::pdtr(x[0], x[1]));
}

#[test]
fn test_pdtrc_f64() {
    xsref::test::<f64, _>("pdtrc", "d_d-d", |x: &[f64]| xsf::pdtrc(x[0], x[1]));
}

#[test]
fn test_pdtri_f64() {
    xsref::test::<f64, _>("pdtri", "p_d-d", |x: &[f64]| xsf::pdtri(x[0] as i32, x[1]));
}

#[test]
fn test_bdtr_f64() {
    xsref::test::<f64, _>("bdtr", "d_p_d-d", |x: &[f64]| {
        xsf::bdtr(x[0], x[1] as i32, x[2])
    });
}

#[test]
fn test_bdtrc_f64() {
    xsref::test::<f64, _>("bdtrc", "d_p_d-d", |x: &[f64]| {
        xsf::bdtrc(x[0], x[1] as i32, x[2])
    });
}

#[test]
fn test_bdtri_f64() {
    xsref::test::<f64, _>("bdtri", "d_p_d-d", |x: &[f64]| {
        xsf::bdtri(x[0], x[1] as i32, x[2])
    });
}

#[test]
fn test_nbdtr_f64() {
    xsref::test::<f64, _>("nbdtr", "p_p_d-d", |x: &[f64]| {
        xsf::nbdtr(x[0] as i32, x[1] as i32, x[2])
    });
}

#[test]
fn test_nbdtrc_f64() {
    xsref::test::<f64, _>("nbdtrc", "p_p_d-d", |x: &[f64]| {
        xsf::nbdtrc(x[0] as i32, x[1] as i32, x[2])
    });
}

#[test]
fn test_itstruve0_f64() {
    xsref::test::<f64, _>("itstruve0", "d-d", |x: &[f64]| xsf::itstruve0(x[0]));
}

#[test]
fn test_it2struve0_f64() {
    xsref::test::<f64, _>("it2struve0", "d-d", |x: &[f64]| xsf::it2struve0(x[0]));
}

#[test]
fn test_itmodstruve0_f64() {
    xsref::test::<f64, _>("itmodstruve0", "d-d", |x: &[f64]| xsf::itmodstruve0(x[0]));
}

#[test]
fn test_struve_h_f64() {
    xsref::test::<f64, _>("struve_h", "d_d-d", |x: &[f64]| xsf::struve_h(x[0], x[1]));
}

#[test]
fn test_struve_l_f64() {
    xsref::test::<f64, _>("struve_l", "d_d-d", |x: &[f64]| xsf::struve_l(x[0], x[1]));
}

#[test]
fn test_sinpi_f64() {
    xsref::test::<f64, _>("sinpi", "d-d", |x: &[f64]| xsf::sinpi(x[0]));
}

#[test]
fn test_sinpi_c64() {
    xsref::test::<Complex<f64>, _>("sinpi", "cd-cd", |x: &[f64]| xsf::sinpi(c64(x[0], x[1])));
}

#[test]
fn test_cospi_f64() {
    xsref::test::<f64, _>("cospi", "d-d", |x: &[f64]| xsf::cospi(x[0]));
}

#[test]
fn test_cospi_c64() {
    xsref::test::<Complex<f64>, _>("cospi", "cd-cd", |x: &[f64]| xsf::cospi(c64(x[0], x[1])));
}

#[test]
fn test_sindg_f64() {
    xsref::test::<f64, _>("sindg", "d-d", |x: &[f64]| xsf::sindg(x[0]));
}

#[test]
fn test_cosdg_f64() {
    xsref::test::<f64, _>("cosdg", "d-d", |x: &[f64]| xsf::cosdg(x[0]));
}

#[test]
fn test_tandg_f64() {
    xsref::test::<f64, _>("tandg", "d-d", |x: &[f64]| xsf::tandg(x[0]));
}

#[test]
fn test_cotdg_f64() {
    xsref::test::<f64, _>("cotdg", "d-d", |x: &[f64]| xsf::cotdg(x[0]));
}

#[test]
fn test_cosm1_f64() {
    xsref::test::<f64, _>("cosm1", "d-d", |x: &[f64]| xsf::cosm1(x[0]));
}

#[test]
fn test_radian_f64() {
    xsref::test::<f64, _>("radian", "d_d_d-d", |x: &[f64]| {
        xsf::radian(x[0], x[1], x[2])
    });
}

#[test]
fn test_wright_bessel_f64() {
    xsref::test::<f64, _>("wright_bessel", "d_d_d-d", |x: &[f64]| {
        xsf::wright_bessel(x[0], x[1], x[2])
    });
}

#[test]
fn test_log_wright_bessel_f64() {
    xsref::test::<f64, _>("log_wright_bessel", "d_d_d-d", |x: &[f64]| {
        xsf::log_wright_bessel(x[0], x[1], x[2])
    });
}

#[test]
fn test_riemann_zeta_f64() {
    xsref::test::<f64, _>("riemann_zeta", "d-d", |x: &[f64]| xsf::riemann_zeta(x[0]));
}

#[test]
fn test_riemann_zeta_c64() {
    xsref::test::<Complex<f64>, _>("riemann_zeta", "cd-cd", |x: &[f64]| {
        xsf::riemann_zeta(c64(x[0], x[1]))
    });
}

#[test]
fn test_zeta_f64() {
    xsref::test::<f64, _>("zeta", "d_d-d", |x: &[f64]| xsf::zeta(x[0], x[1]));
}

#[test]
fn test_zetac_f64() {
    xsref::test::<f64, _>("zetac", "d-d", |x: &[f64]| xsf::zetac(x[0]));
}
