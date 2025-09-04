//! Automatic testing of the `scipy/xsf` rust bindings using the `scipy/xsref` tables.
mod xsref {
    use std::env;
    use std::fs::File;
    use std::path::PathBuf;

    use arrow::array::{
        Array, Float32Array as ArrayF32, Float64Array as ArrayF64, Int32Array as ArrayI32,
        Int64Array as ArrayI64,
    };
    use num_complex::Complex;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

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
            rows.into_iter()
                .map(|row| Complex::new(row[0], row[1]))
                .collect()
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

/// Helper macro to generate the test body
macro_rules! _test {
    ($test_name:ident, $f:ident, $sig:literal, $test_fn:expr, f64) => {
        paste::paste! {
            #[test]
            fn $test_name() {
                xsref::test::<f64, _>(stringify!($f), $sig, $test_fn);
            }
        }
    };
    ($test_name:ident, $f:ident, $sig:literal, $test_fn:expr, Complex<f64>) => {
        paste::paste! {
            #[test]
            fn $test_name() {
                xsref::test::<num_complex::Complex<f64>, _>(stringify!($f), $sig, $test_fn);
            }
        }
    };
}

/// Generate a test function for xsf functions
macro_rules! xsref_test {
    // Single signature (backward compatibility)
    ($f:ident, $sig:tt) => {
        xsref_test!(@single $f, $sig);
    };

    // Multiple signatures (variadic)
    ($f:ident, $($sig:tt),+ $(,)?) => {
        $(
            xsref_test!(@single $f, $sig);
        )+
    };

    // Internal helper for single signature processing
    (@single $f:ident, "d->d") => {
        paste::paste! {
            _test!([<test_ $f _ d>], $f, "d-d", |x: &[f64]| xsf::$f(x[0]), f64);
        }
    };
    (@single $f:ident, "dd->d") => {
        paste::paste! {
            _test!([<test_ $f _ d>], $f, "d_d-d", |x: &[f64]| xsf::$f(x[0], x[1]), f64);
        }
    };
    (@single $f:ident, "id->d") => {
        paste::paste! {
            _test!([<test_ $f _ d>], $f, "p_d-d", |x: &[f64]| xsf::$f(x[0] as i32, x[1]), f64);
        }
    };
    (@single $f:ident, "ddd->d") => {
        paste::paste! {
            _test!([<test_ $f _ d>], $f, "d_d_d-d", |x: &[f64]| xsf::$f(x[0], x[1], x[2]), f64);
        }
    };
    (@single $f:ident, "QQd->d") => {
        paste::paste! {
            _test!(
                [<test_ $f _ d>],
                $f, "d_d_d-d",
                |x: &[f64]| xsf::$f(x[0] as u64, x[1] as u64, x[2]),
                f64
            );
        }
    };
    (@single $f:ident, "qdd->d") => {
        paste::paste! {
            _test!(
                [<test_ $f _ d>],
                $f, "d_d_d-d",
                |x: &[f64]| xsf::$f(x[0] as i64, x[1], x[2]),
                f64
            );
        }
    };
    (@single $f:ident, "did->d") => {
        paste::paste! {
            _test!(
                [<test_ $f _ d>],
                $f,
                "d_p_d-d",
                |x: &[f64]| xsf::$f(x[0], x[1] as i32, x[2]),
                f64
            );
        }
    };
    (@single $f:ident, "iid->d") => {
        paste::paste! {
            _test!(
                [<test_ $f _ d>],
                $f,
                "p_p_d-d",
                |x: &[f64]| xsf::$f(x[0] as i32, x[1] as i32, x[2]),
                f64
            );
        }
    };
    (@single $f:ident, "dddd->d") => {
        paste::paste! {
            _test!(
                [<test_ $f _ d>],
                $f,
                "d_d_d_d-d",
                |x: &[f64]| xsf::$f(x[0], x[1], x[2], x[3]),
                f64
            );
        }
    };
    (@single $f:ident, "D->D") => {
        paste::paste! {
            _test!(
                [<test_ $f _ cd>],
                $f,
                "cd-cd",
                |x: &[f64]| xsf::$f(num_complex::Complex::new(x[0], x[1])),
                Complex<f64>
            );
        }
    };
    (@single $f:ident, "dD->D") => {
        paste::paste! {
            _test!(
                [<test_ $f _ cd>],
                $f,
                "d_cd-cd",
                |x: &[f64]| xsf::$f(x[0], num_complex::Complex::new(x[1], x[2])),
                Complex<f64>
            );
        }
    };
    (@single $f:ident, "Dd->D") => {
        paste::paste! {
            _test!(
                [<test_ $f _ cd>],
                $f,
                "cd_d-cd",
                |x: &[f64]| xsf::$f(num_complex::Complex::new(x[0], x[1]), x[2]),
                Complex<f64>
            );
        }
    };
    (@single $f:ident, "DD->D") => {
        paste::paste! {
            _test!(
                [<test_ $f _ cd>],
                $f,
                "cd_cd-cd",
                |x: &[f64]| xsf::$f(
                    num_complex::Complex::new(x[0], x[1]),
                    num_complex::Complex::new(x[2], x[3]),
                ),
                Complex<f64>
            );
        }
    };
    (@single $f:ident, "ddD->D") => {
        paste::paste! {
            _test!(
                [<test_ $f _ cd>],
                $f,
                "d_d_cd-cd",
                |x: &[f64]| xsf::$f(x[0], x[1], num_complex::Complex::new(x[2], x[3])),
                Complex<f64>
            );
        }
    };
    (@single $f:ident, "dddD->D") => {
        paste::paste! {
            _test!(
                [<test_ $f _ cd>],
                $f,
                "d_d_d_cd-cd",
                |x: &[f64]| xsf::$f(x[0], x[1], x[2], num_complex::Complex::new(x[3], x[4])),
                Complex<f64>
            );
        }
    };
    (@single $f:ident, "Dld->D") => {
        paste::paste! {
            _test!(
                [<test_ $f _ cd>],
                $f,
                "cd_p_d-cd",
                |x: &[f64]| xsf::$f(
                    num_complex::Complex::new(x[0], x[1]),
                    x[2] as std::os::raw::c_long,
                    x[3],
                ),
                Complex<f64>
            );
        }
    };
}

// alg.h
xsref_test!(cbrt, "d->d");

// bessel.h
xsref_test!(cyl_bessel_j0, "d->d");
xsref_test!(cyl_bessel_j1, "d->d");
xsref_test!(cyl_bessel_y0, "d->d");
xsref_test!(cyl_bessel_y1, "d->d");
xsref_test!(cyl_bessel_i0, "d->d");
xsref_test!(cyl_bessel_i0e, "d->d");
xsref_test!(cyl_bessel_i1, "d->d");
xsref_test!(cyl_bessel_i1e, "d->d");
xsref_test!(cyl_bessel_k0, "d->d");
xsref_test!(cyl_bessel_k0e, "d->d");
xsref_test!(cyl_bessel_k1, "d->d");
xsref_test!(cyl_bessel_k1e, "d->d");
xsref_test!(cyl_bessel_j, "dd->d", "dD->D");
xsref_test!(cyl_bessel_je, "dd->d", "dD->D");
xsref_test!(cyl_bessel_y, "dd->d", "dD->D");
xsref_test!(cyl_bessel_ye, "dd->d", "dD->D");
xsref_test!(cyl_bessel_i, "dd->d", "dD->D");
xsref_test!(cyl_bessel_ie, "dd->d", "dD->D");
xsref_test!(cyl_bessel_k, "dd->d", "dD->D");
xsref_test!(cyl_bessel_ke, "dd->d", "dD->D");
xsref_test!(cyl_hankel_1, "dD->D");
xsref_test!(cyl_hankel_1e, "dD->D");
xsref_test!(cyl_hankel_2, "dD->D");
xsref_test!(cyl_hankel_2e, "dD->D");
xsref_test!(besselpoly, "ddd->d");

// beta.h
xsref_test!(beta, "dd->d");
xsref_test!(betaln, "dd->d");

// binom.h
xsref_test!(binom, "dd->d");

// cdflib.h
xsref_test!(gdtrib, "ddd->d");

// digamma.h
xsref_test!(digamma, "d->d", "D->D");

// ellip.h
xsref_test!(ellipk, "d->d");
xsref_test!(ellipkm1, "d->d");
xsref_test!(ellipkinc, "dd->d");
xsref_test!(ellipe, "d->d");
xsref_test!(ellipeinc, "dd->d");

// erf.h
xsref_test!(erf, "d->d", "D->D");
xsref_test!(erfc, "d->d", "D->D");
xsref_test!(erfcx, "d->d", "D->D");
xsref_test!(erfi, "d->d", "D->D");
xsref_test!(dawsn, "d->d", "D->D");
xsref_test!(wofz, "D->D");
xsref_test!(voigt_profile, "ddd->d");

// exp.h
xsref_test!(expm1, "d->d", "D->D");
xsref_test!(exp2, "d->d");
xsref_test!(exp10, "d->d");

// expint.h
xsref_test!(expi, "d->d", "D->D");
xsref_test!(exp1, "d->d", "D->D");
xsref_test!(scaled_exp1, "d->d");

// gamma.h
xsref_test!(gamma, "d->d", "D->D");
xsref_test!(gammainc, "dd->d");
xsref_test!(gammaincc, "dd->d");
xsref_test!(gammaincinv, "dd->d");
xsref_test!(gammainccinv, "dd->d");
xsref_test!(gammaln, "d->d");
xsref_test!(gammasgn, "d->d");

// hyp2f1.h
xsref_test!(hyp2f1, "dddd->d", "dddD->D");

// iv_ratio.h
xsref_test!(iv_ratio, "dd->d");
xsref_test!(iv_ratio_c, "dd->d");

// kelvin.h
xsref_test!(ber, "d->d");
xsref_test!(bei, "d->d");
xsref_test!(ker, "d->d");
xsref_test!(kei, "d->d");
xsref_test!(berp, "d->d");
xsref_test!(beip, "d->d");
xsref_test!(kerp, "d->d");
xsref_test!(keip, "d->d");

// lambertw.h
xsref_test!(lambertw, "Dld->D");

// legendre.h
// no xsref tables?

// log_exp.h
xsref_test!(expit, "d->d");
xsref_test!(exprel, "d->d");
xsref_test!(logit, "d->d");
xsref_test!(log_expit, "d->d");
// xsref_test!(log1mexp, "d->d");  // no xsref table

// log.h
xsref_test!(log1p, "d->d", "D->D");
xsref_test!(log1pmx, "d->d");
xsref_test!(xlogy, "dd->d", "DD->D");
xsref_test!(xlog1py, "dd->d", "DD->D");

// loggamma.h
xsref_test!(loggamma, "d->d", "D->D");
xsref_test!(rgamma, "d->d", "D->D");

// mathieu.h
xsref_test!(cem_cva, "dd->d");
xsref_test!(sem_cva, "dd->d");

// specfun.h
xsref_test!(hyperu, "ddd->d"); // alias of `hypu`
xsref_test!(hyp1f1, "ddD->D"); // no table for ddd->d
xsref_test!(pmv, "qdd->d");

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

// sphd_wave.h
xsref_test!(prolate_segv, "QQd->d");
// xsref_test!(oblate_segv, "QQd->d");  // missing xsref table?

// stats.h
xsref_test!(ndtr, "d->d", "D->D");
xsref_test!(ndtri, "d->d");
// xsref_test!(log_ndtr, "d->d", "D->D");  // no xsref table
xsref_test!(kolmogorov, "d->d");
xsref_test!(kolmogc, "d->d");
xsref_test!(kolmogi, "d->d");
xsref_test!(kolmogci, "d->d");
xsref_test!(kolmogp, "d->d");
xsref_test!(smirnov, "id->d");
xsref_test!(smirnovc, "id->d");
xsref_test!(smirnovi, "id->d");
xsref_test!(smirnovci, "id->d");
xsref_test!(smirnovp, "id->d");
// xsref_test!(tukeylambdacdf, "dd->d");  // no xsref table
xsref_test!(owens_t, "dd->d");
xsref_test!(chdtr, "dd->d");
xsref_test!(chdtrc, "dd->d");
xsref_test!(chdtri, "dd->d");
xsref_test!(fdtr, "ddd->d");
xsref_test!(fdtrc, "ddd->d");
xsref_test!(fdtri, "ddd->d");
xsref_test!(gdtr, "ddd->d");
xsref_test!(gdtrc, "ddd->d");
xsref_test!(pdtr, "dd->d");
xsref_test!(pdtrc, "dd->d");
xsref_test!(pdtri, "id->d");
xsref_test!(bdtr, "did->d");
xsref_test!(bdtrc, "did->d");
xsref_test!(bdtri, "did->d");
xsref_test!(nbdtr, "iid->d");
xsref_test!(nbdtrc, "iid->d");
// xsref_test!(nbdtri, "iid->d");  // no xsref table

// struve.h
xsref_test!(itstruve0, "d->d");
xsref_test!(it2struve0, "d->d");
xsref_test!(itmodstruve0, "d->d");
xsref_test!(struve_h, "dd->d");
xsref_test!(struve_l, "dd->d");

// trig.h
xsref_test!(sinpi, "d->d", "D->D");
xsref_test!(cospi, "d->d", "D->D");
xsref_test!(sindg, "d->d");
xsref_test!(cosdg, "d->d");
xsref_test!(tandg, "d->d");
xsref_test!(cotdg, "d->d");
xsref_test!(cosm1, "d->d");
xsref_test!(radian, "ddd->d");

// wright_bessel.h
xsref_test!(wright_bessel, "ddd->d");
xsref_test!(log_wright_bessel, "ddd->d");

// zeta.h
xsref_test!(riemann_zeta, "d->d", "D->D");
xsref_test!(zeta, "dd->d"); // no complex xsref table
xsref_test!(zetac, "d->d");
