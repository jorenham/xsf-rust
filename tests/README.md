# Testing

Tests for xsf special functions with validation against SciPy reference data.

## Running Tests

```bash
# Run all tests
cargo test

# With reference data (optional)
./tests/setup_test_data.sh
cargo test
```

## Test Coverage

- **gamma** - Gamma function
- **cyl_bessel_j**, **cyl_bessel_j0**, **cyl_bessel_j1** - Bessel functions
- **hyp2f1** - Hypergeometric function

## Reference Data

Tests automatically use SciPy's reference test data when available in `xsref/`. This provides comprehensive validation against high-precision reference values computed with arbitrary precision arithmetic.

If reference data is not available, tests fall back to basic mathematical validation using known function properties.
