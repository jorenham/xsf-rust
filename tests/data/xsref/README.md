# xsref reference tables

The reference tables of the tests (read by `testing/src/xsref.rs`), converted from the parquet
tables in `src/xsref/tables/scipy_special_tests` of
[scipy/xsref](https://github.com/scipy/xsref) at commit
[`9890a90`](https://github.com/scipy/xsref/tree/9890a90db2fb85a50f62d84b02e611526c06e1e7)
(BSD-3-Clause, see [`LICENSE`](LICENSE)). Only the tables that a test uses are included.

They were generated with [`tools/tables/convert.py`](../../../tools/tables/convert.py):

```sh
git clone https://github.com/scipy/xsref /tmp/xsref
git -C /tmp/xsref checkout 9890a90db2fb85a50f62d84b02e611526c06e1e7
uv run --project tools tools/tables/convert.py /tmp/xsref tests/data/xsref
```

[`MANIFEST.tsv`](MANIFEST.tsv) lists the rows, the number of rows in xsref, the dimensions, and
the SHA-256 of each file.

## Contents

Each `<table>/<sig>.bin` holds the rows of xsref's `In_<sig>`, `Out_<sig>`, and
`Err_<sig>_other` tables, in their original order. Every column is cast to `f64` on its own
(`f32` values are widened, nulls become NaN), and the column order is kept. The `fallback` column
of `Out_<sig>` (the reference is xsf's own value) is the row flag, not an output.

`cyl_bessel_i` is subsampled to 20,000 rows per signature: 16,000 rows where all outputs
are finite, 2,000 with a NaN output, and 2,000 with an infinite output (a class with fewer rows
gives its remainder to the finite one), each drawn uniformly without replacement with
`random.Random(0)`, in that order. On top of these, every row with a half-integer order `v`, or
one in `BESSEL_I_OVERFLOWS` (`tools/tables/convert.py`), is kept.

Of the 173,167 rows of the real `airy` table with -1e4 < x < -1e3 (near the zeros of Ai, Ai', Bi,
and Bi'), every 8th is kept (from the first one, in xsref's order).

The `spence` table duplicates the real source output so the real overload and the real-axis
projection of the complex overload can have separate tolerances on the same inputs.

`--no-cap` keeps all rows.

## Generated tables

[`tools/tables/generate.py`](../../../tools/tables/generate.py) replaces xsref's references and
tolerances by our own, for the tables of
[`tools/tables/mapping.py`](../../../tools/tables/mapping.py):

```sh
uv run --project tools tools/tables/convert.py /tmp/xsref /tmp/baseline
uv run --project tools tools/tables/generate.py /tmp/baseline
```

The inputs, the row order, and the row flags stay xsref's. Each output gets a status (see below):

- computed: an mpmath reference of `tools/accuracy/specs.py` with 50 and with 100 digits, which
  must agree to 1e-30. Arb must confirm it where the Arb version of the mapping applies; an output
  without one needs an identity. Its tolerance is the error of xsf (C++), or at most the largest
  error of xsf-rust where the mapping marks xsf as wrong, and at least the floor of the mapping
  where xsf-rust uses another method than xsf, or where the error depends on the platform's libm.
  For an `f32` output, it is at least the `f32` epsilon.
- convention: a rule of the mapping, which xsf-rust must follow already, and which must pass
  xsref's test, unless the mapping marks xsf as wrong there.
- gap: a rule of the mapping with `gap`, where xsf-rust returns something else than the finite
  true value (e.g. NaN where `ellipkinc` is real). The reference is computed as above, but
  neither xsf-rust nor xsf need be near it. The tolerance is `max(2^-26, floor)` for `f64` and
  `max(2^-12, floor)` for `f32`, and the test must fail there, so that it notices when the gap
  closes. The generator errors if the rule's value would pass.
- legacy: where our reference timed out.

## Statuses and scales

Each logical output of each row has a status:

| Status       | Reference                                      | Checked before the tolerance       |
| ------------ | ---------------------------------------------- | ---------------------------------- |
| 0 legacy     | xsref's, with xsref's tolerance                |                                    |
| 1 computed   | the `f64` (or `f32`) nearest to the true value | NaN and the infinities             |
| 2 convention | a rule, e.g. NaN outside of the domain         | NaN, the infinities, and the zeros |
| 3 gap        | the `f64` (or `f32`) nearest to the true value | NaN, the infinities, and the zeros |

A gap is checked like a computed output, but must fail: the test counts it as a known gap, and
fails if it passes (the gap is closed, so replace the rule by a computed reference, with an
`Override` where xsf still differs, and regenerate). Keep gaps to classes (e.g. NaN for a finite
value, or zero for a nonzero one), which do not flip on a rounding error. The zeros are a class
of their own, as for a convention.

An output passes if its error is at most `4 max(eps, tol)`, with `eps` the `f64::EPSILON`. The error
is xsf's `extended_relative_error`. When one component is infinite and the other is finite,
the error also includes the finite component's scalar relative error if larger. Use xsf's
`extended_absolute_error` over the scale (the local
amplitude of an oscillating function, the size of the terms that cancel in a sum near its zeros, 1
for an absolute error, the smallest normal `f64` (or `f32`) for a subnormal reference without one,
or 0 for none) if that is smaller.
[`testing/metric_vectors.tsv`](../../../testing/metric_vectors.tsv) has test vectors.

## Format (version 1)

All numbers are little-endian. A complex output is two `f64`s (real part first) with one
tolerance, scale, and status. Header flag bit 0 announces the statuses and bit 1 the scales
(`s_0` and `s_1` below are 1 with that flag, and 0 without).

| Bytes                                       | Type  | Content                                         |
| ------------------------------------------- | ----- | ----------------------------------------------- |
| 4                                           |       | magic `XSFT`                                    |
| 4                                           | `u32` | version: 1                                      |
| 4                                           | `u32` | `n_rows`                                        |
| 2                                           | `u16` | `n_in`: input `f64`s per row                    |
| 2                                           | `u16` | `n_out`: output `f64`s per row                  |
| 2                                           | `u16` | `n_tol`: logical outputs (tolerances) per row   |
| 2                                           | `u16` | flags: bit 0 and bit 1                          |
| `8 n_rows (n_in + n_out + n_tol (1 + s_1))` | `f64` | inputs, outputs, tolerances, scales (row-major) |
| `n_rows`                                    | `u8`  | row flags: bit 0 is xsref's `fallback`          |
| `n_rows n_tol s_0`                          | `u8`  | statuses (row-major)                            |

So a file is exactly `20 + n_rows (8 (n_in + n_out + n_tol (1 + s_1)) + 1 + n_tol s_0)` bytes
long. Without the statuses, every output is a legacy one. Readers reject other versions, and
unknown flags or statuses.
