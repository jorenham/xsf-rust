# xsref reference tables

The reference tables of the tests (read by `testing/src/xsref.rs`), converted from the parquet
tables in `src/xsref/tables/scipy_special_tests` of
[scipy/xsref](https://github.com/scipy/xsref) at commit
[`e55f882`](https://github.com/scipy/xsref/tree/e55f8825a74246703417e2d5ce17e22b0e03d310)
(BSD-3-Clause, see [`LICENSE`](LICENSE)). Only the tables that a test uses are included.

They were generated with [`tools/tables/convert.py`](../../../tools/tables/convert.py):

```sh
git clone https://github.com/scipy/xsref /tmp/xsref
git -C /tmp/xsref checkout e55f8825a74246703417e2d5ce17e22b0e03d310
uv run tools/tables/convert.py /tmp/xsref tests/data/xsref
```

[`MANIFEST.tsv`](MANIFEST.tsv) lists the rows, the number of rows in xsref, the dimensions, and
the SHA-256 of each file.

## Contents

Each `<table>/<sig>.bin` holds the rows of xsref's `In_<sig>`, `Out_<sig>`, and
`Err_<sig>_other` tables, in their original order. Every column is cast to `f64` on its own
(`f32` values are widened, nulls become NaN), and the column order is kept. The `fallback` column
of `Out_<sig>` (the reference is xsf's own value) is the row flag, not an output.

Only `cyl_bessel_i` is subsampled, to 20,000 rows per signature: 16,000 rows where all outputs
are finite, 2,000 with a NaN output, and 2,000 with an infinite output (a class with fewer rows
gives its remainder to the finite one), each drawn uniformly without replacement with
`random.Random(0)`, in that order. On top of these, every row with a half-integer order `v`, or
one in `BESSEL_I_OVERFLOWS` (`testing/src/xsref.rs`), is kept. `--no-cap` keeps all rows.

## Format (version 1)

All numbers are little-endian. A complex output is two `f64`s (real part first) with one
tolerance.

| Bytes                             | Type  | Content                                       |
| --------------------------------- | ----- | --------------------------------------------- |
| 4                                 |       | magic `XSFT`                                  |
| 4                                 | `u32` | version: 1                                    |
| 4                                 | `u32` | `n_rows`                                      |
| 2                                 | `u16` | `n_in`: input `f64`s per row                  |
| 2                                 | `u16` | `n_out`: output `f64`s per row                |
| 2                                 | `u16` | `n_tol`: logical outputs (tolerances) per row |
| 2                                 | `u16` | flags: 0                                      |
| `8 n_rows (n_in + n_out + n_tol)` | `f64` | inputs, outputs, tolerances (row-major)       |
| `n_rows`                          | `u8`  | row flags: bit 0 is xsref's `fallback`        |

So a file is exactly `20 + n_rows (8 (n_in + n_out + n_tol) + 1)` bytes long. The header flags
are reserved for optional sections: bit 0 for a status byte per logical output after the row
flags, and bit 1 for an `f64` error scale per logical output after the tolerances. Readers reject
other versions, nonzero header flags, and unknown row flags.
