# Docs

| file                         | contents                                                                              |
| ---------------------------- | ------------------------------------------------------------------------------------- |
| [`accuracy.md`](accuracy.md) | relative errors versus mpmath, next to xsf (or SciPy or NumPy) and the Rust std       |
| [`hyppfq.md`](hyppfq.md)     | where `hyppfq` returns `None`, and its errors elsewhere, for $`p \le 4`$, $`q \le 3`$ |
| [`scipy.md`](scipy.md)       | the `scipy.special` functions that SciPy computes partly or entirely without xsf      |

## Tools

| `just` recipe | does                                                                                                     |
| ------------- | -------------------------------------------------------------------------------------------------------- |
| `accuracy`    | regenerates [`accuracy.md`](accuracy.md)                                                                 |
| `pfq-map`     | regenerates the maps of [`hyppfq.md`](hyppfq.md)                                                         |
| `fuzz-diff`   | finds random inputs where the working tree (or `--head`) is less accurate than `--base` (Arb checks)     |
| `domain-diff` | counts per domain the inputs where the working tree (or `--head`) is more or less accurate than `--base` |

Usage: `just fuzz-diff [--base REV] [--head REV] [--only NAME ...] [--cases N]` and
`just domain-diff [--base REV] [--head REV] [--only NAME ...] [--points N]` (exit 1 on a regression,
2 if a comparison or reference is missing); see
[`fuzz_diff.py`](../tools/accuracy/fuzz_diff.py) and [`domain_diff.py`](../tools/accuracy/domain_diff.py).
