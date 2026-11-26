# Docs

| file                         | contents                                                                              |
| ---------------------------- | ------------------------------------------------------------------------------------- |
| [`accuracy.md`](accuracy.md) | relative errors versus mpmath, next to xsf (or SciPy or NumPy) and the Rust std       |
| [`hyppfq.md`](hyppfq.md)     | where `hyppfq` returns `None`, and its errors elsewhere, for $`p \le 4`$, $`q \le 3`$ |
| [`scipy.md`](scipy.md)       | the `scipy.special` functions that SciPy computes partly or entirely without xsf      |

## Tools

| `just` recipe | does                                                                                                 |
| ------------- | ---------------------------------------------------------------------------------------------------- |
| `accuracy`    | regenerates [`accuracy.md`](accuracy.md)                                                             |
| `pfq-map`     | regenerates the maps of [`hyppfq.md`](hyppfq.md)                                                     |
| `diff-fuzz`   | finds random inputs where the working tree (or `--head`) is less accurate than `--base` (Arb checks) |

Usage: `just diff-fuzz [--base REV] [--head REV] [--only NAME ...] [--cases N]`; see
[`diff_fuzz.py`](../tools/accuracy/diff_fuzz.py).
