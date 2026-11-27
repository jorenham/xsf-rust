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
| `domain-diff` | counts per domain the better and worse inputs, and the time ratio, of the working tree vs `--base`       |

Usage: `just fuzz-diff [--base REV] [--head REV] [--only NAME ...] [--cases N]` and
`just domain-diff [--base REV] [--head REV] [--only NAME ...] [--points N]`; see
[`fuzz_diff.py`](../tools/accuracy/fuzz_diff.py) and [`domain_diff.py`](../tools/accuracy/domain_diff.py).

| `domain-diff` exit | meaning (checked in this order)                                                       |
| ------------------ | ------------------------------------------------------------------------------------- |
| 1 `FAIL`           | a regression                                                                          |
| 3 `SLOW`           | a time ratio above 1.25 (`-` if not measured: a panic, timeout, or old base)          |
| 2 `INCONCLUSIVE`   | a comparison or reference is missing                                                  |

CI fails on `SLOW` unless the PR has the label `allow-slowdown`; add it, then re-run the job.

CI runs [`check_loosened_tests.py`](../tools/check_loosened_tests.py), which warns about loosened
tolerances, new property-test exclusions, and removed or ignored tests. It fails the PR unless it has
the label `allow-loosened-tests`, which re-runs the check when added or removed.
