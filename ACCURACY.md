# Accuracy

This table compares the accuracy of the f64 functions of xsf-rust 0.6.0+0.2.2, of the C++
[xsf](https://github.com/scipy/xsf) library that it is a translation of, and, where
available, of the Rust standard library. The xsf version is commit
[`b6cb425`](https://github.com/scipy/xsf/tree/b6cb42595c449c93989b9ad5eecfe797a67ddacb).

**Method.** For each function and domain, N pseudo-random f64 inputs are drawn (with a fixed
seed) uniformly, or log-uniformly where marked with "(log)"; "±" means a random sign. Each
output y is compared to a reference value f, computed with [mpmath](https://mpmath.org/)
1.4.1 at 50 or more significant digits for the exact f64 input; where mpmath is known to
be inaccurate, the reference uses an alternative formula instead (see
`tools/accuracy/specs.py`). The error is the relative error |y − f| / |f|, except where
noted below a table. A result counts as exact if f overflows and y is an infinity of the same
sign, or if f underflows (|f| < 2.2e-308) and so does y. "peak" is the largest error, and
"rms" the root mean square of the finite errors; "inf (k)" means that k results have an
infinite error, e.g. NaN, or a spurious overflow. "bit-identical" is the fraction of the
inputs for which xsf-rust and xsf return the same f64 (rounded down, and counting any two NaNs
as equal). For reference, the machine epsilon is 2.2e-16.

Built on x86_64 Linux with:

- xsf-rust: `rustc 1.99.1 (58cf8f9c3 2026-10-13)`, release profile
- xsf: `g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0`, `-O2 -std=c++17`

To regenerate this file, run `just accuracy`. This needs `cargo`, `g++`, `uv`, and an
internet connection to download xsf. The references are cached in `tools/accuracy/.cache/`.

## Gamma and related functions

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>gamma</code></td><td>x ∈ [−30, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-15</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>1.1e-15</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>gamma</code></td><td>x ∈ [−170, 171.6]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.4e-16</samp></td><td align="right" nowrap><samp>1.8e-16</samp></td><td align="right" nowrap><samp>8.4e-16</samp></td><td align="right" nowrap><samp>1.8e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>gamma</code></td><td>x ∈ [10⁻³⁰⁰, 171] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>4.8e-16</samp></td><td align="right" nowrap><samp>5.7e-17</samp></td><td align="right" nowrap><samp>4.8e-16</samp></td><td align="right" nowrap><samp>5.7e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>gammaln</code></td><td>x ∈ [−30, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.3e-13</samp></td><td align="right" nowrap><samp>3.6e-15</samp></td><td align="right" nowrap><samp>3.3e-13</samp></td><td align="right" nowrap><samp>3.6e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>gammaln</code></td><td>x ∈ [10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.8e-15</samp></td><td align="right" nowrap><samp>7.3e-17</samp></td><td align="right" nowrap><samp>1.8e-15</samp></td><td align="right" nowrap><samp>7.3e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>rgamma</code></td><td>x ∈ [−30, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>9.5e-16</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>9.5e-16</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>rgamma</code></td><td>x ∈ [−200, 200]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>inf (11)</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap><samp>inf (11)</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>digamma</code></td><td>x ∈ [−30, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.4e-13</samp></td><td align="right" nowrap><samp>2.8e-15</samp></td><td align="right" nowrap><samp>8.2e-12</samp></td><td align="right" nowrap><samp>8.3e-14</samp></td><td align="right" nowrap>78.6%</td></tr>
<tr><td><code>digamma</code></td><td>x ∈ [10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.0e-16</samp></td><td align="right" nowrap><samp>4.9e-17</samp></td><td align="right" nowrap><samp>2.0e-16</samp></td><td align="right" nowrap><samp>4.9e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>digammainv</code></td><td>y ∈ [−10, 10]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.4e-15</samp></td><td align="right" nowrap><samp>3.5e-16</samp></td><td align="right" nowrap><samp>inf (10)</samp></td><td align="right" nowrap><samp>3.5e-16</samp></td><td align="right" nowrap>99.9%</td></tr>
<tr><td><code>digammainv</code></td><td>y ∈ [−1000, 700]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.5e-15</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap><samp>inf (1)</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap>99.9%</td></tr>
<tr><td><code>beta</code></td><td>a, b ∈ [0, 20]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.4e-14</samp></td><td align="right" nowrap><samp>3.7e-15</samp></td><td align="right" nowrap><samp>1.4e-14</samp></td><td align="right" nowrap><samp>3.7e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>beta</code></td><td>a, b ∈ [10⁻³⁰⁰, 170] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.1e-14</samp></td><td align="right" nowrap><samp>7.7e-16</samp></td><td align="right" nowrap><samp>6.1e-14</samp></td><td align="right" nowrap><samp>7.7e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>beta</code></td><td>a ∈ ±⁠[10⁻³, 10⁶] (log), b ∈ [10⁻³, 10⁶] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.0e-08</samp></td><td align="right" nowrap><samp>1.6e-10</samp></td><td align="right" nowrap><samp>1.0e+00</samp></td><td align="right" nowrap><samp>1.0e-02</samp></td><td align="right" nowrap>53.7%</td></tr>
<tr><td><code>betaln</code></td><td>a, b ∈ [0, 20]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.2e-12</samp></td><td align="right" nowrap><samp>1.3e-14</samp></td><td align="right" nowrap><samp>1.2e-12</samp></td><td align="right" nowrap><samp>1.3e-14</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>betaln</code></td><td>a, b ∈ [10⁻³⁰⁰, 170] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.5e-15</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>8.5e-15</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>betaln</code></td><td>a ∈ ±⁠[10⁻³, 10⁶] (log), b ∈ [10⁻³, 10⁶] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.5e-09</samp></td><td align="right" nowrap><samp>5.1e-11</samp></td><td align="right" nowrap><samp>inf (2)</samp></td><td align="right" nowrap><samp>3.0e-10</samp></td><td align="right" nowrap>41.6%</td></tr>
<tr><td><code>binom</code></td><td>n ∈ [0, 100], k ∈ [0, n]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.7e-14</samp></td><td align="right" nowrap><samp>7.0e-15</samp></td><td align="right" nowrap><samp>3.7e-14</samp></td><td align="right" nowrap><samp>7.0e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>binom</code></td><td>n ∈ {0, …, 1000}, k ∈ {0, …, n}</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.9e-12</samp></td><td align="right" nowrap><samp>3.3e-13</samp></td><td align="right" nowrap><samp>2.3e-12</samp></td><td align="right" nowrap><samp>3.7e-13</samp></td><td align="right" nowrap>88.8%</td></tr>
<tr><td><code>pow_rising</code></td><td>a ∈ [0, 10], m ∈ [−10, 10]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>4.0e-15</samp></td><td align="right" nowrap><samp>6.6e-16</samp></td><td align="right" nowrap><samp>4.0e-15</samp></td><td align="right" nowrap><samp>6.6e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>pow_rising</code></td><td>a ∈ [10⁻³, 100] (log), m ∈ ±⁠[10⁻³, 100] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.5e-15</samp></td><td align="right" nowrap><samp>6.4e-16</samp></td><td align="right" nowrap><samp>1.7e-13</samp></td><td align="right" nowrap><samp>1.2e-14</samp></td><td align="right" nowrap>79.7%</td></tr>
<tr><td><code>pow_rising</code></td><td>a ∈ ±⁠[10⁻³, 10⁶] (log), m ∈ ±⁠[10⁻³, 100] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.2e-13</samp></td><td align="right" nowrap><samp>9.2e-15</samp></td><td align="right" nowrap><samp>1.2e-08</samp></td><td align="right" nowrap><samp>2.6e-10</samp></td><td align="right" nowrap>55.2%</td></tr>
</tbody>
</table>

`gammaln`, `digamma`, `betaln`: The relative error is ill-conditioned near the zeros of the function, which is where the largest errors come from.

`beta`, `betaln`, `pow_rising`: For negative non-integer a, these are ill-conditioned for large |a|, and near the poles of Γ(a + b) (or Γ(a + m)): a relative change ε of a changes Γ(a + b) / Γ(a) by a factor of about 1 + ε a (ψ(a + b) − ψ(a)). This is where the largest errors come from.

## Error function and normal distribution

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>erf</code></td><td>x ∈ [−6, 6]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.4e-16</samp></td><td align="right" nowrap><samp>5.2e-17</samp></td><td align="right" nowrap><samp>3.4e-16</samp></td><td align="right" nowrap><samp>5.2e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erf</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 10] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.9e-16</samp></td><td align="right" nowrap><samp>8.0e-17</samp></td><td align="right" nowrap><samp>2.9e-16</samp></td><td align="right" nowrap><samp>8.0e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erfc</code></td><td>x ∈ [−6, 27]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.7e-14</samp></td><td align="right" nowrap><samp>1.4e-14</samp></td><td align="right" nowrap><samp>5.7e-14</samp></td><td align="right" nowrap><samp>1.4e-14</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erfc</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 27] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.8e-14</samp></td><td align="right" nowrap><samp>3.0e-16</samp></td><td align="right" nowrap><samp>2.8e-14</samp></td><td align="right" nowrap><samp>3.0e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erfcx</code></td><td>x ∈ [−26, 26]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.7e-14</samp></td><td align="right" nowrap><samp>1.1e-14</samp></td><td align="right" nowrap><samp>5.7e-14</samp></td><td align="right" nowrap><samp>1.1e-14</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erfcx</code></td><td>x ∈ [10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.6e-16</samp></td><td align="right" nowrap><samp>6.3e-17</samp></td><td align="right" nowrap><samp>8.6e-16</samp></td><td align="right" nowrap><samp>6.3e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erfinv</code></td><td>y ∈ [−1, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>7.1e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap><samp>7.1e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erfinv</code></td><td>y ∈ 1 − [10⁻¹⁶, 1] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.1e-16</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap><samp>inf (151)</samp></td><td align="right" nowrap><samp>5.4e-04</samp></td><td align="right" nowrap>51.2%</td></tr>
<tr><td><code>erfcinv</code></td><td>y ∈ [0, 2]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>7.4e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap><samp>7.4e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>erfcinv</code></td><td>y ∈ [10⁻³⁰⁰, 1] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>4.2e-16</samp></td><td align="right" nowrap><samp>1.2e-16</samp></td><td align="right" nowrap><samp>4.2e-16</samp></td><td align="right" nowrap><samp>1.2e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ndtr</code></td><td>x ∈ [−10, 10]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.6e-14</samp></td><td align="right" nowrap><samp>2.8e-15</samp></td><td align="right" nowrap><samp>1.6e-14</samp></td><td align="right" nowrap><samp>2.8e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ndtr</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 38] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.2e-13</samp></td><td align="right" nowrap><samp>1.4e-15</samp></td><td align="right" nowrap><samp>1.2e-13</samp></td><td align="right" nowrap><samp>1.4e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ndtri</code></td><td>p ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>7.3e-16</samp></td><td align="right" nowrap><samp>1.2e-16</samp></td><td align="right" nowrap><samp>7.3e-16</samp></td><td align="right" nowrap><samp>1.2e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ndtri</code></td><td>p ∈ [10⁻³⁰⁰, 0.5] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.8e-16</samp></td><td align="right" nowrap><samp>8.8e-17</samp></td><td align="right" nowrap><samp>3.8e-16</samp></td><td align="right" nowrap><samp>8.8e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>log_ndtr</code></td><td>x ∈ [−20, 6]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.7e-15</samp></td><td align="right" nowrap><samp>7.1e-16</samp></td><td align="right" nowrap><samp>6.7e-15</samp></td><td align="right" nowrap><samp>7.1e-16</samp></td><td align="right" nowrap>98.6%</td></tr>
<tr><td><code>log_ndtr</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 10¹⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.7e-14</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap><samp>6.7e-14</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap>98.5%</td></tr>
</tbody>
</table>

## Elementary functions

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th colspan="2">std</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>expm1</code></td><td>x ∈ [−2, 2]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap><samp>5.1e-17</samp></td><td align="right" nowrap><samp>3.9e-16</samp></td><td align="right" nowrap><samp>8.4e-17</samp></td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap><samp>5.1e-17</samp></td><td align="right" nowrap>72.4%</td></tr>
<tr><td><code>expm1</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 709] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>1.1e-17</samp></td><td align="right" nowrap><samp>2.8e-16</samp></td><td align="right" nowrap><samp>2.2e-17</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>1.1e-17</samp></td><td align="right" nowrap>98.1%</td></tr>
<tr><td><code>log1p</code></td><td>x ∈ [−1, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap><samp>5.1e-17</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap><samp>4.9e-17</samp></td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap><samp>5.1e-17</samp></td><td align="right" nowrap>93.6%</td></tr>
<tr><td><code>log1p</code></td><td>x ∈ [10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>3.5e-17</samp></td><td align="right" nowrap><samp>1.9e-16</samp></td><td align="right" nowrap><samp>3.5e-17</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>3.5e-17</samp></td><td align="right" nowrap>99.8%</td></tr>
<tr><td><code>log1p</code></td><td>x ∈ −⁠[10⁻³⁰⁰, 1] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>1.1e-17</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>1.1e-17</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>1.1e-17</samp></td><td align="right" nowrap>99.9%</td></tr>
<tr><td><code>cbrt</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap><samp>5.0e-17</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap>91.3%</td></tr>
<tr><td><code>exp2</code></td><td>x ∈ [−1, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap><samp>2.0e-16</samp></td><td align="right" nowrap><samp>5.4e-17</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap>88.0%</td></tr>
<tr><td><code>exp2</code></td><td>x ∈ [−1022, 1023]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap><samp>1.8e-16</samp></td><td align="right" nowrap><samp>5.5e-17</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap>87.5%</td></tr>
<tr><td><code>exp10</code></td><td>x ∈ [−1, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap><samp>2.0e-16</samp></td><td align="right" nowrap><samp>5.5e-17</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>87.3%</td></tr>
<tr><td><code>exp10</code></td><td>x ∈ [−307, 308]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>4.7e-17</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>5.6e-17</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>86.2%</td></tr>
<tr><td><code>sinpi</code></td><td>x ∈ [−10, 10]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.4e-16</samp></td><td align="right" nowrap><samp>6.3e-17</samp></td><td align="right" nowrap><samp>2.4e-16</samp></td><td align="right" nowrap><samp>6.3e-17</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>sinpi</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 10¹⁵] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap><samp>6.1e-17</samp></td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap><samp>6.1e-17</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>cospi</code></td><td>x ∈ [−10, 10]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.4e-16</samp></td><td align="right" nowrap><samp>6.2e-17</samp></td><td align="right" nowrap><samp>2.4e-16</samp></td><td align="right" nowrap><samp>6.2e-17</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>cospi</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 10¹⁵] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>1.5e-17</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>1.5e-17</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>cosine_cdf</code></td><td>x ∈ [−π, π]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.3e-16</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap><samp>8.3e-16</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>cosine_invcdf</code></td><td>p ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>9.9e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap><samp>9.9e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>99.7%</td></tr>
<tr><td><code>cosine_invcdf</code></td><td>p ∈ [10⁻³⁰⁰, 0.5] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.2e-16</samp></td><td align="right" nowrap><samp>4.8e-17</samp></td><td align="right" nowrap><samp>6.2e-16</samp></td><td align="right" nowrap><samp>4.8e-17</samp></td><td align="right" nowrap></td><td align="right" nowrap></td><td align="right" nowrap>100%</td></tr>
</tbody>
</table>

## Bessel and Struve functions

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>bessel_j0</code></td><td>x ∈ [0, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>9.9e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap><samp>9.9e-16</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_j0</code></td><td>x ∈ [10⁻³⁰⁰, 10⁸] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>9.1e-16</samp></td><td align="right" nowrap><samp>2.4e-17</samp></td><td align="right" nowrap><samp>9.1e-16</samp></td><td align="right" nowrap><samp>2.4e-17</samp></td><td align="right" nowrap>99.8%</td></tr>
<tr><td><code>bessel_j1</code></td><td>x ∈ [0, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.6e-16</samp></td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap><samp>6.6e-16</samp></td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_j1</code></td><td>x ∈ [10⁻³⁰⁰, 10⁸] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.4e-16</samp></td><td align="right" nowrap><samp>3.2e-17</samp></td><td align="right" nowrap><samp>5.4e-16</samp></td><td align="right" nowrap><samp>3.2e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_y0</code></td><td>x ∈ [0, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.2e-15</samp></td><td align="right" nowrap><samp>3.5e-16</samp></td><td align="right" nowrap><samp>3.2e-15</samp></td><td align="right" nowrap><samp>3.5e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_y0</code></td><td>x ∈ [10⁻³⁰⁰, 10⁸] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.3e-15</samp></td><td align="right" nowrap><samp>9.4e-17</samp></td><td align="right" nowrap><samp>1.3e-15</samp></td><td align="right" nowrap><samp>9.4e-17</samp></td><td align="right" nowrap>99.9%</td></tr>
<tr><td><code>bessel_y1</code></td><td>x ∈ [0, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.8e-15</samp></td><td align="right" nowrap><samp>2.0e-16</samp></td><td align="right" nowrap><samp>2.8e-15</samp></td><td align="right" nowrap><samp>2.0e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_y1</code></td><td>x ∈ [10⁻³⁰⁰, 10⁸] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.4e-15</samp></td><td align="right" nowrap><samp>9.5e-17</samp></td><td align="right" nowrap><samp>1.4e-15</samp></td><td align="right" nowrap><samp>9.5e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_i0</code></td><td>x ∈ [−30, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.7e-16</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap><samp>5.7e-16</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_i0</code></td><td>x ∈ [10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.9e-16</samp></td><td align="right" nowrap><samp>6.6e-17</samp></td><td align="right" nowrap><samp>8.9e-16</samp></td><td align="right" nowrap><samp>6.6e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_i1</code></td><td>x ∈ [−30, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.9e-15</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>1.9e-15</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_i1</code></td><td>x ∈ [10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.5e-16</samp></td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap><samp>8.5e-16</samp></td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_k0</code></td><td>x ∈ [0, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap><samp>1.5e-16</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap><samp>1.5e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_k0</code></td><td>x ∈ [10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.0e-15</samp></td><td align="right" nowrap><samp>9.7e-17</samp></td><td align="right" nowrap><samp>1.0e-15</samp></td><td align="right" nowrap><samp>9.7e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_k1</code></td><td>x ∈ [0, 30]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>7.9e-16</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap><samp>7.9e-16</samp></td><td align="right" nowrap><samp>1.3e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_k1</code></td><td>x ∈ [10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.9e-16</samp></td><td align="right" nowrap><samp>5.0e-17</samp></td><td align="right" nowrap><samp>3.9e-16</samp></td><td align="right" nowrap><samp>5.0e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_j</code></td><td>v ∈ [0, 20], x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.5e-14</samp></td><td align="right" nowrap><samp>1.1e-14</samp></td><td align="right" nowrap><samp>8.5e-14</samp></td><td align="right" nowrap><samp>1.1e-14</samp></td><td align="right" nowrap>99.6%</td></tr>
<tr><td><code>bessel_j</code></td><td>v ∈ [−100, 100], x ∈ [10⁻³⁰⁰, 10⁴] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.0e+00</samp></td><td align="right" nowrap><samp>1.4e-02</samp></td><td align="right" nowrap><samp>inf (20)</samp></td><td align="right" nowrap><samp>1.4e-02</samp></td><td align="right" nowrap>99.7%</td></tr>
<tr><td><code>bessel_y</code></td><td>v ∈ [0, 20], x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.7e-15</samp></td><td align="right" nowrap><samp>7.9e-16</samp></td><td align="right" nowrap><samp>6.7e-15</samp></td><td align="right" nowrap><samp>7.9e-16</samp></td><td align="right" nowrap>99.8%</td></tr>
<tr><td><code>bessel_y</code></td><td>v ∈ [−100, 100], x ∈ [10⁻³⁰⁰, 10⁴] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>inf (5)</samp></td><td align="right" nowrap><samp>8.7e-15</samp></td><td align="right" nowrap><samp>inf (5)</samp></td><td align="right" nowrap><samp>8.7e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_i</code></td><td>v ∈ [0, 20], x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.7e-15</samp></td><td align="right" nowrap><samp>3.7e-16</samp></td><td align="right" nowrap><samp>1.7e-15</samp></td><td align="right" nowrap><samp>3.7e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_i</code></td><td>v ∈ [−100, 100], x ∈ [10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.0e+00</samp></td><td align="right" nowrap><samp>1.0e-02</samp></td><td align="right" nowrap><samp>1.0e+00</samp></td><td align="right" nowrap><samp>1.0e-02</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>bessel_k</code></td><td>v ∈ [0, 20], x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.8e-14</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap><samp>5.8e-14</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap>99.8%</td></tr>
<tr><td><code>bessel_k</code></td><td>v ∈ [0, 100], x ∈ [10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>inf (4)</samp></td><td align="right" nowrap><samp>2.9e-15</samp></td><td align="right" nowrap><samp>inf (4)</samp></td><td align="right" nowrap><samp>2.9e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>struve_h</code></td><td>v ∈ [0, 10], x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.3e-13</samp></td><td align="right" nowrap><samp>3.1e-14</samp></td><td align="right" nowrap><samp>5.3e-13</samp></td><td align="right" nowrap><samp>3.1e-14</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>struve_h</code></td><td>v ∈ [−50, 50], x ∈ [10⁻³⁰⁰, 1000] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.6e-13</samp></td><td align="right" nowrap><samp>1.1e-14</samp></td><td align="right" nowrap><samp>1.6e-13</samp></td><td align="right" nowrap><samp>1.1e-14</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>struve_l</code></td><td>v ∈ −k − ½ ± [10⁻¹², 10⁻³] (log), k ∈ {0, …, 20}, x ∈ [0, 100]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.8e-12</samp></td><td align="right" nowrap><samp>7.7e-14</samp></td><td align="right" nowrap><samp>2.5e-02</samp></td><td align="right" nowrap><samp>3.1e-04</samp></td><td align="right" nowrap>94.6%</td></tr>
</tbody>
</table>

`bessel_j0`, `bessel_j1`, `bessel_y0`, `bessel_y1`, `bessel_j`, `bessel_y`, `struve_h`: In the oscillatory region x > |v|, the error is the absolute error divided by the local amplitude √(Jᵥ(x)² + Yᵥ(x)²), if that is smaller than the relative error.

## Zeta functions

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>zeta</code></td><td>x ∈ [1, 20], q ∈ [0, 10]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>9.0e-16</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap><samp>9.0e-16</samp></td><td align="right" nowrap><samp>2.1e-16</samp></td><td align="right" nowrap>99.9%</td></tr>
<tr><td><code>zeta</code></td><td>x ∈ 1 + [10⁻¹⁵, 1000] (log), q ∈ [10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>4.7e-16</samp></td><td align="right" nowrap><samp>6.3e-17</samp></td><td align="right" nowrap><samp>7.3e-02</samp></td><td align="right" nowrap><samp>7.4e-04</samp></td><td align="right" nowrap>98.7%</td></tr>
<tr><td><code>riemann_zeta</code></td><td>x ∈ [−20, 20]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.8e-13</samp></td><td align="right" nowrap><samp>9.3e-15</samp></td><td align="right" nowrap><samp>5.8e-13</samp></td><td align="right" nowrap><samp>9.3e-15</samp></td><td align="right" nowrap>71.8%</td></tr>
<tr><td><code>riemann_zeta</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 1000] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.4e-14</samp></td><td align="right" nowrap><samp>5.8e-16</samp></td><td align="right" nowrap><samp>3.3e-14</samp></td><td align="right" nowrap><samp>5.8e-16</samp></td><td align="right" nowrap>99.6%</td></tr>
<tr><td><code>zetac</code></td><td>x ∈ [−20, 20]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.3e-13</samp></td><td align="right" nowrap><samp>1.8e-15</samp></td><td align="right" nowrap><samp>1.2e-13</samp></td><td align="right" nowrap><samp>1.8e-15</samp></td><td align="right" nowrap>92.7%</td></tr>
<tr><td><code>zetac</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 1000] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.2e-14</samp></td><td align="right" nowrap><samp>5.4e-16</samp></td><td align="right" nowrap><samp>1.0e+00</samp></td><td align="right" nowrap><samp>3.6e-02</samp></td><td align="right" nowrap>99.7%</td></tr>
</tbody>
</table>

`riemann_zeta`, `zetac`: The relative error is ill-conditioned near the zeros of the function, which is where the largest errors come from.

## Exponential integrals

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>expi</code></td><td>x ∈ [−20, 20]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.5e-14</samp></td><td align="right" nowrap><samp>3.8e-16</samp></td><td align="right" nowrap><samp>2.5e-14</samp></td><td align="right" nowrap><samp>3.8e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>expi</code></td><td>x ∈ ±⁠[10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.7e-14</samp></td><td align="right" nowrap><samp>1.9e-16</samp></td><td align="right" nowrap><samp>1.7e-14</samp></td><td align="right" nowrap><samp>1.9e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>expn</code></td><td>n ∈ {0, …, 10}, x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.8e-15</samp></td><td align="right" nowrap><samp>3.3e-16</samp></td><td align="right" nowrap><samp>1.8e-15</samp></td><td align="right" nowrap><samp>3.3e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>expn</code></td><td>n ∈ {0, …, 1000}, x ∈ [10⁻³⁰⁰, 700] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-13</samp></td><td align="right" nowrap><samp>1.3e-15</samp></td><td align="right" nowrap><samp>1.1e-13</samp></td><td align="right" nowrap><samp>1.3e-15</samp></td><td align="right" nowrap>100%</td></tr>
</tbody>
</table>

## Incomplete gamma and beta functions

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>gammainc</code></td><td>a ∈ [0, 20], x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.6e-14</samp></td><td align="right" nowrap><samp>1.3e-15</samp></td><td align="right" nowrap><samp>1.6e-14</samp></td><td align="right" nowrap><samp>1.3e-15</samp></td><td align="right" nowrap>93.8%</td></tr>
<tr><td><code>gammainc</code></td><td>a, x ∈ [10⁻³⁰⁰, 10⁵] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.1e-13</samp></td><td align="right" nowrap><samp>2.0e-14</samp></td><td align="right" nowrap><samp>1.1e-13</samp></td><td align="right" nowrap><samp>2.0e-14</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>gammainc</code></td><td>a ∈ [10, 10⁸] (log), x ∈ a ± 6√a</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.7e-13</samp></td><td align="right" nowrap><samp>7.8e-15</samp></td><td align="right" nowrap><samp>3.6e-01</samp></td><td align="right" nowrap><samp>2.3e-02</samp></td><td align="right" nowrap>94.3%</td></tr>
<tr><td><code>gammaincc</code></td><td>a ∈ [0, 20], x ∈ [0, 50]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.7e-14</samp></td><td align="right" nowrap><samp>2.7e-15</samp></td><td align="right" nowrap><samp>1.7e-14</samp></td><td align="right" nowrap><samp>2.7e-15</samp></td><td align="right" nowrap>93.2%</td></tr>
<tr><td><code>gammaincc</code></td><td>a, x ∈ [10⁻³⁰⁰, 10⁵] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.1e-13</samp></td><td align="right" nowrap><samp>3.6e-15</samp></td><td align="right" nowrap><samp>2.1e-13</samp></td><td align="right" nowrap><samp>3.6e-15</samp></td><td align="right" nowrap>97.7%</td></tr>
<tr><td><code>gammaincc</code></td><td>a ∈ [10, 10⁸] (log), x ∈ a ± 6√a</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.0e-13</samp></td><td align="right" nowrap><samp>6.2e-15</samp></td><td align="right" nowrap><samp>9.7e-07</samp></td><td align="right" nowrap><samp>2.8e-08</samp></td><td align="right" nowrap>94.8%</td></tr>
<tr><td><code>betainc</code></td><td>a, b ∈ [0, 10], x ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.9e-15</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap><samp>8.9e-15</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>betainc</code></td><td>a, b ∈ [10⁻⁵, 1000] (log), x ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.1e-12</samp></td><td align="right" nowrap><samp>5.5e-14</samp></td><td align="right" nowrap><samp>1.6e-10</samp></td><td align="right" nowrap><samp>1.6e-12</samp></td><td align="right" nowrap>94.0%</td></tr>
<tr><td><code>betainc</code></td><td>a, b ∈ [100, 10⁸] (log), x ∈ μ ± 6σ</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>5.7e-07</samp></td><td align="right" nowrap><samp>2.5e-08</samp></td><td align="right" nowrap><samp>7.9e-02</samp></td><td align="right" nowrap><samp>1.3e-03</samp></td><td align="right" nowrap>28.5%</td></tr>
<tr><td><code>fdtr</code></td><td>a, b ∈ [10⁻³, 10⁶] (log), x ∈ [10⁻³⁰, 10³⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.1e-10</samp></td><td align="right" nowrap><samp>4.5e-12</samp></td><td align="right" nowrap><samp>5.4e+01</samp></td><td align="right" nowrap><samp>3.2e+00</samp></td><td align="right" nowrap>55.2%</td></tr>
<tr><td><code>fdtrc</code></td><td>a, b ∈ [10⁻³, 10⁶] (log), x ∈ [10⁻³⁰, 10³⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.2e-09</samp></td><td align="right" nowrap><samp>1.4e-11</samp></td><td align="right" nowrap><samp>6.4e+01</samp></td><td align="right" nowrap><samp>3.2e+00</samp></td><td align="right" nowrap>55.6%</td></tr>
</tbody>
</table>

`betainc`: μ = a / (a + b) and σ² = μ(1 − μ) / (a + b + 1) are the mean and variance of the beta distribution.

## Elliptic integrals

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>ellipk</code></td><td>m ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.5e-16</samp></td><td align="right" nowrap><samp>6.9e-17</samp></td><td align="right" nowrap><samp>2.5e-16</samp></td><td align="right" nowrap><samp>6.9e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipk</code></td><td>m ∈ 1 − [10⁻¹⁶, 1] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.2e-16</samp></td><td align="right" nowrap><samp>8.8e-17</samp></td><td align="right" nowrap><samp>3.2e-16</samp></td><td align="right" nowrap><samp>8.8e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipk</code></td><td>m ∈ −⁠[10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.1e-16</samp></td><td align="right" nowrap><samp>7.0e-17</samp></td><td align="right" nowrap><samp>3.1e-16</samp></td><td align="right" nowrap><samp>7.0e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipe</code></td><td>m ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap><samp>7.3e-17</samp></td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap><samp>7.3e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipe</code></td><td>m ∈ 1 − [10⁻¹⁶, 1] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap><samp>8.8e-17</samp></td><td align="right" nowrap><samp>2.2e-16</samp></td><td align="right" nowrap><samp>8.8e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipe</code></td><td>m ∈ −⁠[10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.5e-16</samp></td><td align="right" nowrap><samp>5.9e-17</samp></td><td align="right" nowrap><samp>6.5e-16</samp></td><td align="right" nowrap><samp>5.9e-17</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipkinc</code></td><td>φ ∈ [0, π/2], m ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.3e-16</samp></td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap><samp>6.3e-16</samp></td><td align="right" nowrap><samp>1.4e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipkinc</code></td><td>φ ∈ ±⁠[10⁻³, 1000] (log), m ∈ 1 − [10⁻³, 10⁵] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.9e-15</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap><samp>2.9e-15</samp></td><td align="right" nowrap><samp>1.7e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipkinc</code></td><td>φ ∈ ±⁠[10⁻³, 10³⁰⁰] (log), m ∈ 1 − [10⁻³, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.3e-15</samp></td><td align="right" nowrap><samp>1.2e-16</samp></td><td align="right" nowrap><samp>inf (619)</samp></td><td align="right" nowrap><samp>1.2e-16</samp></td><td align="right" nowrap>88.0%</td></tr>
<tr><td><code>ellipeinc</code></td><td>φ ∈ [0, π/2], m ∈ [0, 1]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>8.6e-16</samp></td><td align="right" nowrap><samp>1.6e-16</samp></td><td align="right" nowrap><samp>8.6e-16</samp></td><td align="right" nowrap><samp>1.6e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipeinc</code></td><td>φ ∈ ±⁠[10⁻³, 1000] (log), m ∈ 1 − [10⁻³, 10⁵] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>7.8e-16</samp></td><td align="right" nowrap><samp>1.5e-16</samp></td><td align="right" nowrap><samp>7.8e-16</samp></td><td align="right" nowrap><samp>1.5e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>ellipeinc</code></td><td>φ ∈ ±⁠[10⁻³, 10³⁰⁰] (log), m ∈ 1 − [10⁻³, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>6.9e-16</samp></td><td align="right" nowrap><samp>9.5e-17</samp></td><td align="right" nowrap><samp>inf (616)</samp></td><td align="right" nowrap><samp>9.3e-17</samp></td><td align="right" nowrap>89.6%</td></tr>
</tbody>
</table>

## Other functions

<table>
<thead>
<tr><th rowspan="2">Function</th><th rowspan="2">Domain</th><th rowspan="2">N</th><th colspan="2">xsf-rust</th><th colspan="2">xsf</th><th rowspan="2">bit-identical</th></tr>
<tr><th>peak</th><th>rms</th><th>peak</th><th>rms</th></tr>
</thead>
<tbody>
<tr><td><code>hyp0f1</code></td><td>b, z ∈ ±⁠[10⁻³, 10⁶] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>inf (35)</samp></td><td align="right" nowrap><samp>4.5e-11</samp></td><td align="right" nowrap><samp>inf (1937)</samp></td><td align="right" nowrap><samp>4.6e-02</samp></td><td align="right" nowrap>66.9%</td></tr>
<tr><td><code>lambertw</code></td><td>x ∈ [−1/e, 10]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.4e-15</samp></td><td align="right" nowrap><samp>2.0e-16</samp></td><td align="right" nowrap><samp>1.6e-15</samp></td><td align="right" nowrap><samp>2.6e-16</samp></td><td align="right" nowrap>25.8%</td></tr>
<tr><td><code>lambertw</code></td><td>x ∈ [10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>2.2e-15</samp></td><td align="right" nowrap><samp>1.2e-16</samp></td><td align="right" nowrap><samp>2.2e-15</samp></td><td align="right" nowrap><samp>1.5e-16</samp></td><td align="right" nowrap>68.3%</td></tr>
<tr><td><code>spence</code></td><td>x ∈ [0, 2]</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>3.7e-15</samp></td><td align="right" nowrap><samp>6.4e-16</samp></td><td align="right" nowrap><samp>3.7e-15</samp></td><td align="right" nowrap><samp>6.4e-16</samp></td><td align="right" nowrap>100%</td></tr>
<tr><td><code>spence</code></td><td>x ∈ [10⁻³⁰⁰, 10³⁰⁰] (log)</td><td align="right" nowrap>10000</td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap><samp>1.2e-15</samp></td><td align="right" nowrap><samp>1.1e-16</samp></td><td align="right" nowrap>100%</td></tr>
</tbody>
</table>

`hyp0f1`: In the oscillatory region z < 0 with x = 2√−z > |b − 1|, where the function is Γ(b) (x/2)^(1 − b) J_(b−1)(x), the error is the absolute error divided by the local amplitude |Γ(b)| (x/2)^(1 − b) √(J_(b−1)(x)² + Y_(b−1)(x)²), if that is smaller than the relative error.
