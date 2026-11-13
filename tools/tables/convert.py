# /// script
# requires-python = ">=3.15"
# dependencies = ["polars==2.1.0"]
# ///
"""Convert the xsref tables that the tests use into the flat format of `tests/data/xsref`.

Usage (from the repository root, with a checkout of https://github.com/scipy/xsref):

    uv run tools/tables/convert.py XSREF_DIR tests/data/xsref [--no-cap]

Steps:

1. Find the (table, signature) pairs that a test uses: the entries of the `xsref::tests!` macros
   and the direct `xsref::test(...)` calls in `tests/` and `src/` (not the commented-out ones).
2. Read `In_<sig>`, `Out_<sig>`, and `Err_<sig>_other` of each pair, and cast every column to
   f64 on its own (nulls become NaN), like the parquet loader of the tests did. The `fallback`
   column of `Out_<sig>` is the row flag, not an output.
3. Only keep a stratified sample of the `cyl_bessel_i` rows (see `cap_rows`), unless `--no-cap`.
4. Write one `<table>/<sig>.bin` per pair (the format is described in `tests/data/xsref/README.md`),
   `MANIFEST.tsv`, and xsref's `LICENSE`. The tables that `generate.py` has replaced (those with
   header flags) are kept.

The output is deterministic: running it twice gives identical files.
"""

import argparse
import array
import hashlib
import random
import re
import shutil
import struct
import sys
from pathlib import Path
from typing import Final

import polars as pl

ROOT: Final = Path(__file__).resolve().parents[2]
TABLES: Final = Path("src/xsref/tables/scipy_special_tests")

MAGIC: Final = b"XSFT"
VERSION: Final = 1
HEADER: Final = struct.Struct("<4sIIHHHH")
FALLBACK: Final = 1  # row flag: the reference is xsf's own value
COLUMNS: Final = "table\tsig\trows\tsource_rows\tn_in\tn_out\tn_tol\tsha256"

# per signature of `cyl_bessel_i`: rows with all outputs finite, with a NaN output, and with an
# infinite output (a short class gives its remainder to the finite one)
CAPPED: Final = "cyl_bessel_i"
CAP_NAN: Final = 2_000
CAP_INF: Final = 2_000
CAP: Final = 20_000


def _uncommented(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "\n".join(line for line in lines if not line.lstrip().startswith("//"))


def consumers() -> list[tuple[str, str]]:
    """The (table, signature) of each `xsref::tests!` entry and `xsref::test` call."""
    entry = re.compile(r'\b\w+:\s*"(\w+)",\s*"([\w-]+)",')
    call = re.compile(r'\bxsref::test\(\s*"(\w+)",\s*"([\w-]+)",')
    pairs: list[tuple[str, str]] = []
    for path in sorted([*ROOT.glob("tests/**/*.rs"), *ROOT.glob("src/**/*.rs")]):
        source = _uncommented(path)
        pairs += call.findall(source)
        for block in source.split("xsref::tests! {")[1:]:
            # the entries end at the first unindented `}`
            pairs += entry.findall(block.split("\n}", 1)[0])
    return pairs


def bessel_i_overflows() -> set[float]:
    """The orders `v` of `BESSEL_I_OVERFLOWS` in `testing/src/xsref.rs`."""
    source = (ROOT / "testing/src/xsref.rs").read_text()
    match = re.search(r"const BESSEL_I_OVERFLOWS: \[f64; \d+\] = \[([^\]]*)\];", source)
    assert match, "BESSEL_I_OVERFLOWS not found"
    return {float(v.replace("_", "")) for v in match[1].split(",") if v.strip()}


def read_columns(path: Path) -> list[pl.Series]:
    """Every column cast to f64 on its own, with nulls as NaN."""
    df = pl.read_parquet(path)
    return [col.cast(pl.Float64).fill_null(float("nan")) for col in df.get_columns()]


def cap_rows(inputs: list[pl.Series], outputs: list[pl.Series]) -> pl.Series:
    """The row mask of the stratified sample of a `cyl_bessel_i` table.

    The rows of each reference class are drawn uniformly without replacement with
    `random.Random(0)` (finite, then NaN, then infinite), plus every row with a half-integer
    order `v` or one in `BESSEL_I_OVERFLOWS`.
    """
    frame = pl.DataFrame({f"c{i}": col for i, col in enumerate(outputs)})
    nan = frame.select(pl.any_horizontal(pl.all().is_nan())).to_series().to_list()
    inf = frame.select(pl.any_horizontal(pl.all().is_infinite())).to_series().to_list()
    finite: list[int] = []
    with_nan: list[int] = []
    with_inf: list[int] = []
    for i, (is_nan, is_inf) in enumerate(zip(nan, inf, strict=True)):
        (with_nan if is_nan else with_inf if is_inf else finite).append(i)

    n_nan = min(CAP_NAN, len(with_nan))
    n_inf = min(CAP_INF, len(with_inf))
    n_finite = min(CAP - n_nan - n_inf, len(finite))
    rng = random.Random(0)
    keep = {
        *rng.sample(finite, n_finite),
        *rng.sample(with_nan, n_nan),
        *rng.sample(with_inf, n_inf),
    }
    overflows = bessel_i_overflows()
    orders = inputs[0].to_list()
    keep.update(i for i, v in enumerate(orders) if v % 1 == 0.5 or v in overflows)
    return pl.Series([i in keep for i in range(len(nan))])


def convert(src: Path, table: str, sig: str, *, cap: bool) -> tuple[bytes, int]:
    """The flat file of a table and signature, and its number of source rows."""
    inputs = read_columns(src / f"In_{sig}.parquet")
    out = pl.read_parquet(src / f"Out_{sig}.parquet")
    tols = read_columns(src / f"Err_{sig}_other.parquet")

    # the loader of the tests read the outputs as the leading columns
    assert out.columns[-1] == "fallback", out.columns
    assert out.schema["fallback"] == pl.Boolean, out.schema
    fallback = out.get_column("fallback")
    assert fallback.null_count() == 0
    outputs = [
        col.cast(pl.Float64).fill_null(float("nan")) for col in out.get_columns()[:-1]
    ]

    columns = [*inputs, *outputs, *tols]
    n_rows = len(fallback)
    assert all(len(col) == n_rows for col in columns), (table, sig)

    flags = fallback.cast(pl.UInt8) * FALLBACK
    if cap and table == CAPPED:
        mask = cap_rows(inputs, outputs)
        columns = [col.filter(mask) for col in columns]
        flags = flags.filter(mask)

    # row-major, without a text round trip, so signed zeros and NaN payloads are kept
    frame = pl.DataFrame({f"c{i}": col for i, col in enumerate(columns)})
    values = array.array(
        "d", frame.select(pl.concat_list(pl.all())).to_series().explode().to_list()
    )
    if sys.byteorder == "big":
        values.byteswap()

    header = HEADER.pack(
        MAGIC, VERSION, len(flags), len(inputs), len(outputs), len(tols), 0
    )
    return header + values.tobytes() + bytes(flags.to_list()), n_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("xsref", type=Path, help="a checkout of scipy/xsref")
    parser.add_argument("out", type=Path, help="the output directory")
    parser.add_argument("--no-cap", action="store_true", help=f"keep all {CAPPED} rows")
    args = parser.parse_args()
    xsref: Path = args.xsref
    out: Path = args.out

    # a generated table has header flags
    path = out / "MANIFEST.tsv"
    old = path.read_text().splitlines()[1:] if path.exists() else []
    generated = {
        (table, sig): line
        for line in old
        for table, sig in [line.split("\t")[:2]]
        if (file := out / table / f"{sig}.bin").exists()
        and HEADER.unpack_from(file.read_bytes())[6]
    }

    manifest = [COLUMNS]
    total = total_source = 0
    used = consumers()
    for table, sig in sorted(set(used)):
        if line := generated.get((table, sig)):
            manifest.append(line)
            continue
        data, n_source = convert(
            xsref / TABLES / table, table, sig, cap=not args.no_cap
        )
        (out / table).mkdir(parents=True, exist_ok=True)
        (out / table / f"{sig}.bin").write_bytes(data)

        _, _, n_rows, n_in, n_out, n_tol, _ = HEADER.unpack_from(data)
        digest = hashlib.sha256(data).hexdigest()
        manifest.append(
            f"{table}\t{sig}\t{n_rows}\t{n_source}\t{n_in}\t{n_out}\t{n_tol}\t{digest}"
        )
        total += n_rows
        total_source += n_source

    path.write_text("\n".join(manifest) + "\n")
    shutil.copyfile(xsref / "LICENSE", out / "LICENSE")
    n_converted = len(manifest) - 1 - len(generated)
    print(
        f"{len(used)} tests of {len(manifest) - 1} tables; converted {n_converted}"
        f" ({total} of {total_source} rows), kept {len(generated)} generated ones"
    )


if __name__ == "__main__":
    main()
