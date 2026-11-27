"""Flag the changes of a PR that make the tests weaker.

Usage (from the repository root):

    python3 tools/check_loosened_tests.py [BASE [HEAD]]

Diffs `HEAD` (default: `HEAD`) against its merge base with `BASE` (default: `origin/master`), and
reports with file and line:

- a larger number below 0.1 (a tolerance, e.g. `rtol = 1e-14`, `const RTOL: f64 = 1e-5`,
  `close(a, b, 1e-10)`, `<= 1e-13 * x.abs()`) in a line that replaces a line of the same shape,
  or a tolerance that is not a literal anymore; this can also be another literal below 0.1 (the
  label allows it), but not one in a table row
- a new property-test exclusion (`prop_assume!`, `prop_filter`, `if` or `unless` clause of a
  property, `// known:` comment) in `tests/properties`
- a new `#[ignore]`
- a removed test, which no added test of its file replaces
- removed rows of a test table, whose inputs no added row has
- a removed xsref table, or a removed row or a larger tolerance or scale of one
  (`tests/data/xsref`)

Test code is everything in `tests/` and `testing/`, and the `#[cfg(test)]` part of `src/`. The exit
status is 1 if there are findings. On GitHub Actions, they are also printed as annotations and
appended to the job summary, and `found=true|false` is written to `GITHUB_OUTPUT` (exit status 0).
"""

import os
import re
import struct
import subprocess
import sys
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from functools import cache
from pathlib import Path

EPSILON = 2.0**-52
# the largest tolerance that a comparison `.abs() < NUM` is taken to be
SMALL = 0.1
XSREF = "tests/data/xsref/"
PROPERTIES = "tests/properties/"

NUM = r"(?>\d[\d_]*(?:\.[\d_]*)?(?:[eE][+-]?\d+)?)"
# `NUM * EPSILON` and `EPSILON * NUM`
EPSILON_MULTIPLE = re.compile(
    rf"(?:(?P<before>{NUM})\s*\*\s*(?:f64::)?EPSILON\b|\bEPSILON\s*\*\s*(?P<after>{NUM}))"
)
# a tolerance that is a named literal, or a name that is no literal
TOLERANCE_NAME = r"\b\w*tol\w*\b(?:\s*:\s*\w+)?\s*=\s*"
TOLERANCE_LITERAL = re.compile(TOLERANCE_NAME + NUM, re.IGNORECASE)
TOLERANCE_VARIABLE = re.compile(TOLERANCE_NAME + "[a-z_]", re.IGNORECASE)
# the `if` guard of `property!(check, name(args) if guard)` or `name(args) if guard`
PROPERTY_GUARD = r"^\s*(?:property!\(\w+,\s*)?[a-z_]\w*\(.*\)\s+if\b"
EXCLUSION = re.compile(
    rf"prop_assume!|prop_filter|\bunless\b|{PROPERTY_GUARD}|//\s*known:"
)
ATTRIBUTE = re.compile(r"\s*(?:#\[|//)")
TEST_FN = re.compile(r"^\s*(?:pub\s+)?fn\s+(\w+)")
# an entry of `xsref::tests!`
TEST_ENTRY = re.compile(r"^\s*(\w+):\s*\"")
# the words of a table row that are numbers
ROW_WORD = re.compile(r"f64::[A-Z_]+|_f64|inf|nan|NAN|INFINITY")
ROW = re.compile(r"\s*[(\[][\d\s.,_+\-eE()\[\]]*\d[\d\s.,_+\-eE()\[\]]*[)\]],?\s*")
TOKEN = re.compile(NUM)


@dataclass(frozen=True)
class Finding:
    path: str
    line: int
    message: str


@dataclass(frozen=True)
class Line:
    """A line of a hunk: `kind` is `-` or `+` (context lines are dropped), `no` is in the old or new file"""

    kind: str
    no: int
    text: str
    tested: bool  # `#[test]` is among the attributes right above it
    block: int  # the run of changed lines that it is part of


# unquoted paths
GIT = ["git", "-c", "core.quotePath=false"]


def git(*args: str) -> bytes | None:
    """The output of a git command, or `None` if it fails"""
    out = subprocess.run(GIT + list(args), check=False, capture_output=True)
    return out.stdout if out.returncode == 0 else None


def git_text(*args: str) -> str:
    out = subprocess.run(GIT + list(args), check=True, capture_output=True)
    return out.stdout.decode(errors="replace")


def after_test(tested: bool, text: str) -> bool:
    """Whether `#[test]` is among the attributes right above the line after `text`"""
    if text.strip() == "#[test]":
        return True
    return tested and ATTRIBUTE.match(text) is not None


def hunks(diff: str) -> Iterator[tuple[str, bool, list[Line]]]:
    """The `(path, whether it is new, changed lines)` of each hunk of a unified diff"""
    path = ""
    is_new = False
    lines: list[Line] = []
    old = new = 0
    tested_old = tested_new = False
    header = True
    block = 0
    # not `splitlines`, which also splits at U+2028 and others
    for text in diff.split("\n"):
        if text.startswith("diff --git "):
            header = True
            if lines:
                yield path, is_new, lines
            lines = []
            is_new = False
            path = text.rsplit(" b/", 1)[1]
        elif text.startswith("new file mode"):
            is_new = True
        elif text.startswith("@@"):
            header = False
            if lines:
                yield path, is_new, lines
            lines = []
            m = re.match(r"@@ -(\d+)(?:,\d+)? \+(\d+)", text)
            assert m
            old, new = int(m[1]), int(m[2])
            tested_old = tested_new = False
            block = 0
        elif header:
            continue
        elif text.startswith("-"):
            lines.append(Line("-", old, text[1:], tested_old, block))
            tested_old = after_test(tested_old, text[1:])
            old += 1
        elif text.startswith("+"):
            lines.append(Line("+", new, text[1:], tested_new, block))
            tested_new = after_test(tested_new, text[1:])
            new += 1
        elif text.startswith(" "):
            block += 1
            tested_old = tested_new = after_test(tested_old, text[1:])
            old += 1
            new += 1
    if lines:
        yield path, is_new, lines


@cache
def test_start(rev: str, path: str) -> int:
    """The first line of the test code of a file, which is the whole file in a test directory"""
    if path.startswith(("tests/", "testing/")):
        return 0
    text = (git("show", f"{rev}:{path}") or b"").decode(errors="replace")
    for i, line in enumerate(text.splitlines(), 1):
        if line.strip() == "#[cfg(test)]":
            return i
    return sys.maxsize


def is_row(text: str) -> bool:
    return ROW.fullmatch(ROW_WORD.sub("", text.split("//")[0])) is not None


def numbers(text: str) -> tuple[float, ...]:
    """The number literals of a line, where a multiple of `EPSILON` is its value"""

    def multiple(m: re.Match[str]) -> str:
        return repr(float((m["before"] or m["after"]).replace("_", "")) * EPSILON)

    text = EPSILON_MULTIPLE.sub(multiple, text.split("//")[0])
    return tuple(float(t.replace("_", "")) for t in TOKEN.findall(text))


def test_name(text: str, tested: bool) -> re.Match[str] | None:
    """The name of a test: a `#[test] fn`, or an entry of `xsref::tests!`"""
    if m := TEST_ENTRY.match(text):
        return m
    return TEST_FN.match(text) if tested else None


def shape(text: str) -> str:
    """A line without its numbers, to find the line that replaces it"""
    return re.sub(NUM, "#", "".join(text.split()))


def loosened(path: str, removed: list[Line], added: list[Line]) -> Iterator[Finding]:
    """The larger small numbers (tolerances) of the lines that replace a line of the same shape"""
    old = [(shape(ln.text), numbers(ln.text), ln) for ln in removed]
    for ln in added:
        for i, (old_shape, old_nums, old_ln) in enumerate(old):
            if old_ln.block == ln.block and old_shape == shape(ln.text):
                del old[i]
                if is_row(ln.text):
                    break
                new_nums = numbers(ln.text)
                for x, y in zip(old_nums, new_nums, strict=False):
                    if 0 < x < SMALL and y > x * (1 + 1e-9):
                        yield Finding(path, ln.no, f"tolerance {y:.3g} > {x:.3g}")
                break
        else:
            # a tolerance that became a variable, whose value is unknown
            if TOLERANCE_VARIABLE.search(ln.text) and any(
                old_ln.block == ln.block and TOLERANCE_LITERAL.search(old_ln.text)
                for _, _, old_ln in old
            ):
                yield Finding(path, ln.no, "tolerance is not a literal anymore")


def orphans(removed: list[Line]) -> Iterator[Line]:
    """The removed `#[test]` attributes whose function is not removed"""
    by_no = {ln.no: ln for ln in removed}
    for ln in removed:
        if ln.text.strip() != "#[test]":
            continue
        no = ln.no + 1
        while (nxt := by_no.get(no)) and ATTRIBUTE.match(nxt.text):
            no += 1
        if not nxt or not TEST_FN.match(nxt.text):
            yield ln


def check_hunks(diff: str, base: str, head: str) -> list[Finding]:
    all_hunks = list(hunks(diff))
    added_tests: Counter[str] = Counter()  # by file
    added_names: set[str] = set()
    # a row with the same inputs as a removed one only gets other outputs
    added_rows: set[tuple[float, ...]] = set()
    for path, _, ls in all_hunks:
        for ln in ls:
            if ln.kind == "+":
                added_rows.add(numbers(ln.text)[:-1])
                if m := test_name(ln.text, ln.tested):
                    added_tests[path] += 1
                    added_names.add(m[1])

    gone_tests: dict[str, list[Finding]] = {}
    out: list[Finding] = []
    for path, is_new, ls in all_hunks:
        removed = [
            ln for ln in ls if ln.kind == "-" and ln.no >= test_start(base, path)
        ]
        added = [ln for ln in ls if ln.kind == "+" and ln.no >= test_start(head, path)]

        out.extend(loosened(path, removed, added))

        for ln in added:
            if "#[ignore" in ln.text:
                out.append(Finding(path, ln.no, "ignored test"))
            if not is_new and path.startswith(PROPERTIES) and EXCLUSION.search(ln.text):
                out.append(Finding(path, ln.no, "property-test exclusion"))

        gone_tests.setdefault(path, []).extend(
            Finding(path, ln.no, f"removed test `{m[1]}`")
            for ln in removed
            if (m := test_name(ln.text, ln.tested)) and m[1] not in added_names
        )

        out.extend(
            Finding(path, ln.no, "removed `#[test]` of a function")
            for ln in orphans(removed)
        )

        gone = [
            ln
            for ln in removed
            if is_row(ln.text) and numbers(ln.text)[:-1] not in added_rows
        ]
        n_rows = sum(1 for ln in added if is_row(ln.text))
        n_removed = sum(1 for ln in removed if is_row(ln.text))
        if gone and n_removed > n_rows:
            out.append(
                Finding(
                    path,
                    gone[0].no,
                    f"removed {plural(len(gone), 'row')} of a test table",
                )
            )
    # a test that is renamed in its file is not removed
    for path, found in gone_tests.items():
        out.extend(found[added_tests[path] :])
    return out


def parse_table(data: bytes) -> dict[bytes, tuple[float, ...]] | None:
    """The tolerances and scales of each row of an xsref table, by the bytes of its inputs

    The format is described in `tests/data/xsref/README.md`.
    """
    if data[:4] != b"XSFT" or len(data) < 20:
        return None
    n_rows, n_in, n_out, n_tol, flags = struct.unpack_from("<IHHHH", data, 8)
    width = n_in + n_out + n_tol * (1 + bool(flags & 2))
    if len(data) < 20 + n_rows * width * 8:
        return None
    values = struct.unpack_from(f"<{n_rows * width}d", data, 20)
    rows: dict[bytes, tuple[float, ...]] = {}
    for r in range(n_rows):
        row = values[r * width : (r + 1) * width]
        # the tolerances, then the scales (if any)
        tols = row[n_in + n_out : n_in + n_out + n_tol]
        # like the reader of the tables, without scales there are none (0)
        scales = row[n_in + n_out + n_tol :] or (0.0,) * n_tol
        rows[struct.pack(f"<{n_in}d", *row[:n_in])] = tols + scales
    return rows


def check_tables(base: str, head: str, names: list[str]) -> list[Finding]:
    out: list[Finding] = []
    for name in names:
        old_data = git("show", f"{base}:{name}")
        if old_data is None:
            msg = f"{name} is not in {base}"
            raise RuntimeError(msg)
        new_data = git("show", f"{head}:{name}")
        old = parse_table(old_data)
        new = parse_table(new_data) if new_data else None
        if old is None:
            continue
        if new is None:
            out.append(Finding(name, 1, "removed xsref table"))
            continue
        missing = sum(k not in new for k in old)
        looser = sum(k in new and larger(new[k], tol) for k, tol in old.items())
        if missing:
            out.append(
                Finding(name, 1, f"removed {plural(missing, 'row')} of an xsref table")
            )
        if looser:
            out.append(
                Finding(
                    name,
                    1,
                    f"larger tolerance or scale in {plural(looser, 'xsref row')}",
                )
            )
    return out


def plural(n: int, noun: str) -> str:
    return f"{n} {noun}" + "s" * (n != 1)


def larger(new: tuple[float, ...], old: tuple[float, ...]) -> bool:
    """Whether a tolerance or scale grew"""
    return any(n > o for n, o in zip(new, old, strict=False))


def find(base: str, head: str) -> list[Finding]:
    rng = f"{base}...{head}"
    merge_base = git_text("merge-base", base, head).strip()
    head = git_text("rev-parse", head).strip()
    diff = git_text(
        "diff", "--no-color", "--no-renames", "-U3", rng, "--", ":(top)*.rs"
    )
    changed = [
        line.split("\t")
        for line in git_text(
            "diff", "--name-status", "--no-renames", rng, "--", f":(top){XSREF}"
        ).split("\n")
        if line
    ]
    # a table that is replaced by another `.bin` of its function is not removed
    added = {Path(n).parent for s, n in changed if s == "A" and n.endswith(".bin")}
    tables = [
        n
        for s, n in changed
        if n.endswith(".bin") and s != "A" and (s != "D" or Path(n).parent not in added)
    ]
    findings = check_hunks(diff, merge_base, head)
    return findings + check_tables(merge_base, head, tables)


def escape(text: str) -> str:
    return text.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")


def main(argv: list[str]) -> int:
    base = argv[1] if len(argv) > 1 else "origin/master"
    head = argv[2] if len(argv) > 2 else "HEAD"
    findings = sorted(find(base, head), key=lambda f: (f.path, f.line))
    for f in findings:
        print(f"{f.path}:{f.line}: {f.message}")
        if os.environ.get("GITHUB_ACTIONS"):
            print(f"::warning file={f.path},line={f.line}::{escape(f.message)}")
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(summary).open("a", encoding="utf-8") as fh:
            fh.write("### Loosened tests\n\n")
            if findings:
                fh.write("| Location | Change |\n| --- | --- |\n")
                fh.writelines(
                    f"| `{f.path}:{f.line}` | {f.message} |\n" for f in findings
                )
            else:
                fh.write("None.\n")
    # in CI, the findings are an output (so that only they, not a crash, can be allowed)
    if output := os.environ.get("GITHUB_OUTPUT"):
        with Path(output).open("a", encoding="utf-8") as fh:
            fh.write(f"found={str(bool(findings)).lower()}\n")
        return 0
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
