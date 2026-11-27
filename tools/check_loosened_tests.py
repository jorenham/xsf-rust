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
- fewer tests (`#[test]` or an entry of `xsref::tests!`) than before
- removed rows of a test table, whose inputs no added row has
- a removed xsref table, or a removed row or a larger tolerance or scale of one
  (`tests/data/xsref`)

Test code is everything in `tests/` and `testing/`, and the `#[cfg(test)]` part of `src/`. The exit
status is 1 if there are findings. On GitHub Actions, they are also printed as annotations, and
`found=true|false` is written to `GITHUB_OUTPUT` (exit status 0).
"""

import os
import re
import struct
import subprocess
import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
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
PROPERTY_GUARD = r"^\s*(?:property!\([\w:]+,\s*)?[a-z_]\w*\(.*\)\s+if\b"
EXCLUSION = re.compile(
    rf"prop_assume!|prop_filter|\bunless\b|{PROPERTY_GUARD}|//\s*known:"
)
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
    no: int  # in the old file for a removed line, else in the new file
    text: str


@dataclass
class Hunk:
    """A run of changed lines"""

    path: str
    is_new: bool
    removed: list[Line] = field(default_factory=list)
    added: list[Line] = field(default_factory=list)


def git(*args: str) -> bytes:
    # unquoted paths
    cmd = ["git", "-c", "core.quotePath=false", *args]
    return subprocess.run(cmd, check=True, capture_output=True).stdout


def hunks(diff: str) -> list[Hunk]:
    """The hunks of a `git diff -U0`"""
    out: list[Hunk] = []
    path, is_new, hunk = "", False, None
    old = new = 0
    # not `splitlines`, which also splits at U+2028 and others
    for text in diff.split("\n"):
        if text.startswith("diff --git "):
            path, is_new, hunk = text.rsplit(" b/", 1)[-1], False, None
        elif text.startswith("new file mode"):
            is_new = True
        elif m := re.match(r"@@ -(\d+)(?:,\d+)? \+(\d+)", text):
            old, new = int(m[1]), int(m[2])
            hunk = Hunk(path, is_new)
            out.append(hunk)
        elif hunk and text.startswith("-"):
            hunk.removed.append(Line(old, text[1:]))
            old += 1
        elif hunk and text.startswith("+"):
            hunk.added.append(Line(new, text[1:]))
            new += 1
    return out


@cache
def test_start(rev: str, path: str) -> int:
    """The first line of the test code of a file, which is the whole file in a test directory"""
    if path.startswith(("tests/", "testing/")):
        return 0
    text = git("show", f"{rev}:{path}").decode(errors="replace")
    for i, line in enumerate(text.split("\n"), 1):
        if line.strip() == "#[cfg(test)]":
            return i
    return sys.maxsize


def test_lines(rev: str, path: str, lines: list[Line]) -> list[Line]:
    """The lines that are test code"""
    return [ln for ln in lines if ln.no >= test_start(rev, path)] if lines else []


def is_row(text: str) -> bool:
    return ROW.fullmatch(ROW_WORD.sub("", text.split("//")[0])) is not None


def is_test(text: str) -> bool:
    return text.strip() == "#[test]" or TEST_ENTRY.match(text) is not None


def numbers(text: str) -> tuple[float, ...]:
    """The number literals of a line, where a multiple of `EPSILON` is its value"""

    def multiple(m: re.Match[str]) -> str:
        return repr(float((m["before"] or m["after"]).replace("_", "")) * EPSILON)

    text = EPSILON_MULTIPLE.sub(multiple, text.split("//")[0])
    return tuple(float(t.replace("_", "")) for t in TOKEN.findall(text))


def shape(text: str) -> str:
    """A line without its numbers, to find the line that replaces it"""
    return re.sub(NUM, "#", "".join(text.split()))


def loosened(path: str, removed: list[Line], added: list[Line]) -> Iterator[Finding]:
    """The larger small numbers (tolerances) of the lines that replace a line of the same shape"""
    old = list(removed)
    for ln in added:
        match = next((o for o in old if shape(o.text) == shape(ln.text)), None)
        if match:
            old.remove(match)
            if not is_row(ln.text):
                for x, y in zip(numbers(match.text), numbers(ln.text), strict=False):
                    if 0 < x < SMALL and y > x * (1 + 1e-9):
                        yield Finding(path, ln.no, f"tolerance {y:.3g} > {x:.3g}")
        elif TOLERANCE_VARIABLE.search(ln.text) and any(
            TOLERANCE_LITERAL.search(o.text) for o in old
        ):
            # a tolerance that became a variable, whose value is unknown
            yield Finding(path, ln.no, "tolerance is not a literal anymore")


def check_code(diff: str, base: str, head: str) -> list[Finding]:
    out: list[Finding] = []
    gone_tests: list[tuple[str, Line]] = []
    n_added_tests = 0
    gone_rows: list[tuple[str, Line]] = []
    # a row with the same inputs as a removed one only gets other outputs
    added_rows: set[tuple[float, ...]] = set()
    for h in hunks(diff):
        removed = test_lines(base, h.path, h.removed)
        added = test_lines(head, h.path, h.added)
        out.extend(loosened(h.path, removed, added))
        for ln in added:
            if "#[ignore" in ln.text:
                out.append(Finding(h.path, ln.no, "ignored test"))
            # a new file is a new test, which is not weakened
            if (
                not h.is_new
                and h.path.startswith(PROPERTIES)
                and EXCLUSION.search(ln.text)
            ):
                out.append(Finding(h.path, ln.no, "property-test exclusion"))
        gone_tests += [(h.path, ln) for ln in removed if is_test(ln.text)]
        n_added_tests += sum(is_test(ln.text) for ln in added)
        gone_rows += [(h.path, ln) for ln in removed if is_row(ln.text)]
        added_rows |= {numbers(ln.text)[:-1] for ln in added if is_row(ln.text)}

    # a renamed or moved test is not removed
    if (n := len(gone_tests) - n_added_tests) > 0:
        path, ln = gone_tests[0]
        out.append(Finding(path, ln.no, f"removed {plural(n, 'test')}"))
    if gone := [r for r in gone_rows if numbers(r[1].text)[:-1] not in added_rows]:
        path, ln = gone[0]
        out.append(
            Finding(path, ln.no, f"removed {plural(len(gone), 'row')} of a test table")
        )
    return out


def parse_table(data: bytes) -> dict[bytes, tuple[float, ...]]:
    """The tolerances and scales of each row of an xsref table, by the bytes of its inputs

    The format is described in `tests/data/xsref/README.md`.
    """
    if data[:4] != b"XSFT":
        msg = "not an xsref table"
        raise ValueError(msg)
    n_rows, n_in, n_out, n_tol, flags = struct.unpack_from("<IHHHH", data, 8)
    width = n_in + n_out + n_tol * (1 + bool(flags & 2))
    values = struct.unpack_from(f"<{n_rows * width}d", data, 20)
    rows: dict[bytes, tuple[float, ...]] = {}
    for r in range(n_rows):
        row = values[r * width : (r + 1) * width]
        tols = row[n_in + n_out : n_in + n_out + n_tol]
        # like the reader of the tables, without scales there are none (0)
        scales = row[n_in + n_out + n_tol :] or (0.0,) * n_tol
        rows[struct.pack(f"<{n_in}d", *row[:n_in])] = tols + scales
    return rows


def check_tables(base: str, head: str, changed: list[tuple[str, str]]) -> list[Finding]:
    """The findings of the changed xsref tables, as `(git status, path)`"""
    out: list[Finding] = []
    # a table whose bytes reappear in an added `.bin` was moved, not removed
    added = {
        git("show", f"{head}:{name}")
        for status, name in changed
        if status == "A" and name.endswith(".bin")
    }
    for status, name in changed:
        if status == "A" or not name.endswith(".bin"):
            continue
        old_data = git("show", f"{base}:{name}")
        if status == "D":
            if old_data not in added:
                out.append(Finding(name, 1, "removed xsref table"))
            continue
        old, new = parse_table(old_data), parse_table(git("show", f"{head}:{name}"))
        if missing := sum(k not in new for k in old):
            out.append(
                Finding(name, 1, f"removed {plural(missing, 'row')} of an xsref table")
            )
        # a larger tolerance or scale
        if looser := sum(
            k in new and any(n > o for n, o in zip(new[k], v, strict=True))
            for k, v in old.items()
        ):
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


def find(base: str, head: str) -> list[Finding]:
    merge_base = git("merge-base", base, head).decode().strip()
    diff = git("diff", "--text", "--no-renames", "-U0", merge_base, head, "--", "*.rs")
    changed = [
        (status, path)
        for line in git(
            "diff", "--name-status", "--no-renames", merge_base, head, "--", XSREF
        )
        .decode()
        .splitlines()
        for status, path in [line.split("\t")]
    ]
    return check_code(diff.decode(errors="replace"), merge_base, head) + check_tables(
        merge_base, head, changed
    )


def main(argv: list[str]) -> int:
    base = argv[1] if len(argv) > 1 else "origin/master"
    head = argv[2] if len(argv) > 2 else "HEAD"
    findings = sorted(find(base, head), key=lambda f: (f.path, f.line))
    for f in findings:
        print(f"{f.path}:{f.line}: {f.message}")
        if os.environ.get("GITHUB_ACTIONS"):
            print(f"::warning file={f.path},line={f.line}::{f.message}")
    # in CI, the findings are an output (so that only they, not a crash, can be allowed)
    if output := os.environ.get("GITHUB_OUTPUT"):
        with Path(output).open("a", encoding="utf-8") as fh:
            fh.write(f"found={str(bool(findings)).lower()}\n")
        return 0
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
