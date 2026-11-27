"""Tests of `check_loosened_tests.py`, from the repository root:

python3 -m unittest discover -s tools -p 'test_check*.py'
"""

import math
import os
import struct
import subprocess
import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

import check_loosened_tests as c

# the PR that loosened tests of `hyppfq` for `2F0`, and PRs that did not
PR_1817 = "ae623e7e2"
ORDINARY = ["5850f5158", "e5b1564bc", "ffc0c30b9", "60641a787", "93eb196b7"]

TEST = """\
#[test]
fn test_f() {{
    xsref::np_assert_allclose!(&[f(1.0)], &[1.0], rtol = {rtol});
    for (x, y) in [
        (1.0, 2.0),
        (3.0, 4.0),
    ] {{
        assert!((f(x) - y).abs() <= {tol} * y.abs());
    }}
}}

#[test]
fn test_g() {{}}
"""
F = "tests/x.rs"
BEFORE = TEST.format(rtol="1e-14", tol="1e-13")
BIN = "tests/data/xsref/f/d-d.bin"
ROWS: list[tuple[float, ...]] = [(1.0, 2.0, 1e-15, 0.0), (2.0, 3.0, 1e-15, 0.0)]
# a line over 120 characters, with U+2028 (which is not a line break) above `#[cfg(test)]`
LONG_LINE = 'const S: &str = "\u2028\u2028";' + " " * 100 + "\n"

# file contents by path, where `None` deletes the file
type Files = dict[str, str | bytes | None]


def have(rev: str) -> bool:
    return subprocess.run(["git", "cat-file", "-e", rev], check=False).returncode == 0


def table(rows: Sequence[tuple[float, ...]], *, scaled: bool = True) -> bytes:
    """An xsref table of rows `(input, output, tolerance[, scale])`"""
    header = b"XSFT" + struct.pack("<IIHHHH", 1, len(rows), 1, 1, 1, 2 * scaled)
    values = b"".join(struct.pack(f"<{len(row)}d", *row) for row in rows)
    return header + values + bytes(len(rows))


def statuses(rows: Sequence[tuple[float, ...]], status: int, *, ruled: bool) -> bytes:
    """A table of rows `(input, output, tolerance, scale[, rule value])` of one status"""
    flags = 3 | 4 * ruled
    header = b"XSFT" + struct.pack("<IIHHHH", 1, len(rows), 1, 1, 1, flags)
    values = b"".join(struct.pack(f"<{len(row)}d", *row) for row in rows)
    return header + values + bytes(len(rows)) + bytes([status]) * len(rows)


class Repo:
    """A temporary git repository with a base commit and a changed one"""

    def __init__(self, files: Files) -> None:
        tmp = tempfile.TemporaryDirectory()
        self.root = Path(tmp.name)
        self.cleanup = tmp.cleanup
        self.git("init", "-q")
        self.commit(files)

    def git(self, *args: str) -> None:
        env = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t"}
        env |= {"GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}
        subprocess.run(
            ["git", *args], cwd=self.root, check=True, capture_output=True, env=env
        )

    def commit(self, files: Files) -> None:
        for name, content in files.items():
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            if content is None:
                path.unlink()
            elif isinstance(content, str):
                path.write_text(content)
            else:
                path.write_bytes(content)
        self.git("add", "-A")
        self.git("commit", "-q", "-m", "c")


class Synthetic(unittest.TestCase):
    def check(self, before: Files, after: Files) -> list[str]:
        repo = Repo(before)
        self.addCleanup(repo.cleanup)
        repo.commit(after)
        cwd = Path.cwd()
        os.chdir(repo.root)
        try:
            return [f"{f.path}: {f.message}" for f in c.find("HEAD~", "HEAD")]
        finally:
            os.chdir(cwd)

    def one(self, before: str, after: str, expected: list[str], path: str = F) -> None:
        self.assertEqual(self.check({path: before}, {path: after}), expected)

    def test_numbers_with_epsilon_multiples(self) -> None:
        self.assertEqual(
            c.numbers("let tol = 4.0 * f64::EPSILON; // 1e-3"), (4.0 * c.EPSILON,)
        )
        self.assertEqual(c.numbers("f(1_000.5, 2e-3)"), (1000.5, 2e-3))

    def test_loosened(self) -> None:
        for old, new in [
            ("rtol = 1e-14);", "rtol = 1e-13);"),
            ("const RTOL: f64 = 1e-5;", "const RTOL: f64 = 1e-4;"),
            ("close(ab, ba, 1e-10);", "close(ab, ba, 1e-2);"),
            ("e <= 2.0 * f64::EPSILON * y", "e <= 4.0 * f64::EPSILON * y"),
        ]:
            with self.subTest(old):
                (found,) = self.check({F: old + "\n"}, {F: new + "\n"})
                self.assertRegex(found, rf"^{F}: tolerance [\d.e-]+ > [\d.e-]+$")
                self.one(new + "\n", old + "\n", [])

    def test_tolerances_that_are_not_loosened(self) -> None:
        self.one(BEFORE, BEFORE + "\n// c\n", [])
        self.one(BEFORE, TEST.format(rtol="1e-15", tol="1e-14"), [])
        # the output of a row is not a tolerance
        self.one("(1.0, 1e-5),\n", "(1.0, 1e-3),\n", [])

    def test_tolerance_variable(self) -> None:
        self.one(
            "const RTOL: f64 = 1e-12;\n",
            "const RTOL: f64 = DEFAULT_RTOL;\n",
            [f"{F}: tolerance is not a literal anymore"],
        )

    def test_src_tests_only(self) -> None:
        code = (
            LONG_LINE
            + "fn f() {{ let rtol = {}; }}\n#[cfg(test)]\nmod t {{ fn t() {{ let rtol = {}; }} }}\n"
        )
        self.one(
            code.format("1e-14", "1e-14"),
            code.format("1e-3", "1e-13"),
            ["src/a.rs: tolerance 1e-13 > 1e-14"],
            "src/a.rs",
        )

    def test_unicode_path(self) -> None:
        path = "tests/é.rs"
        self.one(
            "let rtol = 1e-14;\n",
            "let rtol = 1e-13;\n",
            [f"{path}: tolerance 1e-13 > 1e-14"],
            path,
        )

    def test_removed_tests(self) -> None:
        no_g = BEFORE.replace("\n#[test]\nfn test_g() {}\n", "")
        self.one(BEFORE, no_g, [f"{F}: removed test `test_g`"])
        # only the attribute
        self.one(
            BEFORE,
            BEFORE.replace("#[test]\nfn test_g", "fn test_g"),
            [f"{F}: removed test `test_g`"],
        )
        # a reordering of attributes and doc comments
        doc = "#[test]\n#[should_panic]\nfn t() {}\n"
        self.one(doc, "#[should_panic]\n/// doc\n#[test]\nfn t() {}\n", [])
        # a new test in the same file excuses a removed one (a rename), also an unrelated one
        two = "#[test]\nfn a() {}\n\n#[test]\nfn b() {}\n"
        self.one(two, two.replace("fn b", "fn c"), [])
        # but not a test that only lost its attribute, even with a new test
        orphan = two.replace("#[test]\nfn b() {}", "fn b() {}\n\n#[test]\nfn c() {}")
        self.one(two, orphan, [f"{F}: removed test `b`"])
        # one new test excuses one removed test
        self.one(two, "#[test]\nfn c() {}\n", [f"{F}: removed test `b`"])
        # a move to another file
        moved: Files = {F: no_g, "tests/y.rs": "#[test]\nfn test_g() {}\n"}
        self.assertEqual(self.check({F: BEFORE}, moved), [])
        # a test of that name in another file that already had it is no move
        smoke = "#[test]\nfn smoke() {}\n"
        a, b = "tests/a.rs", "tests/b.rs"
        self.assertEqual(
            self.check({a: smoke, b: smoke}, {a: "", b: smoke}),
            [f"{a}: removed test `smoke`"],
        )
        # equal names in two modules of one file
        mods = "mod m {\n#[test]\nfn t() {}\n}\nmod n {\n#[test]\nfn t() {}\n}\n"
        one = mods.replace("#[test]\nfn t() {}\n", "", 1)
        self.one(mods, one, [f"{F}: removed test `t`"])
        # a deleted file
        self.assertEqual(
            self.check({F: BEFORE}, {F: None, "tests/y.rs": ""}),
            [
                f"{F}: removed test `test_f`",
                f"{F}: removed test `test_g`",
                f"{F}: removed 2 rows of a test table",
            ],
        )
        # an entry of `xsref::tests!`, keyed by its name: a changed body is no removal
        self.one('a: "x",\nb: "y",\n', 'a: "x",\n', [f"{F}: removed test `b`"])
        self.one('a: "x", f(1)\n', 'a: "x", f(2)\n', [])

    def test_ignored(self) -> None:
        self.one(
            BEFORE,
            BEFORE.replace("#[test]\nfn test_g", "#[test]\n#[ignore]\nfn test_g"),
            [f"{F}: ignored test"],
        )
        # also when moved to a new file
        ignored: Files = {F: "", "tests/new.rs": "#[test]\n#[ignore]\nfn t() {}\n"}
        self.assertEqual(
            self.check({F: "#[test]\nfn t() {}\n"}, ignored),
            ["tests/new.rs: ignored test"],
        )

    def test_removed_rows(self) -> None:
        self.one(
            BEFORE,
            BEFORE.replace("        (3.0, 4.0),\n", ""),
            [f"{F}: removed 1 row of a test table"],
        )
        # replaced by a row with other inputs
        self.one(
            BEFORE,
            BEFORE.replace("(3.0, 4.0)", "(5.0, 4.0)"),
            [f"{F}: removed 1 row of a test table"],
        )
        # a rewritten input is flagged, which is accepted
        row = "(0.6, 1.0, 3.0),\n"
        self.one(
            row, row.replace("0.6", "0.65"), [f"{F}: removed 1 row of a test table"]
        )
        # a new output for the same input
        self.one(BEFORE, BEFORE.replace("(3.0, 4.0)", "(3.0, 5.0)"), [])

    def test_exclusions(self) -> None:
        f = "tests/properties/x.rs"
        base = "property!(no_panic, f(a in positive()));\n"
        found = [f"{f}: property-test exclusion"]
        guards = [
            "property!(no_panic, f(a in positive()) if a < 1e300);\n",
            "property!(checks::no_panic, f(x in positive()) if x > 1.0);\n",
            "property!(no_panic, f(x) if x > 0.0);\n",
        ]
        for guard in guards:
            with self.subTest(guard):
                self.one(base, base + guard, found, f)
        # a match arm is not a guard, and a new file is a new test
        self.one(base, base + "    Some(w) if w > 0 => 1,\n", [], f)
        self.assertEqual(self.check({"README.md": ""}, {f: base + guard}), [])

    def test_xsref(self) -> None:
        loose = [(1.0, 2.0, 1e-14, 0.0), ROWS[1]]
        scaled = [(1.0, 2.0, 1e-15, 1.0), ROWS[1]]
        row = f"{BIN}: larger tolerance or scale, or weaker status in 1 xsref row"
        cases: list[tuple[list[tuple[float, ...]], list[str]]] = [
            (ROWS[::-1], []),
            (loose, [row]),
            (scaled, [row]),
            (ROWS[:1], [f"{BIN}: removed 1 row of an xsref table"]),
        ]
        for rows, expected in cases:
            self.assertEqual(
                self.check({BIN: table(ROWS)}, {BIN: table(rows)}), expected
            )
        # a scale column that is added is a scale
        plain = table([r[:3] for r in ROWS], scaled=False)
        self.assertEqual(
            self.check({BIN: plain}, {BIN: table([(1.0, 2.0, 1e-15, 1e3), ROWS[1]])}),
            [row],
        )
        self.assertEqual(self.check({BIN: plain}, {BIN: table(ROWS)}), [])

    def test_xsref_gap(self) -> None:
        nan = math.nan
        convention = statuses([(1.0, nan, 0.0, 0.0)], 2, ruled=False)
        computed = statuses([(1.0, 2.0, 1e-15, 0.0)], 1, ruled=False)
        found = f"{BIN}: larger tolerance or scale, or weaker status in 1 xsref row"

        def gap(rule: float, tol: float = 1e-8) -> bytes:
            return statuses([(1.0, 2.0, tol, 1.0, rule)], 3, ruled=True)

        # a convention that becomes a gap with the same rule value, whatever its tolerance
        self.assertEqual(self.check({BIN: convention}, {BIN: gap(nan)}), [])
        self.assertEqual(self.check({BIN: gap(nan, 1e-12)}, {BIN: gap(nan, 1e-6)}), [])
        # another rule value, a computed output, or a gap that is a convention for another value
        self.assertEqual(self.check({BIN: gap(nan)}, {BIN: gap(0.0)}), [found])
        self.assertEqual(self.check({BIN: computed}, {BIN: gap(nan)}), [found])
        self.assertEqual(self.check({BIN: convention}, {BIN: gap(0.0)}), [found])
        self.assertEqual(
            self.check(
                {BIN: gap(nan)}, {BIN: statuses([(1.0, 0.0, 0.0, 0.0)], 2, ruled=False)}
            ),
            [found],
        )
        # a computed output that becomes a convention for another reference
        self.assertEqual(
            self.check(
                {BIN: computed}, {BIN: statuses([(1.0, 3.0, 0.0, 0.0)], 2, ruled=False)}
            ),
            [found],
        )
        # a gap without a rule value accepts any wrong output, so a rule value is stricter
        old = statuses([(1.0, 2.0, 1e-8, 1.0)], 3, ruled=False)
        self.assertEqual(self.check({BIN: old}, {BIN: gap(nan)}), [])
        # a gap that closes, or becomes a convention for its rule value
        self.assertEqual(self.check({BIN: gap(nan)}, {BIN: computed}), [])
        loose = statuses([(1.0, 2.0, 1e-7, 0.0)], 1, ruled=False)
        self.assertEqual(self.check({BIN: gap(nan)}, {BIN: loose}), [found])
        self.assertEqual(
            self.check(
                {BIN: gap(nan)}, {BIN: statuses([(1.0, nan, 0.0, 0.0)], 2, ruled=False)}
            ),
            [],
        )

    def test_xsref_gap_subsets(self) -> None:
        nan, zero = math.nan, 0.0
        found = f"{BIN}: larger tolerance or scale, or weaker status in 1 xsref row"
        convention, computed, legacy = 2, 1, 0

        def one(row: tuple[float, ...], status: int, *, ruled: bool = False) -> bytes:
            return statuses([row], status, ruled=ruled)

        def gap(
            ref: float, rule: float, tol: float = 1e-8, scale: float = 0.0
        ) -> bytes:
            return one((1.0, ref, tol, scale, rule), 3, ruled=True)

        cases: list[tuple[bytes, bytes, list[str]]] = [
            # a gap that loses its assertion, or closes at another reference
            (gap(2.0, nan), one((1.0, 2.0, 1e-8, 0.0), legacy), [found]),
            (gap(3.0, nan), one((1.0, 42.0, 1e-8, 0.0), computed), [found]),
            (gap(3.0, nan), one((1.0, 3.0, 1e-8, 0.0), computed), []),
            # a smaller tolerance shrinks what the gap excludes
            (gap(2.0, 1.0, 0.125), gap(2.0, 1.0, 0.0625), [found]),
            (gap(2.0, 1.0, 0.0625), gap(2.0, 1.0, 0.125), []),
            (gap(2.0, 1.0, scale=1.0), gap(2.0, 1.0, scale=0.5), [found]),
            (
                gap(1.0, 1 + 2.0**-24 + 2.0**-52, 2.0**-26),
                gap(1.0, 1 + 2.0**-24 + 2.0**-52, 2.0**-27),
                [found],
            ),
            # a rule-less gap has no rule to lose, but its reference must stay
            (one((1.0, 2.0, 0.125, 0.0), 3), gap(2.0, 1.0, 0.0625), [found]),
            (one((1.0, 2.0, 0.125, 0.0), 3), gap(3.0, 1.0, 0.125), [found]),
            # a wider convention than the rule check of a gap
            (gap(3.0, 2.0), one((1.0, 2.0, 0.01, 0.0), convention), [found]),
            (gap(3.0, 2.0), one((1.0, 2.0, 0.0, 1.0), convention), [found]),
            (gap(3.0, 2.0), one((1.0, 2.0, 0.0, 0.0), convention), []),
            # an output that is checked for more, on the same reference
            (
                one((1.0, zero, 1e-15, 0.0), computed),
                one((1.0, zero, 0.0, 0.0), convention),
                [],
            ),
            (
                one((1.0, nan, 0.0, 0.0), convention),
                gap(1.0, -nan),
                [],
            ),
            (one((1.0, zero, 0.0, 0.0), legacy), gap(1.0, zero), []),
            (one((1.0, zero, 0.0, 0.0), legacy), gap(1.0, -zero), [found]),
        ]
        for old, new, expected in cases:
            with self.subTest(old=old[20:], new=new[20:]):
                self.assertEqual(self.check({BIN: old}, {BIN: new}), expected)

    def test_parse_table(self) -> None:
        row = struct.pack("<5d", 1.0, 2.0, 1e-8, 0.0, math.nan)
        gap = statuses([(1.0, 2.0, 1e-8, 0.0, math.nan)], 3, ruled=True)
        self.assertEqual(c.parse_table(gap)[struct.pack("<d", 1.0)][0].status, 3)
        # rule values without statuses
        header = b"XSFT" + struct.pack("<IIHHHH", 1, 1, 1, 1, 1, 2 | 4)
        invalid = [
            header + row + bytes(1),
            # statuses that are missing, or one byte too many
            gap[:-1],
            gap + bytes(1),
            table(ROWS)[:-1],
            table(ROWS) + bytes(1),
        ]
        for data in invalid:
            with self.subTest(data=data[8:]), self.assertRaises(ValueError):
                c.parse_table(data)

    def test_xsref_removed(self) -> None:
        removed = [f"{BIN}: removed xsref table"]
        new = "tests/data/xsref/f/d-d_d.bin"
        cases: list[tuple[Files, list[str]]] = [
            ({BIN: None}, removed),
            # only a table with the same rows replaces a table
            ({BIN: None, "tests/data/xsref/f/README.md": "x"}, removed),
            ({BIN: None, new: table(ROWS[:1])}, removed),
            ({BIN: None, new: table(ROWS)}, []),
        ]
        for after, expected in cases:
            self.assertEqual(self.check({BIN: table(ROWS)}, after), expected)

    def test_unreadable_head_table(self) -> None:
        # a changed table whose blob can't be read is an error, not a removal
        repo = Repo({BIN: table(ROWS)})
        self.addCleanup(repo.cleanup)
        repo.commit({BIN: table(ROWS[:1])})
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(repo.root)
        with self.assertRaises(subprocess.CalledProcessError):
            c.check_tables("HEAD~", "no-such-rev", [("M", BIN)])


class History(unittest.TestCase):
    # these need the full history (they skip without it)
    def test_1817(self) -> None:
        if not have(PR_1817):
            self.skipTest("no history")
        found = c.find(f"{PR_1817}^", PR_1817)
        self.assertIn("tolerance is not a literal anymore", [f.message for f in found])
        self.assertTrue(any(f.path == "tests/ext/pfq.rs" for f in found))

    def test_ordinary(self) -> None:
        for rev in ORDINARY:
            if not have(rev):
                self.skipTest("no history")
            with self.subTest(rev=rev):
                self.assertEqual(c.find(f"{rev}^", rev), [])


if __name__ == "__main__":
    unittest.main()
