"""Tests of `check_loosened_tests.py`, from the repository root:

python3 -m unittest discover -s tools -p 'test_check*.py'
"""

import os
import struct
import subprocess
import tempfile
import unittest
from collections.abc import Sequence
from pathlib import Path

import check_loosened_tests as c

# the PR that loosened tests of `hyppfq` for `2F0`, and test-only or ordinary PRs that did not
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


# file contents by path, where `None` deletes the file
type Files = dict[str, str | bytes | None]


def have(rev: str) -> bool:
    return subprocess.run(["git", "cat-file", "-e", rev], check=False).returncode == 0


def table(rows: Sequence[tuple[float, ...]], *, scaled: bool = True) -> bytes:
    """An xsref table of rows `(input, output, tolerance[, scale])`"""
    header = b"XSFT" + struct.pack("<IIHHHH", 1, len(rows), 1, 1, 1, 2 * scaled)
    values = b"".join(struct.pack(f"<{len(row)}d", *row) for row in rows)
    return header + values + bytes(len(rows))


class Repo:
    """A temporary git repository with a base commit `base` and a changed one"""

    def __init__(self, files: Files) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
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

    def findings(self) -> list[str]:
        cwd = Path.cwd()
        os.chdir(self.root)
        try:
            return [f"{f.path}: {f.message}" for f in c.find("HEAD~", "HEAD")]
        finally:
            os.chdir(cwd)


class Synthetic(unittest.TestCase):
    def check(self, before: Files, after: Files) -> list[str]:
        repo = Repo(before)
        self.addCleanup(repo.tmp.cleanup)
        repo.commit(after)
        return repo.findings()

    def test_numbers(self) -> None:
        self.assertEqual(c.numbers("rtol = 1e-14);"), (1e-14,))
        self.assertEqual(
            c.numbers("let tol = 4.0 * f64::EPSILON; // 1e-3"), (4.0 * c.EPSILON,)
        )
        self.assertEqual(c.numbers("f(1_000.5, 2e-3)"), (1000.5, 2e-3))

    def test_unchanged(self) -> None:
        f = "tests/x.rs"
        same = TEST.format(rtol="1e-14", tol="1e-13")
        self.assertEqual(self.check({f: same}, {f: same + "\n// c\n"}), [])
        self.assertEqual(
            self.check({f: same}, {f: TEST.format(rtol="1e-15", tol="1e-14")}), []
        )

    def test_loosened(self) -> None:
        f = "tests/x.rs"
        found = self.check(
            {f: TEST.format(rtol="1e-14", tol="1e-13")},
            {f: TEST.format(rtol="1e-13", tol="1e-12")},
        )
        self.assertEqual(len(found), 2, found)
        self.assertTrue(all("tolerance" in x for x in found), found)

    def test_other_tolerance_forms(self) -> None:
        f = "tests/x.rs"
        for old, new in [
            ("const RTOL: f64 = 1e-5;", "const RTOL: f64 = 1e-4;"),
            ("close(ab, ba, 1e-10);", "close(ab, ba, 1e-2);"),
            (
                "assert!(e <= 2.0 * f64::EPSILON * y);",
                "assert!(e <= 4.0 * f64::EPSILON * y);",
            ),
        ]:
            with self.subTest(old):
                self.assertEqual(len(self.check({f: old + "\n"}, {f: new + "\n"})), 1)
                self.assertEqual(self.check({f: new + "\n"}, {f: old + "\n"}), [])

    def test_rows_are_not_tolerances(self) -> None:
        f = "tests/x.rs"
        self.assertEqual(self.check({f: "(1e-5, 2.0),\n"}, {f: "(1e-3, 2.0),\n"}), [])

    def test_tolerance_variable(self) -> None:
        f = "tests/x.rs"
        found = self.check(
            {f: "const RTOL: f64 = 1e-12;\n"}, {f: "const RTOL: f64 = DEFAULT_RTOL;\n"}
        )
        self.assertEqual(found, [f"{f}: tolerance is not a literal anymore"])

    def test_src_tests_only(self) -> None:
        code = "fn f() {{ let rtol = {}; }}\n\n#[cfg(test)]\nmod tests {{\n    fn t() {{ let rtol = {}; }}\n}}\n"
        found = self.check(
            {"src/a.rs": code.format("1e-14", "1e-14")},
            {"src/a.rs": code.format("1e-3", "1e-13")},
        )
        self.assertEqual(found, ["src/a.rs: tolerance 1e-13 > 1e-14"])

    def test_ignored_and_removed(self) -> None:
        f = "tests/x.rs"
        before = TEST.format(rtol="1e-14", tol="1e-13")
        after = before.replace(
            "#[test]\nfn test_g() {}", "#[test]\n#[ignore]\nfn test_g() {}"
        )
        self.assertEqual(self.check({f: before}, {f: after}), [f"{f}: ignored test"])
        after = before.replace("\n#[test]\nfn test_g() {}\n", "")
        self.assertEqual(
            self.check({f: before}, {f: after}), [f"{f}: removed test `test_g`"]
        )
        # a rename and a move to another file
        after = before.replace("test_g", "test_h")
        self.assertEqual(self.check({f: before}, {f: after}), [])
        gone = before.replace("\n#[test]\nfn test_g() {}\n", "")
        moved: Files = {
            f: gone,
            "tests/y.rs": "#[test]\nfn test_g() {}\n",
        }
        self.assertEqual(self.check({f: before}, moved), [])

    def test_rows(self) -> None:
        f = "tests/x.rs"
        before = TEST.format(rtol="1e-14", tol="1e-13")
        self.assertEqual(
            self.check({f: before}, {f: before.replace("        (3.0, 4.0),\n", "")}),
            [f"{f}: removed 1 row of a test table"],
        )
        # a new output for the same input
        self.assertEqual(
            self.check({f: before}, {f: before.replace("(3.0, 4.0)", "(3.0, 5.0)")}), []
        )

    def test_exclusions(self) -> None:
        f = "tests/properties/ranges.rs"
        before = "property!(no_panic, f(a in positive()));\n"
        after = "property!(no_panic, f(a in positive()) if a < 1e300);\n"
        self.assertEqual(
            self.check({f: before}, {f: after}), [f"{f}: property-test exclusion"]
        )
        # a guard of a property with default strategies
        guard = "property!(no_panic, f(x) if x > 0.0);\n"
        self.assertEqual(
            self.check({f: before}, {f: before + guard}),
            [f"{f}: property-test exclusion"],
        )
        # a match arm is not a guard
        arm = "    Some(w) if w > 0 => 1,\n"
        self.assertEqual(self.check({f: before}, {f: before + arm}), [])
        # a new file is a new test, which is not weakened, but its ignored tests are
        self.assertEqual(self.check({"x": ""}, {f: after}), [])
        ignored = "#[test]\n#[ignore]\nfn t() {}\n"
        self.assertEqual(
            self.check({"x": ""}, {"tests/new.rs": ignored}),
            ["tests/new.rs: ignored test"],
        )

    def test_removed_attributes(self) -> None:
        f = "tests/x.rs"
        before = '// c\n#[test]\n#[should_panic(expected = "x")]\n/// doc\nfn t() {}\n\nfn a() {}\n'
        self.assertEqual(
            self.check({f: before}, {f: "fn a() {}\n"}), [f"{f}: removed test `t`"]
        )
        # the function stays
        kept = before.replace("#[test]\n", "")
        self.assertEqual(
            self.check({f: before}, {f: kept}),
            [f"{f}: removed `#[test]` of a function"],
        )

    def test_unicode(self) -> None:
        f = "tests/\u00e9.rs"
        line = 'let x = "\u2028"; let rtol = {};\n'
        found = self.check({f: line.format("1e-14")}, {f: line.format("1e-13")})
        self.assertEqual(found, [f"{f}: tolerance 1e-13 > 1e-14"])

    def test_xsref(self) -> None:
        f = "tests/data/xsref/f/d-d.bin"
        rows = [(1.0, 2.0, 1e-15, 0.0), (2.0, 3.0, 1e-15, 0.0)]
        loose = [(1.0, 2.0, 1e-14, 0.0), (2.0, 3.0, 1e-15, 0.0)]
        scaled = [(1.0, 2.0, 1e-15, 1.0), (2.0, 3.0, 1e-15, 0.0)]
        self.assertEqual(self.check({f: table(rows)}, {f: table(rows[::-1])}), [])
        self.assertEqual(
            self.check({f: table(rows)}, {f: table(loose)}),
            [f"{f}: larger tolerance or scale in 1 xsref row"],
        )
        self.assertEqual(
            self.check({f: table(rows)}, {f: table(scaled)}),
            [f"{f}: larger tolerance or scale in 1 xsref row"],
        )
        self.assertEqual(
            self.check({f: table(rows)}, {f: table(rows[:1])}),
            [f"{f}: removed 1 row of an xsref table"],
        )
        self.assertEqual(
            self.check({f: table(rows)}, {f: None}), [f"{f}: removed xsref table"]
        )
        # a scale column that is added is a scale
        plain = table([r[:3] for r in rows], scaled=False)
        self.assertEqual(
            self.check({f: plain}, {f: table([(1.0, 2.0, 1e-15, 1e3), rows[1]])}),
            [f"{f}: larger tolerance or scale in 1 xsref row"],
        )
        self.assertEqual(self.check({f: plain}, {f: table(rows)}), [])
        # only another table replaces a table
        other = {"tests/data/xsref/f/README.md": "x"}
        self.assertEqual(
            self.check({f: table(rows)}, {f: None, **other}),
            [f"{f}: removed xsref table"],
        )
        replaced = {"tests/data/xsref/f/d-d_d.bin": table(rows)}
        self.assertEqual(self.check({f: table(rows)}, {f: None, **replaced}), [])
        # not another, different table
        different = {"tests/data/xsref/f/d-d_d.bin": table(rows[:1])}
        self.assertEqual(
            self.check({f: table(rows)}, {f: None, **different}),
            [f"{f}: removed xsref table"],
        )

    def test_unreadable_head_table(self) -> None:
        # a changed table whose blob can't be read is an error, not a removal
        repo = Repo({"tests/data/xsref/f/d-d.bin": table([(1.0, 2.0, 1e-15, 0.0)])})
        self.addCleanup(repo.tmp.cleanup)
        repo.commit({"tests/data/xsref/f/d-d.bin": table([(1.0, 2.0, 1e-14, 0.0)])})
        cwd = Path.cwd()
        os.chdir(repo.root)
        try:
            names = [["M", "tests/data/xsref/f/d-d.bin"]]
            with self.assertRaises(RuntimeError):
                c.check_tables("HEAD~", "no-such-rev", names)
        finally:
            os.chdir(cwd)

    def test_qualified_property_check(self) -> None:
        f = "tests/properties/x.rs"
        before = "property!(checks::no_panic, f(x in positive()));\n"
        guard = "property!(checks::no_panic, f(x in positive()) if x > 1.0);\n"
        self.assertEqual(
            self.check({f: before}, {f: before + guard}),
            [f"{f}: property-test exclusion"],
        )

    def test_unicode_line_separator_in_src(self) -> None:
        f = "src/x.rs"
        before = 'const S: &str = "\u2028";\n#[cfg(test)]\nmod t {\n    const RTOL: f64 = 1e-5;\n}\n'
        found = self.check({f: before}, {f: before.replace("1e-5", "1e-4")})
        self.assertEqual(found, [f"{f}: tolerance 0.0001 > 1e-05"])


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
