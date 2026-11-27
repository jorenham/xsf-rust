"""The specs of `specs.py` whose Rust sources depend on changed files of `src/`.

A coarse, conservative approximation: the file of a function depends on the files whose
top-level public items it uses by name (a call that is not a method call, or a capitalized
name), and a spec is affected if its function is defined in a changed file or in a file that
it depends on (transitively). A spec whose function is not defined by name (e.g. `pmv`) depends
on every file. Files of only `mod`/`pub use` lines (re-exports) are skipped.
"""

import re
from pathlib import Path

import specs

import accuracy

SRC = accuracy.ROOT / "src"
ITEM = r"(?:const\s+|unsafe\s+)*(?:fn|const|static|struct|enum|trait|type)\s+(\w+)"
# a call (also of a generic function) that is not a method call, or a type or constant
USE = re.compile(r"(?<![.\w])(\w+)\b(?=\s*(?:::<[^>]*>)?\()|\b([A-Z]\w*)\b")


def graph() -> tuple[dict[str, set[Path]], dict[Path, set[Path]]]:
    """The files that define each function name, and the files that each file depends on."""
    texts = {p: p.read_text() for p in sorted(SRC.rglob("*.rs"))}
    # (a re-export hub has no function)
    texts = {
        p: t
        for p, t in texts.items()
        if re.search(r"\bfn\b", t) or p.name not in {"mod.rs", "lib.rs"}
    }
    defs: dict[str, set[Path]] = {}
    exports: dict[str, set[Path]] = {}
    for p, t in texts.items():
        # (`#[doc(alias = "...")]` names the functions that have no Rust name, like `pmv`)
        for m in re.finditer(r'\bfn\s+(\w+)|\balias\s*=\s*"(\w+)"', t):
            defs.setdefault(m[1] or m[2], set()).add(p)
        for m in re.finditer(rf"^pub(?:\(\w+\))?\s+{ITEM}", t, re.MULTILINE):
            exports.setdefault(m[1], set()).add(p)
    deps = {}
    for p, t in texts.items():
        used = {a or b for a, b in USE.findall(re.sub(r"//.*", "", t))}
        deps[p] = {q for n in used for q in exports.get(n, ())} - {p}
    return defs, deps


def affected(changed: list[Path]) -> list[str]:
    """The names of the functions of the specs that the changed files can affect."""
    defs, deps = graph()
    changed_set = set(changed)
    names: list[str] = []
    for name in dict.fromkeys(s.func for s in specs.SPECS):
        roots = defs.get(re.split(r"[.\[]", name)[0], set(deps))
        seen, todo = set(roots), list(roots)
        while todo:
            for q in deps[todo.pop()] - seen:
                seen.add(q)
                todo.append(q)
        if seen & changed_set:
            names.append(name)
    return names
