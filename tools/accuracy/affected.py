"""The specs of `specs.py` whose Rust sources depend on changed files of `src/`.

A coarse, conservative approximation: the file of a function depends on the files whose
top-level public items it names (other than as a method or field after a `.`), and a spec is
affected if its function is defined or re-exported in a changed file or in a file that it
depends on (transitively). A spec whose function is not defined by name (e.g. `pmv`) depends
on every file, and a changed file of `src/` that is not a Rust source of the revision (e.g. a
deleted one), or a changed module declaration or re-export, affects every spec.
"""

import io
import re
import subprocess
import tarfile
from pathlib import Path

import specs

import accuracy

SRC = accuracy.ROOT / "src"
ITEM = r"(?:const\s+|unsafe\s+)*(?:fn|const|static|struct|enum|trait|type)\s+(\w+)"
PUB = r"\bpub(?:\([^)]*\))?\s+"
REEXPORT = re.compile(rf"{PUB}use\b[^;]*;")
MOD = re.compile(rf"(?:{PUB})?mod\s+\w+\s*;")
# a name, also of a function passed as a value, that is not a method or field
USE = re.compile(r"(?<![.\w'])([A-Za-z_]\w*)")


def sources(rev: str | None) -> dict[Path, str]:
    """The Rust sources of `src/` in the working tree (`rev` None), or in a commit, without
    comments."""
    if rev is None:
        return {
            p: re.sub(r"//.*", "", p.read_text()) for p in sorted(SRC.rglob("*.rs"))
        }
    archive = ["git", "-C", str(accuracy.ROOT), "archive", rev, "--", "src"]
    try:
        tar = subprocess.check_output(archive, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        # (a revision without `src/`)
        return {}
    texts = {}
    with tarfile.open(fileobj=io.BytesIO(tar)) as t:
        for m in t.getmembers():
            # (None if not a file)
            if (f := t.extractfile(m)) and m.name.endswith(".rs"):
                texts[accuracy.ROOT / m.name] = re.sub(r"//.*", "", f.read().decode())
    return texts


def bindings(text: str) -> list[str]:
    """The module declarations and re-exports, without whitespace."""
    return [
        re.sub(r"\s", "", m[0])
        for m in re.finditer(f"{MOD.pattern}|{REEXPORT.pattern}", text)
    ]


def graph(
    texts: dict[Path, str],
) -> tuple[dict[str, set[Path]], dict[Path, set[Path]]]:
    """The files that define or re-export each name, and the files that each file depends on."""
    defs: dict[str, set[Path]] = {}
    exports: dict[str, set[Path]] = {}
    for p, t in texts.items():
        # (`#[doc(alias = "...")]` names the functions that have no Rust name, like `pmv`)
        for m in re.finditer(r'\bfn\s+(\w+)|\balias\s*=\s*"(\w+)"', t):
            defs.setdefault(m[1] or m[2], set()).add(p)
        for m in REEXPORT.finditer(t):
            for name in re.findall(r"\w+", m[0]):
                defs.setdefault(name, set()).add(p)
        for m in re.finditer(rf"^{PUB}{ITEM}", t, re.MULTILINE):
            exports.setdefault(m[1], set()).add(p)
    deps = {}
    for p, t in texts.items():
        # (a module or its re-exports don't use what they name, or a module would depend on all
        # its items)
        used = set(USE.findall(MOD.sub("", REEXPORT.sub("", t))))
        deps[p] = {q for n in used for q in exports.get(n, ())} - {p}
    return defs, deps


def affected(changed: list[Path], base: str, head: str | None = None) -> list[str]:
    """The names of the functions of the specs that the changed files of `head` can affect."""
    texts = sources(head)
    defs, deps = graph(texts)
    funcs = list(dict.fromkeys(s.func for s in specs.SPECS))
    changed_set = set(changed)
    if changed_set - set(deps):
        return funcs
    # (the graph has no edges from the users of a name to its re-exports)
    old = sources(base)
    if any(bindings(old.get(p, "")) != bindings(texts[p]) for p in changed_set):
        return funcs
    names: list[str] = []
    for name in funcs:
        roots = defs.get(re.split(r"[.\[]", name)[0], set(deps))
        seen, todo = set(roots), list(roots)
        while todo:
            for q in deps[todo.pop()] - seen:
                seen.add(q)
                todo.append(q)
        if seen & changed_set:
            names.append(name)
    return names
