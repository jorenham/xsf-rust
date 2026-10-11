"""Check that `Cargo.toml` excludes exactly the xsf headers that `build.rs` doesn't need.

The needed headers are those in `build.rs`, and (recursively) their `#include "..."`s.
Preprocessor conditions are ignored, so headers that are only needed on some platforms are kept.
"""

import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INCLUDE = ROOT / "xsf" / "include"


def find(name: str, *dirs: Path) -> Path | None:
    """Like `#include "name"`; returns `None` for system headers, e.g. "assert.h"."""
    return next(((d / name).resolve() for d in dirs if (d / name).is_file()), None)


def main() -> int:
    # e.g. "cephes/jv.h" in `WRAPPER_INCLUDES`, and "xsf/third_party/kokkos/mdspan.hpp"
    names = re.findall(r'"([\w/]+\.h(?:pp)?)"', (ROOT / "build.rs").read_text())
    todo = [find(name, INCLUDE / "xsf", INCLUDE) for name in names]
    needed = set()
    while todo:
        path = todo.pop()
        if path and path not in needed:
            needed.add(path)
            includes = re.findall(
                r'^\s*#\s*include\s*"([^"]+)"', path.read_text(), re.MULTILINE
            )
            todo += [find(name, path.parent, INCLUDE) for name in includes]

    unused = {
        f"!{path.relative_to(ROOT).as_posix()}"
        for path in INCLUDE.rglob("*")
        if path.is_file() and path not in needed
    }
    excluded = {
        p
        for p in tomllib.loads((ROOT / "Cargo.toml").read_text())["package"]["include"]
        if p.startswith("!xsf/include/")
    }

    if excluded == unused:
        print(f"ok: {len(needed)} headers included, {len(unused)} excluded")
        return 0

    print("`package.include` in Cargo.toml is out of date:")
    for p in sorted(unused - excluded):
        print(f'  add:    "{p}",')
    for p in sorted(excluded - unused):
        print(f'  remove: "{p}",')
    return 1


if __name__ == "__main__":
    sys.exit(main())
