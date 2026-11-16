"""Plot where `hyppfq` returns `None`, and its accuracy elsewhere, for p <= 4 and q <= 3.

Usage (from the repository root; or `just pfq-map`):

    uv run --project tools tools/accuracy/pfq_map.py [--no-cache] [--jobs N] [--per-bin N]

Writes `docs/img/hyppfq-error.svg` and `docs/img/hyppfq-none.svg`, with a panel per p and q.
Each panel bins random inputs by z (on a signed log scale) and by the scale s of the parameters,
which are uniform in [-s, s]. Their `-int` variants are the same for p >= 2, but with
a_2 = a_1 + k, for a random integer k in [-3, 3]. The references are those of `specs.py`, cached
as in `accuracy.py`.
"""

import argparse
import itertools
import math
import random
from pathlib import Path
from typing import Final

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import specs
from matplotlib import colormaps
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from matplotlib.patches import Patch

import accuracy

# the spec (and so the reference) of each p and q
SPECS: Final = {
    (0, 0): "hyp0f0",
    (1, 0): "hyp1f0",
    (0, 1): "hyp0f1",
    (1, 1): "hyp1f1",
    (2, 0): "hyppfq#2f0",
    (2, 1): "hyp2f1",
    (0, 2): "hyppfq#0f2",
    (1, 2): "hyp1f2",
    (2, 2): "hyp2f2",
    (3, 2): "hyp3f2",
    **{(p, 3): f"hyppfq#{p}f3" for p in range(5)},
}
P_MAX: Final = 4
Q_MAX: Final = 3

# |z| in [10^Z_MIN, 10^Z_MAX], in Z_BINS bins on each side of 0; x = sign(z) (log10|z| - Z_MIN)
Z_MIN: Final = -3
Z_MAX: Final = 6
Z_BINS: Final = 27
X_MAX: Final = Z_MAX - Z_MIN
# s in [10^S_MIN, 10^S_MAX], in S_BINS bins
S_MIN: Final = -2
S_MAX: Final = 2
S_BINS: Final = 12

# a_2 = a_1 + k in the maps of a_i that differ by integers, with |k| <= K_MAX
K_MAX: Final = 3

# the domain index in the cache file names of the references, without and with a_2 = a_1 + k
MAP_DOMAINS: Final = {False: 100, True: 101}
EPS: Final = 2.0**-52
# the color levels of the error, in eps; at 1e16 no correct digits are left
LEVELS: Final = [1, 4, 16, 64, 256, 1e3, 1e4, 1e6, 1e8, 1e16]
NONE_COLOR: Final = "0.8"
NOREF_COLOR: Final = "tan"

type Grid = np.ndarray[tuple[int, int], np.dtype[np.float64]]


def inputs(
    p: int, q: int, per_bin: int, int_diff: bool
) -> list[tuple[int, int, accuracy.Args]]:
    """`per_bin` random inputs (a_1, ..., a_p, b_1, ..., b_q, z) of each bin (i, j) of z and
    s, except where z >= 1 for p = q + 1, and where z > 0 for 2F0, where it is complex. If
    `int_diff`, then a_2 = a_1 + k for a random integer k in [-K_MAX, K_MAX]."""
    rng = random.Random(f"{p}F{q}-int" if int_diff else f"{p}F{q}")
    out: list[tuple[int, int, accuracy.Args]] = []
    dx, ds = X_MAX / Z_BINS, (S_MAX - S_MIN) / S_BINS
    for i in range(2 * Z_BINS):
        x0 = -X_MAX + i * dx
        if (p == q + 1 and x0 >= -Z_MIN) or ((p, q) == (2, 0) and x0 >= 0):
            continue
        for j in range(S_BINS):
            for _ in range(per_bin):
                x = rng.uniform(x0, x0 + dx)
                z = math.copysign(10 ** (abs(x) + Z_MIN), x)
                s = 10 ** rng.uniform(S_MIN + j * ds, S_MIN + (j + 1) * ds)
                ab = [rng.uniform(-s, s) for _ in range(p + q)]
                if int_diff:
                    # (a_1 on a grid of 2^-46, so that a_1 + k is exact for |a_1 + k| < 2^7)
                    ab[0] = round(ab[0] * 2**46) / 2**46
                    ab[1] = ab[0] + rng.randint(-K_MAX, K_MAX)
                out.append((i, j, (*ab, z)))
    return out


def grids(
    pts: list[tuple[int, int, accuracy.Args]],
    ys: list[complex],
    refs: list[accuracy.Ref | str],
) -> tuple[Grid, Grid]:
    """The geometric mean of the errors in eps of each bin, as in `accuracy.stats` (NaN where it
    has no inputs, all are `None`, or the others have no reference or an infinite error), and the share of `None` (NaN where it has no inputs)."""
    errs: dict[tuple[int, int], list[float]] = {}
    nones: dict[tuple[int, int], list[bool]] = {}
    for (i, j, _), y, ref in zip(pts, ys, refs, strict=True):
        # (the dump program writes NaN for `None`, and `hyppfq` is only NaN where complex)
        nones.setdefault((i, j), []).append(math.isnan(y.real))
        if not math.isnan(y.real) and not isinstance(ref, str):
            errs.setdefault((i, j), []).append(accuracy.error(y, ref))
    err = np.full((S_BINS, 2 * Z_BINS), np.nan)
    none = np.full((S_BINS, 2 * Z_BINS), np.nan)
    for (i, j), es in errs.items():
        stats = accuracy.stats(es)
        if stats and stats[1] is not None:
            err[j, i] = stats[1] / EPS
    for (i, j), ns in nones.items():
        none[j, i] = sum(ns) / len(ns)
    return err, none


def plot(
    panels: dict[tuple[int, int], tuple[Grid, Grid]],
    ps: range,
    none: bool,
    title: str,
    path: Path,
) -> None:
    """A panel per p in `ps` (rows) and q (columns): the share of `None`, or the gmean of the error.
    A bin without inputs is blank: there it is complex."""
    fig, axes = plt.subplots(
        len(ps),
        Q_MAX + 1,
        figsize=(12, 2.4 * len(ps)),
        sharex=True,
        sharey=True,
        layout="constrained",
    )
    keys = [(NONE_COLOR, "None"), ("white", "complex")]
    if none:
        fig.suptitle(f"hyppfq: share of None{title}")
        norm, cmap, label = Normalize(0, 1), colormaps["Reds"], "share of None"
    else:
        fig.suptitle(f"hyppfq: gmean of the relative error{title}")
        norm = BoundaryNorm(LEVELS, len(LEVELS) - 1)
        cmap, label = colormaps["viridis"].resampled(len(LEVELS) - 1), "eps"
        keys.insert(1, (NOREF_COLOR, "no reference"))
    kw = {"extent": (-X_MAX, X_MAX, S_MIN, S_MAX), "origin": "lower", "aspect": "auto"}
    im = None
    for (row, q), ax in np.ndenumerate(axes):
        p = ps[row]
        ax.set_title(f"{p}F{q}")
        if (p, q) not in panels:
            # (unless it is a polynomial)
            ax.set_facecolor(NONE_COLOR)
            ax.text(0, (S_MIN + S_MAX) / 2, "None", ha="center", va="center")
            continue
        err, share = panels[p, q]
        if none:
            im = ax.imshow(share, cmap=cmap, norm=norm, **kw)
        else:
            for color, mask in [
                (NONE_COLOR, share == 1),
                (NOREF_COLOR, np.isnan(err) & (share < 1)),
            ]:
                ax.imshow(
                    np.where(mask, 1.0, np.nan), cmap=ListedColormap([color]), **kw
                )
            # (saturated at the outer levels)
            im = ax.imshow(
                np.clip(err, LEVELS[0], LEVELS[-1]), cmap=cmap, norm=norm, **kw
            )
    ticks = [-X_MAX, Z_MIN, 0, -Z_MIN, X_MAX]
    labels = [
        f"$-10^{{{Z_MAX}}}$",
        "$-1$",
        f"$\\pm 10^{{{Z_MIN}}}$",
        "$1$",
        f"$10^{{{Z_MAX}}}$",
    ]
    for ax in axes[-1]:
        ax.set_xticks(ticks, labels)
        ax.set_xlabel("z")
    for ax in axes[:, 0]:
        ax.set_yticks(
            range(S_MIN, S_MAX + 1), [f"$10^{{{k}}}$" for k in range(S_MIN, S_MAX + 1)]
        )
        ax.set_ylabel("s")
    assert im is not None
    bar = fig.colorbar(im, ax=axes, label=label, shrink=0.5)
    if not none:
        bar.set_ticks(
            LEVELS,
            labels=[
                f"{v:g}" if v < 1e4 else f"$10^{{{math.log10(v):.0f}}}$" for v in LEVELS
            ],
        )
    handles = [Patch(facecolor=c, edgecolor="0.5", label=t) for c, t in keys]
    fig.legend(handles=handles, loc="outside lower center", ncols=len(keys))
    fig.savefig(path, metadata={"Date": None})


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None
    )
    ap.add_argument("--no-cache", action="store_true", help="recompute all references")
    ap.add_argument("--jobs", type=int, default=8, help="reference worker processes")
    ap.add_argument("--per-bin", type=int, default=4, help="random inputs per bin")
    opts = ap.parse_args()
    mp.mp.dps = 30
    # the figure pairs: with a_2 = a_1 + k or not, the file and title suffixes, and the rows p
    figures = [
        (False, "", "", range(P_MAX + 1)),
        (True, "-int", ", $a_2 = a_1 + k$", range(2, P_MAX + 1)),
    ]
    pts = {
        (pq, int_diff): inputs(*pq, opts.per_bin, int_diff)
        for int_diff, *_, ps in figures
        for pq in SPECS
        if pq[0] in ps
    }
    lines = [(f"hyppfq{p}", args) for ((p, _), _), xs in pts.items() for *_, args in xs]
    ys = iter(accuracy.run([*accuracy.build_rust(), "xsf"], lines))
    jobs = [
        accuracy.Job(
            specs.BY_NAME[SPECS[pq]], MAP_DOMAINS[int_diff], [a for *_, a in xs]
        )
        for (pq, int_diff), xs in pts.items()
    ]
    refs = accuracy.compute_refs(jobs, not opts.no_cache, opts.jobs)
    panels = {
        key: grids(xs, list(itertools.islice(ys, len(xs))), rs)
        for (key, xs), rs in zip(pts.items(), refs, strict=True)
    }
    plt.rcParams["svg.hashsalt"] = "hyppfq"
    for int_diff, suffix, title, ps in figures:
        sub = {pq: v for (pq, i), v in panels.items() if i == int_diff}
        for none, name in [(False, "error"), (True, "none")]:
            path = accuracy.ROOT / "docs" / "img" / f"hyppfq-{name}{suffix}.svg"
            plot(sub, ps, none, title, path)


if __name__ == "__main__":
    main()
