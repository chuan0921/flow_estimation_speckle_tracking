"""Shared look for every figure in the pipeline.

One place to change type sizes, palette and spacing so all figures stay a
set. Sizes are picked so text stays readable when a figure is dropped into a
slide at half width -- the axis label is about 3% of the figure width.

The palette is carried over from the MATLAB +viz figures on purpose, so a
ported figure can be compared against its predecessor without the colour
changing under you.
"""

from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- type sizes (pt) -------------------------------------------------------
FS_TICK = 19
FS_LABEL = 22
FS_TITLE = 24
FS_LEGEND = 18
FS_ANNOT = 17

# ---- palette (same values as the MATLAB figures) ---------------------------
C_RAW = "#3373BF"      # [0.20 0.45 0.75]  raw / estimate
C_SG = "#D9541A"       # [0.85 0.33 0.10]  Savitzky-Golay
C_TRUTH = "#000000"    #                   ground truth, dashed
C_POINT = "#B3BFD1"    # [0.70 0.75 0.82]  faint per-slice scatter
C_BAR = "#99A3B3"      # [0.60 0.64 0.70]  count bars
C_POINT2 = "#9EA9B8"   # [0.62 0.66 0.72]  per-slice dots in the summaries
C_LOA = "#BF3333"      # [0.75 0.20 0.20]  limits of agreement
C_LINE = "#262626"     # [0.15 0.15 0.15]  arrows, bias line
C_GRID = "#D9D9D9"
C_MISS = "#B8B8B8"     # points with no answer

# matplotlib has no parula; viridis is the closest perceptual equivalent.
CMAP_SEQ = "viridis"
CMAP_DIV = "bwr"

BAND_ALPHA = 0.18
DPI = 150

# ---- default figure sizes (inches) -----------------------------------------
FIG_SINGLE = (11.0, 7.5)
FIG_STACKED = (11.0, 9.0)
FIG_TALL = (11.0, 11.0)
FIG_SQUARE = (9.5, 9.5)


def new_figure(size=FIG_SINGLE):
    """Blank white figure with constrained layout (labels never clip)."""
    return plt.figure(figsize=size, facecolor="w", layout="constrained")


def style_axes(ax, grid=True):
    """Tick sizes, spine weight and a grid that stays behind the data."""
    ax.tick_params(labelsize=FS_TICK, direction="out", length=6, width=1.2, pad=6)
    for s in ax.spines.values():
        s.set_linewidth(1.2)
    if grid:
        ax.grid(True, color=C_GRID, linewidth=0.9)
        ax.set_axisbelow(True)
    return ax


def label_axes(ax, xlabel=None, ylabel=None, title=None):
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=FS_LABEL, labelpad=10)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=FS_LABEL, labelpad=10)
    if title is not None:
        set_title(ax, title)
    return ax


def set_title(ax, text, fontsize=FS_TITLE, min_fontsize=13, pad=14):
    """Bold title, shrunk just enough to fit the figure width.

    Dataset names are long and arbitrary, so a fixed size would eventually
    run off the canvas -- measure the rendered text and step down instead.
    """
    t = ax.set_title(text, fontsize=fontsize, fontweight="bold", pad=pad)
    fig = ax.figure
    fig.canvas.draw()
    while fontsize > min_fontsize:
        # The title is centred on the axes, not on the figure, so a width
        # test is not enough -- compare both ends against the canvas.
        bb = t.get_window_extent(fig.canvas.get_renderer())
        fb = fig.get_window_extent()
        if bb.x0 >= fb.x0 + 2 and bb.x1 <= fb.x1 - 2:
            break
        fontsize -= 1
        t.set_fontsize(fontsize)
        fig.canvas.draw()
    return t


def legend(ax, loc="lower left", handles=None, labels=None):
    args = () if handles is None else (handles, labels)
    return ax.legend(*args, loc=loc, fontsize=FS_LEGEND, frameon=False,
                     handlelength=2.4, borderaxespad=0.6, labelspacing=0.5)


def band(ax, x, m, s, color, alpha=BAND_ALPHA, zorder=1):
    """+/-1 std ribbon, skipped where the mean or the spread is undefined."""
    x = np.asarray(x, float)
    m = np.asarray(m, float)
    s = np.asarray(s, float)
    ok = np.isfinite(m) & np.isfinite(s)
    if not ok.any():
        return None
    return ax.fill_between(x[ok], (m - s)[ok], (m + s)[ok], color=color,
                           alpha=alpha, linewidth=0, zorder=zorder)


def collapse(r, v):
    """Across-slice mean, std and count at each radius.

    Mirrors the `collapse` helper the MATLAB figures used: radii are matched
    at 1e-6 mm, NaNs are skipped, and a radius no slice could measure comes
    back as NaN rather than 0. std uses the N-1 normalisation and is 0 -- not
    NaN -- where a single slice contributed, which is what MATLAB's
    std(x,'omitnan') returns and what keeps the ribbon closed there.

    Returns (radius, mean, std, n_finite), all sorted by radius.
    """
    d = pd.DataFrame({"r": np.round(np.asarray(r, float), 6),
                      "v": np.asarray(v, float)})
    g = d.groupby("r", sort=True)["v"]
    n = g.count()
    m = g.mean()
    s = g.std(ddof=1).where(n != 1, 0.0)
    return (n.index.to_numpy(), m.to_numpy(), s.to_numpy(),
            n.to_numpy().astype(float))


def set_suptitle(fig, text, fontsize=FS_TITLE, min_fontsize=13):
    """Figure-level title, shrunk just enough to fit the canvas width."""
    t = fig.suptitle(text, fontsize=fontsize, fontweight="bold")
    fig.canvas.draw()
    while fontsize > min_fontsize:
        bb = t.get_window_extent(fig.canvas.get_renderer())
        fb = fig.get_window_extent()
        if bb.x0 >= fb.x0 + 2 and bb.x1 <= fb.x1 - 2:
            break
        fontsize -= 1
        t.set_fontsize(fontsize)
        fig.canvas.draw()
    return t


def span(v):
    """Finite data range, widened to a non-degenerate interval."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (0.0, 1.0)
    lo, hi = float(v.min()), float(v.max())
    if hi <= lo:
        return (lo - 0.5, hi + 0.5)
    return (lo, hi)


# ---- vessel segments -------------------------------------------------------
# The pipeline classifies each slice from geometry and truth disturbance, which
# yields up to six labels; post-stenosis in particular can be a single slice.
# Figures collapse those onto the three flow states that actually differ:
# entrance flow that has not developed, the accelerated core through the
# throat, and the disturbed wake downstream. The boundary is geometric -- once
# the lumen starts widening again the slice is downstream.
SEGMENT_OF = {
    "pre-stenosis": "upstream",
    "narrowing": "upstream",
    "stenosis throat": "stenosis",
    "post-stenosis": "downstream",
    "disturbed flow": "downstream",
    "recovery": "downstream",
}
SEGMENT_STYLE = {
    "upstream": ("#B3BDC7", "Upstream (developing)"),
    "stenosis": ("#DB402E", "Stenosis (laminar)"),
    "downstream": ("#7A529E", "Downstream (turbulent)"),
}
SEGMENT_FALLBACK = ("#B8BCBF", "Uniform vessel")


def segment(region_name):
    """(key, colour, label) for a pipeline region label."""
    key = SEGMENT_OF.get(str(region_name).strip().lower())
    if key is None:
        return None, *SEGMENT_FALLBACK
    return key, *SEGMENT_STYLE[key]


def representative_seed(T):
    """Keep one realisation per (dataset, slice): the lowest seed_base.

    Mirrors representative_seed_points in main.m, so a Python figure and its
    MATLAB predecessor pick the same rows. Spatial figures need this -- every
    seed lands on the same lattice, so plotting them all just overprints the
    same coordinates. Tables without a seed_base column pass through
    unchanged.
    """
    if "seed_base" not in T.columns:
        return T
    keys = [c for c in ("dataset", "slice_idx") if c in T.columns]
    if not keys:
        return T[T["seed_base"] == T["seed_base"].min()]
    return T[T["seed_base"] == T.groupby(keys)["seed_base"].transform("min")]


def seed_note(T):
    """' (N slices x M seeds)' when a table carries repeats, else ''."""
    if "seed_base" not in T.columns or "slice_idx" not in T.columns:
        return ""
    per_slice = T.groupby("slice_idx")["seed_base"].nunique()
    if per_slice.max() <= 1:
        return ""
    return f" ({per_slice.size} slices x {per_slice.max()} seeds)"


def order_datasets(col):
    """Dataset names ordered by peak speed, with the speed in cm/s.

    The speed is the trailing number in the folder name (sim_v080cms -> 80).
    If any name lacks one the whole set falls back to alphabetical, so the
    ordering never mixes two rules. Shared by every summary figure -- the
    MATLAB version carried four copies of this.
    """
    names = list(pd.unique(pd.Series(col)))
    speed = []
    for nm in names:
        m = re.search(r"(\d+)\s*cms?$", str(nm))
        speed.append(float(m.group(1)) if m else np.nan)
    speed = np.array(speed, float)
    if np.all(np.isfinite(speed)):
        k = np.argsort(speed, kind="stable")
    else:
        k = np.argsort([str(n) for n in names], kind="stable")
    return [names[i] for i in k], speed[k]


def speed_labels(names, speed):
    """Axis labels for the ordered datasets: '80 cm/s', or the raw name."""
    return [f"{s:g} cm/s" if np.isfinite(s) else str(n)
            for n, s in zip(names, speed)]


def dataset_colors(n):
    """One distinguishable colour per dataset (MATLAB's `lines` stand-in)."""
    cyc = plt.get_cmap("tab10").colors
    return [cyc[i % len(cyc)] for i in range(max(n, 1))]


def save(fig, path):
    """Write the PNG and close the figure. Returns the path."""
    out_dir = os.path.dirname(os.path.abspath(path))
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(path, dpi=DPI, facecolor="w")
    plt.close(fig)
    return path
