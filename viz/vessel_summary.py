"""Figure made from metrics_slice.csv (+ vessel_regions.csv): the vessel summary.

Ports the current viz.plot_metrics_vs_slice, which main() saves as
vessel_summary.png. Five stacked panels against elevational position, every
one of them an estimation result:

    1  NRMSE over the valid ROI mask, raw vs 2-D SG
    2  NRMSE over the full lumen (no-slip fill vs SG)
    3  flow error
    4  bias
    5  estimator yield

Error bars are the across-seed standard deviation and appear only where
n_seeds > 1. Regions are still read, but only to shade and label the three
flow segments behind the panels; the geometry and truth-disturbance curves
they come from are not plotted -- this figure is about the estimate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import style as st

# A six-panel stack needs a smaller step than the single-panel figures.
FS_TICK = 16
FS_LABEL = 18
FS_LEGEND = 16
FS_ANNOT = 14

C_RAW = "#2966AD"       # [0.16 0.40 0.68]
C_SG = "#D13D26"        # [0.82 0.24 0.15]
C_GRAY = "#59616B"      # [0.35 0.38 0.42]
C_GOLD = "#DB8514"      # [0.86 0.52 0.08]

# Shading and labels come from style.segment: the pipeline's finer regions are
# collapsed onto upstream / stenosis / downstream.


def plot_vessel_summary(msl, regions=None, title=""):
    m = msl.sort_values("slice_pos_mm")
    y = m["slice_pos_mm"].to_numpy(float)
    reg = _normalise_regions(m, regions)
    n_seeds = _col(m, "n_seeds", default=1.0)

    fig = st.new_figure((12.5, 15.0))
    axes = fig.subplots(5, 1, sharex=True)

    filled = []          # per panel: did anything finite get drawn?

    # 1 NRMSE over the valid mask
    ax = axes[0]
    got = _series(ax, y, 100 * _col(m, "nrmse_raw"),
                  100 * _col(m, "nrmse_raw_seed_std"), n_seeds, C_RAW, "o", "Raw")
    got |= _series(ax, y, 100 * _col(m, "nrmse_sg"),
                   100 * _col(m, "nrmse_sg_seed_std"), n_seeds, C_SG, "s",
                   "After SG (2-D)")
    filled.append(got)
    _axis(ax, "NRMSE [%]")
    ax.legend(loc="upper left", fontsize=FS_LEGEND, frameon=False, ncols=2)
    # Extra pad leaves room for the region strip between title and axes.
    st.set_title(ax, f"{title}  vessel-region error summary".strip(), pad=38)

    # 2 NRMSE over the full lumen
    ax = axes[1]
    got = _series(ax, y, 100 * _col(m, "nrmse_fill_full"),
                  100 * _col(m, "nrmse_fill_full_seed_std"), n_seeds, C_GOLD,
                  "o", "No-slip fill")
    got |= _series(ax, y, 100 * _col(m, "nrmse_sg_full"),
                   100 * _col(m, "nrmse_sg_full_seed_std"), n_seeds, C_SG, "s",
                   "After SG (2-D)")
    filled.append(got)
    _axis(ax, "Full-lumen NRMSE [%]")
    ax.legend(loc="upper left", fontsize=FS_LEGEND, frameon=False, ncols=2)

    # 3 flow error
    ax = axes[2]
    ax.axhline(0, color="#262626", lw=1.4, zorder=2)
    got = _series(ax, y, _col(m, "flow_error_fill_pct"),
                  _col(m, "flow_error_fill_pct_seed_std"), n_seeds, C_GOLD,
                  "o", None)
    got |= _series(ax, y, _col(m, "flow_error_sg_pct"),
                   _col(m, "flow_error_sg_pct_seed_std"), n_seeds, C_SG, "s",
                   None)
    filled.append(got)
    _axis(ax, "Flow error [%]")

    # 4 bias
    ax = axes[3]
    ax.axhline(0, color="#262626", lw=1.4, zorder=2)
    got = _series(ax, y, _col(m, "bias_raw_mms"), _col(m, "bias_raw_seed_std_mms"),
                  n_seeds, C_RAW, "o", None)
    got |= _series(ax, y, _col(m, "bias_sg_mms"), _col(m, "bias_sg_seed_std_mms"),
                   n_seeds, C_SG, "s", None)
    filled.append(got)
    _axis(ax, "Bias [mm/s]")

    # 5 estimator yield
    ax = axes[4]
    filled.append(_series(ax, y, 100 * _col(m, "valid_frac"),
                          100 * _col(m, "valid_frac_seed_std"), n_seeds,
                          C_GRAY, "o", None))
    _annotate_seed_counts(ax, y, n_seeds)
    ax.set_ylim(0, 112)
    _axis(ax, "Valid ROI [%]")
    ax.set_xlabel("Slice position [mm]", fontsize=FS_LABEL, labelpad=10)

    runs = _region_runs(reg)
    for ax, got in zip(axes, filled):
        _shade_regions(ax, runs)
        _note_if_empty(ax, got)
    _label_regions(axes[0], runs)
    axes[0].set_xlim(runs[0][0], runs[-1][1])
    return fig


def _col(T, name, default=np.nan):
    """Column as a float array, or a constant when the table predates it."""
    if name in T.columns:
        return T[name].to_numpy(float)
    return np.full(len(T), default, float)


def _series(ax, x, v, spread, n_seeds, colour, marker, label):
    """Draw one metric; returns whether it had anything finite to show."""
    ax.plot(x, v, f"-{marker}", color=colour, lw=2.0, ms=7, label=label,
            zorder=3)
    ok = (n_seeds > 1) & np.isfinite(v) & np.isfinite(spread)
    if ok.any():
        ax.errorbar(x[ok], v[ok], spread[ok], fmt="none", ecolor=colour,
                    elinewidth=1.6, capsize=7, zorder=4)
    return bool(np.isfinite(v).any())


def _axis(ax, ylabel):
    st.style_axes(ax)
    ax.tick_params(labelsize=FS_TICK)
    ax.set_ylabel(ylabel, fontsize=FS_LABEL, labelpad=10)


def _normalise_regions(msl, regions):
    """Position and segment label per slice, sorted; only what shading needs.

    Falls back to the columns carried on msl itself, and finally to a single
    unlabelled band, so the figure still draws for a vessel the pipeline never
    classified.
    """
    need = ["slice_pos_mm", "vessel_region"]
    if regions is None or len(regions) == 0:
        if all(c in msl.columns for c in need):
            regions = msl[need].drop_duplicates()
        else:
            regions = pd.DataFrame({"slice_pos_mm": msl["slice_pos_mm"],
                                    "vessel_region": "uniform vessel"})
    return regions.sort_values("slice_pos_mm").copy()


def _region_runs(reg):
    """[(left, right, colour, label)] per contiguous run of one segment.

    Runs are formed after collapsing the pipeline's regions onto the three
    flow segments, so neighbouring regions that share a segment shade as one
    band instead of a stripe per label.

    Edges sit halfway between neighbouring slices, so the shading tiles the
    axis without gaps or overlap.
    """
    y = reg["slice_pos_mm"].to_numpy(float)
    if y.size == 1:
        edges = np.array([y[0] - 0.5, y[0] + 0.5])
    else:
        mid = (y[:-1] + y[1:]) / 2
        edges = np.concatenate([[2 * y[0] - mid[0]], mid, [2 * y[-1] - mid[-1]]])

    seg = np.array([st.segment(v)[0] for v in reg["vessel_region"]], dtype=object)
    change = np.concatenate([[True], seg[1:] != seg[:-1]])
    starts = np.flatnonzero(change)
    stops = np.append(starts[1:] - 1, len(seg) - 1)
    out = []
    for s, e in zip(starts, stops):
        _, colour, label = st.segment(reg["vessel_region"].iloc[s])
        out.append((edges[s], edges[e + 1], colour, label))
    return out


def _annotate_seed_counts(ax, y, n_seeds):
    """How many repeats each slice had.

    Stated once when every slice ran the same number of seeds -- per-slice
    labels collide into an unreadable strip past a handful of slices -- and
    per slice only for the ones that differ from that norm.
    """
    rep = n_seeds[np.isfinite(n_seeds)]
    if rep.size == 0 or rep.max() <= 1:
        return
    common = np.bincount(rep.astype(int)).argmax()
    ax.text(0.995, 0.04, f"n = {common:.0f} seeds per slice",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=FS_ANNOT, color="0.25")
    for k in np.flatnonzero(n_seeds != common):
        ax.text(y[k], 103, f"n={n_seeds[k]:.0f}", ha="center", va="bottom",
                fontsize=FS_ANNOT, color="0.2")


def _note_if_empty(ax, has_data):
    """Say so when a panel's metric is missing, instead of showing a void.

    Results written before a metric existed simply have no column, and an
    empty panel reads as "the value is zero" rather than "not recorded".
    """
    if has_data:
        return
    ax.text(0.5, 0.5, "not in this results file", transform=ax.transAxes,
            ha="center", va="center", fontsize=FS_ANNOT, style="italic",
            color="0.45")


def _shade_regions(ax, runs):
    for left, right, colour, _ in runs:
        ax.axvspan(left, right, color=colour, alpha=0.10, lw=0, zorder=0)


def _label_regions(ax, runs):
    """Region names in the strip above the top panel.

    Above rather than inside: inside they collide with the legend, and the
    panel is where the data has to be readable. Runs too narrow to hold their
    own name are left unlabelled -- the shading still identifies them -- and
    the outermost labels are pulled inside the axes so they cannot be clipped.
    """
    lo, hi = runs[0][0], runs[-1][1]
    span = hi - lo
    for left, right, _, label in runs:
        if (right - left) < 0.06 * span:
            continue
        x = min(max((left + right) / 2, lo + 0.07 * span), hi - 0.07 * span)
        ax.text(x, 1.012, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=FS_ANNOT, fontweight="bold",
                color="#2E2E2E", clip_on=False, zorder=5)
