"""Figures made from several metrics_slice.csv files: one metric, every speed.

The per-dataset vessel_summary stacks five metrics for one flow condition. This
module transposes that: one figure per metric, with every speed condition drawn
in it, so a metric can be read against speed instead of against its neighbours.

Two of the panels are in physical units and carry the truth as well as the
estimate (velocity, flow rate); the errors come both in cm/s and normalised.
Both are here on purpose: NRMSE divides by that condition's own mean truth, so
it compares a method against itself but not two geometries against each other --
a stenosis and a straight pipe have mean velocities 25% apart, which moves NRMSE
without anything measuring better. Read RMSE across conditions, NRMSE within
one. Each returns its
own figure, so nothing shares a canvas -- a metric that needs a closer look can
be opened on its own.

    P = collect_speeds(['results/sim_v020cms_5seed/sim_v020cms', ...])
    for key in metric_keys():
        style.save(plot_speed_metric(P, key), f'{key}.png')

Speed is an ordered quantity, not an identity, so the four conditions are drawn
on a single-hue lightness ramp (light = slow) rather than on categorical hues.
Colour is the speed condition and nothing else. The raw and the Savitzky-Golay
result get their own figures rather than sharing one by line style: with four
conditions on the canvas a second encoding put eight curves in front of the
reader, and the two estimators are read against each other far less often than
each is read against speed. The two figures in physical units keep the truth
alongside the estimate, since that comparison is the point of them.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from . import style as st

# One hue per speed, in speed order. Four overlaid lines have to be told apart
# pair by pair, which a single-hue ramp cannot manage at four steps, so identity
# beats order here. Validated against a white surface on all pairs, not just
# adjacent ones: lightness band, chroma floor, CVD separation (worst pair dE 9.1
# protan / 8.6 tritan), normal-vision floor (16.0) and 3:1 contrast all pass.
# Assigned in this fixed order and never cycled -- a fifth condition would need
# the list extended and revalidated, not a generated hue.
SPEED_RAMP = ("#0072B2", "#C77F00", "#009E73", "#C0392B")

C_ZERO = "#262626"
FS_TICK = 18
FS_LABEL = 21
FS_LEGEND = 16
FS_ANNOT = 15
FIG_SIZE = (13.0, 7.8)


# ---------------------------------------------------------------------------
# panel definitions
#
# Each series is (value, spread, dash, marker, kind). value and spread name a
# column or take (table -> array); kind labels the series in the style legend.
# ---------------------------------------------------------------------------

def _est(bias_col):
    """Mean estimate over the valid mask: the truth plus that estimator's bias,
    which is what bias is defined as. Avoids a second pass over points.csv just
    to recover a number metrics_slice.csv already implies."""
    return lambda m: _col(m, "mean_true_mms") + _col(m, bias_col)


def _relative(bias_col):
    """Bias as a percentage of the condition's peak speed."""
    return lambda m: 100 * _col(m, bias_col) / _vmax(m)


PANELS = {
    "velocity_raw": dict(
        ylabel="Mean velocity over valid ROIs [cm/s]",
        headline="mean velocity, truth against raw estimate",
        zero=False, scale=0.1,
        series=[("mean_true_mms", None, (0, (5, 3)), None, "Truth (transit)"),
                (_est("bias_raw_mms"), "bias_raw_seed_std_mms", "-", "o",
                 "Raw estimate")],
    ),
    "velocity_sg": dict(
        ylabel="Mean velocity over valid ROIs [cm/s]",
        headline="mean velocity, truth against SG estimate",
        zero=False, scale=0.1,
        series=[("mean_true_mms", None, (0, (5, 3)), None, "Truth (transit)"),
                (_est("bias_sg_mms"), "bias_sg_seed_std_mms", "-", "o",
                 "After SG (2-D)")],
    ),
    "flow_fill": dict(
        ylabel="Flow rate [ml/min]",
        headline="flow rate, truth against no-slip fill",
        zero=False,
        series=[("flow_true_ml_min", None, (0, (5, 3)), None,
                 "Truth (row1 plane)"),
                ("flow_fill_ml_min", "flow_fill_ml_min_seed_std", "-", "o",
                 "No-slip fill")],
    ),
    "flow_sg": dict(
        ylabel="Flow rate [ml/min]",
        headline="flow rate, truth against SG estimate",
        zero=False,
        series=[("flow_true_ml_min", None, (0, (5, 3)), None,
                 "Truth (row1 plane)"),
                ("flow_sg_ml_min", "flow_sg_ml_min_seed_std", "-", "o",
                 "After SG (2-D)")],
    ),
    "rmse_raw": dict(
        ylabel="RMSE [cm/s]",
        headline="point RMSE, raw",
        zero=False, scale=0.1, floor=True,
        series=[("rmse_raw_mms", "rmse_raw_seed_std_mms", "-", "o", None)],
    ),
    "rmse_sg": dict(
        ylabel="RMSE [cm/s]",
        headline="point RMSE, after SG",
        zero=False, scale=0.1, floor=True,
        series=[("rmse_sg_mms", "rmse_sg_seed_std_mms", "-", "o", None)],
    ),
    "rmse_full_fill": dict(
        ylabel="Full-lumen RMSE [cm/s]",
        headline="full-lumen RMSE, no-slip fill",
        zero=False, scale=0.1, floor=True,
        series=[("rmse_fill_full_mms", "rmse_fill_full_mms_seed_std", "-", "o",
                 None)],
    ),
    "rmse_full_sg": dict(
        ylabel="Full-lumen RMSE [cm/s]",
        headline="full-lumen RMSE, after SG",
        zero=False, scale=0.1, floor=True,
        series=[("rmse_sg_full_mms", "rmse_sg_full_mms_seed_std", "-", "o",
                 None)],
    ),
    "nrmse_raw": dict(
        ylabel="NRMSE over valid ROIs [%]",
        headline="point NRMSE, raw",
        zero=False, scale=100.0, floor=True,
        series=[("nrmse_raw", "nrmse_raw_seed_std", "-", "o", None)],
    ),
    "nrmse_sg": dict(
        ylabel="NRMSE over valid ROIs [%]",
        headline="point NRMSE, after SG",
        zero=False, scale=100.0, floor=True,
        series=[("nrmse_sg", "nrmse_sg_seed_std", "-", "o", None)],
    ),
    "nrmse_full_fill": dict(
        ylabel="Full-lumen NRMSE [%]",
        headline="full-lumen NRMSE, no-slip fill",
        zero=False, scale=100.0, floor=True,
        series=[("nrmse_fill_full", "nrmse_fill_full_seed_std", "-", "o", None)],
    ),
    "nrmse_full_sg": dict(
        ylabel="Full-lumen NRMSE [%]",
        headline="full-lumen NRMSE, after SG",
        zero=False, scale=100.0, floor=True,
        series=[("nrmse_sg_full", "nrmse_sg_full_seed_std", "-", "o", None)],
    ),
    "flow_error_fill": dict(
        ylabel="Flow error [%]",
        headline="flow error, no-slip fill",
        zero=True,
        series=[("flow_error_fill_pct", "flow_error_fill_pct_seed_std", "-", "o",
                 None)],
    ),
    "flow_error_sg": dict(
        ylabel="Flow error [%]",
        headline="flow error, after SG",
        zero=True,
        series=[("flow_error_sg_pct", "flow_error_sg_pct_seed_std", "-", "o",
                 None)],
    ),
    "bias_raw": dict(
        ylabel="Bias [cm/s]",
        headline="bias, raw",
        zero=True, scale=0.1,
        series=[("bias_raw_mms", "bias_raw_seed_std_mms", "-", "o", None)],
    ),
    "bias_sg": dict(
        ylabel="Bias [cm/s]",
        headline="bias, after SG",
        zero=True, scale=0.1,
        series=[("bias_sg_mms", "bias_sg_seed_std_mms", "-", "o", None)],
    ),
    "bias_relative_raw": dict(
        ylabel="Bias / peak speed [%]",
        headline="bias relative to peak speed, raw",
        zero=True,
        series=[(_relative("bias_raw_mms"), None, "-", "o", None)],
    ),
    "bias_relative_sg": dict(
        ylabel="Bias / peak speed [%]",
        headline="bias relative to peak speed, after SG",
        zero=True,
        series=[(_relative("bias_sg_mms"), None, "-", "o", None)],
    ),
    "valid_frac": dict(
        ylabel="Valid ROIs [%]",
        headline="estimator yield",
        zero=False, scale=100.0, ylim=(0, 105),
        series=[("valid_frac", "valid_frac_seed_std", "-", "o", None)],
    ),
}


def metric_keys():
    """Panel keys in the order they are meant to be read."""
    return list(PANELS)


# ---------------------------------------------------------------------------


def collect_speeds(dataset_dirs):
    """One metrics_slice table for several dataset folders, tagged and ordered.

    Each folder is a results/<run>/<dataset> directory. The dataset column is
    filled from the folder name when the CSV predates it, and the speed column
    is the peak speed in cm/s that style.order_datasets reads off that name.
    """
    parts = []
    for d in dataset_dirs:
        csv = os.path.join(d, "metrics_slice.csv")
        if not os.path.isfile(csv):
            print(f"skip: no metrics_slice.csv in {d}")
            continue
        T = pd.read_csv(csv)
        if "dataset" not in T.columns:
            T["dataset"] = os.path.basename(os.path.abspath(d))
        parts.append(T)
    if not parts:
        raise ValueError("none of the folders held a metrics_slice.csv")
    P = pd.concat(parts, ignore_index=True)
    names, speed = st.order_datasets(P["dataset"])
    P["speed_cms"] = P["dataset"].map(dict(zip(names, speed)))
    return P


def plot_speed_metric(P, key, title=""):
    """One metric, one figure, every speed condition drawn in it."""
    if key not in PANELS:
        raise ValueError(f"unknown metric {key!r}; have {metric_keys()}")
    spec = PANELS[key]
    names, speed = st.order_datasets(P["dataset"])
    labels = _speed_labels(names, speed)
    colours = speed_colors(len(names))
    scale = spec.get("scale", 1.0)

    fig = st.new_figure(FIG_SIZE)
    ax = fig.subplots()
    if spec["zero"]:
        ax.axhline(0, color=C_ZERO, lw=1.4, zorder=2)

    ends = []            # (x, y, label, [(is_truth, value), ...]) per condition
    drew = False
    for j, (nm, lab, colour) in enumerate(zip(names, labels, colours)):
        m = P[P["dataset"] == nm].sort_values("slice_pos_mm")
        if m.empty:
            continue
        x = m["slice_pos_mm"].to_numpy(float)
        n_seeds = _col(m, "n_seeds", 1.0)
        last = []
        for k, (value, spread, dash, marker, _) in enumerate(spec["series"]):
            v = _resolve(m, value) * scale
            if not np.isfinite(v).any():
                continue
            drew = True
            ok = np.flatnonzero(np.isfinite(v))
            last.append((marker is None, v[ok[-1]]))   # dashed truth has no marker
            # The spread is a band rather than a bar per slice: 41 slices times
            # up to eight series makes bars the loudest thing on the canvas,
            # which is backwards -- the seed spread is context, not the subject.
            sd = _resolve(m, spread) * scale if spread is not None else None
            if sd is not None and (n_seeds > 1).any():
                st.band(ax, x, v, sd, colour, alpha=0.13, zorder=1)
            # Markers identify a series without claiming every slice is a
            # measured event worth its own dot; offset so they interleave.
            ax.plot(x, v, linestyle=dash, marker=marker, color=colour, lw=2.0,
                    ms=7, mew=0, markevery=(2 * j + 4 * k, 8), zorder=3)
        v = _resolve(m, spec["series"][-1][0]) * scale
        ok = np.flatnonzero(np.isfinite(v))
        if ok.size:
            ends.append((x[ok[-1]], v[ok[-1]], lab, last))

    st.style_axes(ax)
    ax.tick_params(labelsize=FS_TICK)
    ax.set_xlabel("Slice position [mm]", fontsize=FS_LABEL, labelpad=10)
    ax.set_ylabel(spec["ylabel"], fontsize=FS_LABEL, labelpad=10)
    if "ylim" in spec:
        ax.set_ylim(*spec["ylim"])
    elif spec.get("floor"):
        # A magnitude reads against zero; letting the axis start mid-range
        # exaggerates whatever variation is left.
        ax.set_ylim(bottom=0)
    # The legends stack in the strip between the axes and the title -- speed
    # above, line style below -- so neither can run into the other however wide
    # its entries get, and the title is padded past both.
    st.set_title(ax, f"{title}  {spec['headline']}".strip(), pad=72)

    _speed_legend(ax, labels, colours)
    _style_legend(ax, spec["series"])
    _end_labels(ax, ends)
    _seed_note(ax, P)
    if not drew:
        ax.text(0.5, 0.5, "not in these results files", transform=ax.transAxes,
                ha="center", va="center", fontsize=FS_ANNOT, style="italic",
                color="0.45")
    return fig


def speed_colors(n):
    """The first n hues, in order. Beyond the list the caller is out of
    validated colours: raise rather than cycle or invent one, because a repeated
    hue silently relabels a condition."""
    if n > len(SPEED_RAMP):
        raise ValueError(
            f"{n} conditions but only {len(SPEED_RAMP)} validated hues; extend "
            "SPEED_RAMP and re-run the palette check, or facet the figure")
    return list(SPEED_RAMP[:n])


# ---------------------------------------------------------------------------


def _col(T, name, default=np.nan):
    if name in T.columns:
        return T[name].to_numpy(float)
    return np.full(len(T), default, float)


def _resolve(T, value):
    """A series definition is either a column name or a function of the table."""
    if callable(value):
        return np.asarray(value(T), float)
    return _col(T, value)


def _speed_labels(names, speed):
    """Condition labels in cm/s, matching the velocity axes in this set.

    The pipeline works in mm/s throughout and every CSV column is named for it,
    but these figures are read as cm/s, so the conversion happens once here and
    in the panel scale rather than in the reader's head.
    """
    return [f"{s:g} cm/s" if np.isfinite(s) else str(n)
            for n, s in zip(names, speed)]


def _vmax(T):
    """Peak speed in mm/s for each row, from the dataset name."""
    names, speed = st.order_datasets(T["dataset"])
    lut = {n: s * 10.0 for n, s in zip(names, speed)}   # cm/s -> mm/s
    return T["dataset"].map(lut).to_numpy(float)


def _speed_legend(ax, labels, colours):
    """Colour is speed. Kept as a real legend even though the curves are also
    labelled at their ends: the ends can only take one label per condition, and
    a panel with two estimators draws two curves per colour."""
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=c, lw=3.0, label=l)
               for c, l in zip(colours, labels)]
    leg = ax.legend(handles=handles, bbox_to_anchor=(0.0, 1.062), loc="lower left",
                    fontsize=FS_LEGEND, frameon=False, ncols=len(handles),
                    title="Peak speed (vmax)", title_fontsize=FS_LEGEND,
                    alignment="left", borderaxespad=0.0)
    ax.add_artist(leg)
    return leg


def _style_legend(ax, series):
    """Line style is which estimator or which of truth/estimate, drawn in
    neutral ink so it cannot be mistaken for a fifth speed. Dropped when the
    panel has only one series."""
    from matplotlib.lines import Line2D
    named = [(dash, marker, kind) for _, _, dash, marker, kind in series if kind]
    if len(named) < 2:
        return None
    handles = [Line2D([], [], color="#4A4A4A", lw=2.0, linestyle=dash,
                      marker=marker, ms=7, mew=0, label=kind)
               for dash, marker, kind in named]
    return ax.legend(handles=handles, bbox_to_anchor=(0.0, 1.004),
                     loc="lower left", fontsize=FS_LEGEND, frameon=False,
                     ncols=len(handles), borderaxespad=0.0)


def _end_labels(ax, ends):
    """Condition and value printed at the end of each curve.

    The value is spelled out because the alternative is reading it off the axis:
    the label names the condition, so without a number beside it the reader
    still has to trace the line back to a gridline. It is the value at that
    x position, not a mean over the vessel, since that is where the text sits.

    Dropped entirely as soon as two curves finish closer together than a label
    is tall: overlapping text is unreadable, and nudging it apart would put each
    number beside a neighbour's line and read as data. The legend still carries
    identity in that case.
    """
    if not ends:
        return
    ylo, yhi = ax.get_ylim()
    step = 0.10 * (yhi - ylo)          # two text lines per label
    order = sorted(ends, key=lambda e: e[1])
    gaps = [b[1] - a[1] for a, b in zip(order, order[1:])]
    if gaps and min(gaps) < step:
        return
    # Freeze the ticks before widening, or the locator adds one out in the
    # margin and the axis claims slices that were never measured.
    lo, hi = ax.get_xlim()
    ax.set_xticks([t for t in ax.get_xticks() if lo - 1e-9 <= t <= hi + 1e-9])
    ax.set_xlim(lo, hi + 0.24 * (hi - lo))
    for x, y, label, last in order:
        ax.annotate(f"{label}\n{_value_text(last)}", (x, y),
                    textcoords="offset points", xytext=(12, 0),
                    ha="left", va="center", fontsize=FS_ANNOT,
                    color="#2E2E2E", annotation_clip=False, linespacing=1.35)


def _value_text(last):
    """The numbers for one condition: the estimate, and the truth it is measured
    against when the panel carries one."""
    truth = [v for is_truth, v in last if is_truth]
    est = [v for is_truth, v in last if not is_truth]
    if est and truth:
        return f"{_fmt(est[-1])}  (true {_fmt(truth[-1])})"
    if est:
        return _fmt(est[-1])
    return _fmt(truth[-1]) if truth else ""


def _fmt(v):
    """Enough digits to be useful, no more than the value supports."""
    a = abs(v)
    if a >= 100:
        return f"{v:.0f}"
    if a >= 10:
        return f"{v:.1f}"
    return f"{v:.2f}"


def _seed_note(ax, P):
    n = _col(P, "n_seeds")
    n = n[np.isfinite(n)]
    if n.size == 0 or n.max() <= 1:
        return
    common = int(np.bincount(n.astype(int)).argmax())
    ax.text(1.0, -0.145, f"shaded band: SD over {common} seeds",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=FS_ANNOT, color="0.35")
