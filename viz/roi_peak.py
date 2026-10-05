"""ROI peak-lag diagnostics for one slice, from points.csv.

The transit estimate is a lag in frames, and every downstream number is that
lag in disguise. These figures show the lag itself, ROI by ROI: the full
cross-section, plus 3x3 zooms at the axis, mid-radius and near the wall.

Everything is derived from what main() already writes:

    lag_true_frames  = lag_peak_frames * Vy_mms / VyTrue_mms
    lag_error_frames = lag_peak_frames - lag_true_frames
    roi_id           = row order within the slice

All four figures of a slice share one symmetric colour scale, so a zoom is
literally a window onto the full map.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

from . import style as st

CMAP = "turbo"
PAD_TILE = 0.60          # frame padding around a zoom, in tile pitches
TILE_IN = 0.92           # full-map tile size [inch]; text has to fit inside
C_NAN = "#D1D1D1"

# radius fractions the three zooms are centred on
ZOOMS = (("center", 0.0), ("middle", 0.5), ("nearwall", 0.82))


def roi_peak_frame(P, slice_idx=None):
    """One row per ROI for a single slice, with the lag columns derived.

    Takes the points table (any number of slices/seeds), keeps one slice and
    one realisation, and returns the columns the figures need.
    """
    d = st.representative_seed(P)
    if slice_idx is not None:
        d = d[d["slice_idx"] == slice_idx]
    if "slice_idx" in d.columns and d["slice_idx"].nunique() > 1:
        raise ValueError("pass a single slice: got "
                         f"{d['slice_idx'].nunique()} in the table")
    if d.empty:
        raise ValueError("no rows for that slice")

    d = d.reset_index(drop=True)
    lag_peak = d["lag_peak_frames"].to_numpy(float)
    ratio = d["Vy_mms"].to_numpy(float) / d["VyTrue_mms"].to_numpy(float)
    lag_true = lag_peak * ratio
    return pd.DataFrame({
        "roi_id": np.arange(1, len(d) + 1),
        "x_mm": d["x_mm"].to_numpy(float),
        "z_mm": d["z_mm"].to_numpy(float),
        "r_mm": d["r_mm"].to_numpy(float),
        "cc_peak": d["ccA"].to_numpy(float),
        "vy_mms": d["Vy_mms"].to_numpy(float),
        "lag_peak_frames": lag_peak,
        "lag_true_frames": lag_true,
        "lag_error_frames": lag_peak - lag_true,
    })


def lag_error_limit(d, pct=90.0):
    """Symmetric colour limit shared by a slice's four figures.

    A percentile, not the maximum: lag_true = d_row/(VyTrue*dt) is singular
    where the true velocity approaches zero, so a handful of near-wall ROIs in
    a recirculating slice reach hundreds of frames and would flatten every
    other tile to the same green. Tiles past the limit saturate, and the
    colorbar says so; their printed est/true values are untouched.
    """
    err = d["lag_error_frames"].to_numpy(float)
    err = np.abs(err[np.isfinite(err)])
    if err.size == 0:
        return 1.0
    return max(float(np.percentile(err, pct)), 1.0)


def velocity_limit(d, pct=98.0):
    """Colour ceiling for the truth-free velocity mode."""
    vy = d["vy_mms"].to_numpy(float)
    vy = np.abs(vy[np.isfinite(vy)])
    if vy.size == 0:
        return 1.0
    return max(float(np.percentile(vy, pct)), 1.0)


def plot_roi_peak_map(d, title="", lim=None, labels=True, color="error",
                      z0=None):
    """Every ROI of the slice as a tile, coloured by lag error.

    Keeps the per-tile text, which is what makes this a lookup table rather
    than just a picture -- at the cost of a large canvas.

    color="velocity" is the truth-free mode for experimental data: tiles are
    coloured by the estimated velocity and the true-lag line is dropped.
    """
    dx, dz = _pitch(d)
    if lim is None:
        lim = velocity_limit(d) if color == "velocity" else lag_error_limit(d)
    ux, uz = np.unique(np.round(d["x_mm"], 6)), np.unique(np.round(d["z_mm"], 6))
    ax_w = TILE_IN * (len(ux) + 2 * PAD_TILE)
    ax_h = ax_w * (len(uz) + 2 * PAD_TILE) / (len(ux) + 2 * PAD_TILE) * (dz / dx)
    fs_roi = max(5.0, TILE_IN * 72 / 7.0)
    return _draw(d, dx, dz, lim, ax_w, ax_h, fs_roi, labels,
                 f"{title}  ROI peak lag map".strip(), color, z0)


def plot_roi_peak_3x3(d, target, name, title="", lim=None, color="error",
                      z0=None):
    """The 3x3 block centred on the ROI closest to r/R = target."""
    dx, dz = _pitch(d)
    if lim is None:
        lim = velocity_limit(d) if color == "velocity" else lag_error_limit(d)
    keep = _pick_3x3(d, dx, dz, target)
    if not keep.any():
        raise ValueError(f"no 3x3 block found at r/R = {target}")
    return _draw(d[keep], dx, dz, lim, 7.6, None, 25.0, True,
                 f"{title}  {name} 3x3 ROI peak check".strip(), color, z0)


# ---------------------------------------------------------------------------


def _pitch(d):
    """Lattice spacing in x and z [mm]."""
    ux = np.unique(np.round(d["x_mm"], 6))
    uz = np.unique(np.round(d["z_mm"], 6))
    dx = np.median(np.diff(ux)) if ux.size > 1 else 0.25
    dz = np.median(np.diff(uz)) if uz.size > 1 else 0.25
    return float(dx), float(dz)


def _pick_3x3(d, dx, dz, target):
    """Boolean mask of the 3x3 block around the ROI closest to r/R = target."""
    r = d["r_mm"].to_numpy(float)
    R = np.nanmax(r)
    c0 = int(np.nanargmin(r)) if target == 0 else int(np.nanargmin(np.abs(r / R - target)))
    x0, z0 = d["x_mm"].iloc[c0], d["z_mm"].iloc[c0]
    return ((np.abs(d["x_mm"] - x0) < 1.5 * dx)
            & (np.abs(d["z_mm"] - z0) < 1.5 * dz)).to_numpy()


def _draw(d, dx, dz, lim, ax_w, ax_h, fs_roi, labels, title, color="error",
          z0=None):
    # z0: absolute depth of the vessel centre [mm]. When given, the z axis
    # is relabelled as absolute image depth (tick positions unchanged).
    velocity = color == "velocity"
    x = d["x_mm"].to_numpy(float)
    z = d["z_mm"].to_numpy(float)
    xlim = (x.min() - PAD_TILE * dx, x.max() + PAD_TILE * dx)
    zlim = (z.min() - PAD_TILE * dz, z.max() + PAD_TILE * dz)
    if ax_h is None:
        ax_h = ax_w * (zlim[1] - zlim[0]) / (xlim[1] - xlim[0])

    mar_l, mar_b, mar_t = 1.8, 1.25, 0.9
    gap_cb, cb_w, mar_r = 0.45, 0.26, 1.45
    fig_w = mar_l + ax_w + gap_cb + cb_w + mar_r
    fig_h = mar_b + ax_h + mar_t

    fig = st.new_figure((fig_w, fig_h))
    fig.set_layout_engine("none")           # the layout is computed here
    ax = fig.add_axes([mar_l / fig_w, mar_b / fig_h, ax_w / fig_w, ax_h / fig_h])

    cmap = st.plt.get_cmap(CMAP)
    norm = Normalize(vmin=0, vmax=lim) if velocity else Normalize(vmin=-lim, vmax=lim)
    for i in range(len(d)):
        val = (d["vy_mms"] if velocity else d["lag_error_frames"]).iloc[i]
        face = cmap(norm(val)) if np.isfinite(val) else (0.82, 0.82, 0.82, 1.0)
        ax.add_patch(Rectangle((x[i] - 0.46 * dx, z[i] - 0.46 * dz),
                               0.92 * dx, 0.92 * dz, facecolor=face,
                               edgecolor="k", linewidth=1.2, zorder=2))
        if not labels:
            continue
        est, tru, cc = (d["lag_peak_frames"].iloc[i], d["lag_true_frames"].iloc[i],
                        d["cc_peak"].iloc[i])
        if velocity:
            vy = d["vy_mms"].iloc[i]
            txt = (f"ROI {d['roi_id'].iloc[i]:.0f}\n"
                   f"{vy:.0f} mm/s\nlag {est:.0f}\nCC {cc:.2f}"
                   if np.isfinite(est) else
                   f"ROI {d['roi_id'].iloc[i]:.0f}\nno est\n—\nCC NaN")
        else:
            txt = (f"ROI {d['roi_id'].iloc[i]:.0f}\n"
                   f"est {est:.0f}\ntrue {tru:.0f}\nCC {cc:.2f}"
                   if np.isfinite(est) else
                   f"ROI {d['roi_id'].iloc[i]:.0f}\nest NaN\ntrue {tru:.0f}\nCC NaN")
        lum = 0.299 * face[0] + 0.587 * face[1] + 0.114 * face[2]
        ax.text(x[i], z[i], txt, ha="center", va="center", zorder=3,
                fontsize=fs_roi, fontweight="bold",
                color="w" if lum < 0.45 else "k", linespacing=1.3)

    ax.set_xlim(*xlim)
    ax.set_ylim(zlim[1], zlim[0])           # z grows downwards
    ax.set_aspect("equal")
    ax.set_xticks(np.unique(np.round(x, 6))[::max(1, len(np.unique(x)) // 12)])
    ax.set_yticks(np.unique(np.round(z, 6))[::max(1, len(np.unique(z)) // 12)])
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.2f}")
    if z0 is None:
        ax.yaxis.set_major_formatter(lambda v, _: f"{v:.2f}")
        zlabel = "Row1 ROI centre z [mm]"
    else:
        ax.yaxis.set_major_formatter(lambda v, _: f"{v + z0:.2f}")
        zlabel = "Image depth z [mm]"
    st.style_axes(ax, grid=False)
    st.label_axes(ax, "Row1 ROI centre x [mm]", zlabel)
    st.set_title(ax, title)

    cb_h = 0.70 * ax_h
    cax = fig.add_axes([(mar_l + ax_w + gap_cb) / fig_w,
                        (mar_b + 0.5 * (ax_h - cb_h)) / fig_h,
                        cb_w / fig_w, cb_h / fig_h])
    cb = fig.colorbar(st.plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    if velocity:
        vy = d["vy_mms"].to_numpy(float)
        worst = np.nanmax(np.abs(vy)) if np.isfinite(vy).any() else 0.0
        clipped = f", clipped at {lim:.0f} of {worst:.0f}" if worst > 1.01 * lim else ""
        cb.set_label(f"estimated $v_y$ [mm/s]{clipped}",
                     fontsize=st.FS_LABEL, labelpad=14)
    else:
        err = d["lag_error_frames"].to_numpy(float)
        worst = np.nanmax(np.abs(err)) if np.isfinite(err).any() else 0.0
        clipped = f", clipped at $\\pm${lim:.0f} of {worst:.0f}" if worst > 1.01 * lim else ""
        cb.set_label(f"lag error: est - true [frames]{clipped}",
                     fontsize=st.FS_LABEL, labelpad=14)
    cb.ax.tick_params(labelsize=st.FS_TICK)
    return fig
