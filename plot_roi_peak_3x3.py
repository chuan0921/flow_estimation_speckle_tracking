#!/usr/bin/env python3
"""Draw representative 3x3 ROI peak-lag diagnostic maps.

Reads the ROI peak table (CSV) written by the MATLAB pipeline plus the slice
.mat file (for the absolute ROI centre xc/zc) and writes three PNGs:
center, middle, nearwall.

The layout is laid out explicitly in inches so the tiles, the fonts and the
colorbar keep the same proportions no matter what the data range is:

    | left | axes (square, data aspect) | gap | cbar | cbar labels |

Usage
-----
python3 plot_roi_peak_3x3.py \
    --csv     results/sim_stenosis50_80cms/figures/slice_028_roi_peak_tile_map_autolag.csv \
    --slice   /home/chihyu/psf_sim/data/sim_stenosis50_80cms/slices/slice_028.mat \
    --prefix  slice_028_roi
"""

from __future__ import annotations

import argparse
import os
import re

import h5py
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

REQUIRED = [
    "roi_id", "x_mm", "z_mm", "r_mm", "cc_peak",
    "lag_peak_frames", "lag_true_frames", "lag_error_frames",
]

# ---- type sizes (pt). Everything readable at a glance, colorbar kept slim ---
FS_TICK = 22
FS_LABEL = 26
FS_TITLE = 30
FS_ROI = 25
FS_CB_TICK = 20
FS_CB_LABEL = 23

# ---- layout (inches) -------------------------------------------------------
AX_W = 7.6          # plot box width; height follows the data aspect ratio
MAR_L = 1.80        # y tick labels + y axis label
MAR_B = 1.25        # x tick labels + x axis label
MAR_T = 0.85        # title
GAP_CB = 0.45       # plot box -> colorbar
CB_W = 0.26         # colorbar width (slim on purpose)
MAR_R = 1.45        # colorbar tick labels + colorbar label
PAD_TILE = 0.60     # frame padding around the 3x3 block, in tile pitches


def read_slice_centre(slice_file):
    """Return (xc, zc) in metres from a v7.3 (HDF5) or older .mat file."""
    try:
        with h5py.File(slice_file, "r") as f:
            return float(np.array(f["xc"]).ravel()[0]), float(np.array(f["zc"]).ravel()[0])
    except OSError:
        from scipy.io import loadmat
        m = loadmat(slice_file, variable_names=["xc", "zc"])
        return float(np.ravel(m["xc"])[0]), float(np.ravel(m["zc"])[0])


def pick_3x3(d, ux, uz, target):
    """Indices of the 3x3 block centred on the ROI closest to r/R = target."""
    if target == 0:
        c0 = int(np.nanargmin(d["r_mm"].to_numpy()))
    else:
        c0 = int(np.nanargmin(np.abs(d["rr"].to_numpy() - target)))
    ix0 = int(np.argmin(np.abs(ux - d["xr"].iloc[c0])))
    iz0 = int(np.argmin(np.abs(uz - d["zr"].iloc[c0])))
    xs = ux[max(ix0 - 1, 0):ix0 + 2]
    zs = uz[max(iz0 - 1, 0):iz0 + 2]
    keep = d["xr"].isin(xs) & d["zr"].isin(zs) & np.isfinite(d["r_mm"])
    return d.index[keep]


def draw_one(d, idx, dx, dz, lim, out_png, title):
    x = d.loc[idx, "x"].to_numpy()
    z = d.loc[idx, "z"].to_numpy()
    est = d.loc[idx, "lag_peak_frames"].to_numpy()
    tru = d.loc[idx, "lag_true_frames"].to_numpy()
    err = d.loc[idx, "lag_error_frames"].to_numpy()
    cc = d.loc[idx, "cc_peak"].to_numpy()
    rid = d.loc[idx, "roi_id"].to_numpy()

    xlim = (x.min() - PAD_TILE * dx, x.max() + PAD_TILE * dx)
    zlim = (z.min() - PAD_TILE * dz, z.max() + PAD_TILE * dz)

    ax_h = AX_W * (zlim[1] - zlim[0]) / (xlim[1] - xlim[0])
    fig_w = MAR_L + AX_W + GAP_CB + CB_W + MAR_R
    fig_h = MAR_B + ax_h + MAR_T

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([MAR_L / fig_w, MAR_B / fig_h, AX_W / fig_w, ax_h / fig_h])

    cmap = plt.get_cmap("turbo")
    norm = Normalize(vmin=-lim, vmax=lim)

    for k in range(len(idx)):
        face = cmap(norm(err[k])) if np.isfinite(err[k]) else (0.82, 0.82, 0.82, 1.0)
        ax.add_patch(Rectangle(
            (x[k] - 0.46 * dx, z[k] - 0.46 * dz), 0.92 * dx, 0.92 * dz,
            facecolor=face, edgecolor="k", linewidth=2.0, zorder=2))

        if np.isfinite(est[k]):
            label = f"ROI {rid[k]:.0f}\nest {est[k]:.0f}\ntrue {tru[k]:.0f}\nCC {cc[k]:.2f}"
        else:
            label = f"ROI {rid[k]:.0f}\nest NaN\ntrue {tru[k]:.0f}\nCC NaN"
        lum = 0.299 * face[0] + 0.587 * face[1] + 0.114 * face[2]
        ax.text(x[k], z[k], label, ha="center", va="center", zorder=3,
                fontsize=FS_ROI, fontweight="bold",
                color="w" if lum < 0.45 else "k", linespacing=1.35)

    ax.set_xlim(*xlim)
    ax.set_ylim(zlim[1], zlim[0])          # z grows downwards
    ax.set_aspect("equal")
    ax.set_xticks(np.unique(np.round(x, 6)))
    ax.set_yticks(np.unique(np.round(z, 6)))
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.2f}")
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:.2f}")
    ax.tick_params(labelsize=FS_TICK, direction="out", length=7, width=1.4, pad=8)
    for s in ax.spines.values():
        s.set_linewidth(1.4)
    ax.set_xlabel("Row1 ROI centre x [mm]", fontsize=FS_LABEL, labelpad=12)
    ax.set_ylabel("Row1 ROI centre z [mm]", fontsize=FS_LABEL, labelpad=12)
    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", pad=16)

    # slim colorbar, 70% of the plot height, docked next to the plot box
    cb_h = 0.70 * ax_h
    cax = fig.add_axes([
        (MAR_L + AX_W + GAP_CB) / fig_w,
        (MAR_B + 0.5 * (ax_h - cb_h)) / fig_h,
        CB_W / fig_w,
        cb_h / fig_h,
    ])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label("lag error: est - true [frames]", fontsize=FS_CB_LABEL, labelpad=14)
    cb.ax.tick_params(labelsize=FS_CB_TICK, length=6, width=1.2)
    cb.outline.set_linewidth(1.2)

    fig.savefig(out_png, dpi=150, facecolor="w")
    plt.close(fig)
    return out_png


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--slice", required=True)
    ap.add_argument("--outdir", default="")
    ap.add_argument("--prefix", default="")
    ap.add_argument("--title-prefix", default="")
    ap.add_argument("--radius-targets", type=float, nargs=3, default=(0.0, 0.5, 0.82))
    args = ap.parse_args()

    out_dir = args.outdir or os.path.dirname(os.path.abspath(args.csv))
    os.makedirs(out_dir, exist_ok=True)
    prefix = args.prefix or os.path.splitext(os.path.basename(args.csv))[0]

    ttl_prefix = args.title_prefix
    if not ttl_prefix:
        m = re.search(r"[Ss]lice[_\-\s]*0*(\d+)", prefix)
        ttl_prefix = f"Slice {m.group(1)}" if m else ""

    d = pd.read_csv(args.csv)
    missing = [c for c in REQUIRED if c not in d.columns]
    if missing:
        raise SystemExit(f"CSV is missing required column(s): {', '.join(missing)}")

    xc, zc = read_slice_centre(args.slice)
    d["x"] = d["x_mm"] + xc * 1e3
    d["z"] = d["z_mm"] + zc * 1e3
    d["xr"] = d["x"].round(6)
    d["zr"] = d["z"].round(6)
    d["rr"] = d["r_mm"] / np.nanmax(d["r_mm"])

    ux = np.unique(d["xr"])
    uz = np.unique(d["zr"])
    dx = np.median(np.diff(ux)) if ux.size > 1 else 0.25
    dz = np.median(np.diff(uz)) if uz.size > 1 else 0.25

    lim = np.nanmax(np.abs(d["lag_error_frames"]))
    if not np.isfinite(lim) or lim == 0:
        lim = 1.0

    out = []
    for name, target in zip(("center", "middle", "nearwall"), args.radius_targets):
        idx = pick_3x3(d, ux, uz, target)
        title = f"{ttl_prefix} {name} 3x3 ROI peak check".strip()
        out.append(draw_one(d, idx, dx, dz, lim,
                            os.path.join(out_dir, f"{prefix}_3x3_{name}.png"), title))
    for p in out:
        print(p)


if __name__ == "__main__":
    main()
