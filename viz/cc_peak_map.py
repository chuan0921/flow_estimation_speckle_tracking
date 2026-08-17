"""Figure made from points.csv: the through-plane CC peak map, slice by slice.

Ports viz.plot_cc_peak_map. One tile per slice in the x-z plane; marker colour
is the selected peak correlation and grey crosses mark ROIs where no peak
passed the gate, so a tile shows both how good the peaks were and where the
estimator gave up.
"""

from __future__ import annotations

import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from . import style as st

# Beyond this the tiles are too small to read; the excess slices are sampled
# out evenly and the drop is reported, never silent.
DEFAULT_MAX_SLICES = 16


def plot_cc_peak_map(P, title="", cc_col="ccA", hide_outliers=False, marker_size=44,
                     max_slices=DEFAULT_MAX_SLICES):
    for need in ("x_mm", "z_mm", cc_col):
        if need not in P.columns:
            raise ValueError(f"points table must contain a {need} column")
    # Every seed samples the same lattice, so repeats would overprint each
    # other at identical coordinates and inflate the per-tile n.
    P = st.representative_seed(P)

    if "slice_idx" in P.columns:
        slices = list(dict.fromkeys(P["slice_idx"]))
    else:
        P = P.assign(slice_idx=1)
        slices = [1]
    if np.isfinite(max_slices) and len(slices) > max_slices:
        keep = np.unique(np.round(np.linspace(0, len(slices) - 1, int(max_slices))).astype(int))
        print(f"cc_peak_map: showing {keep.size} of {len(slices)} slices")
        slices = [slices[i] for i in keep]

    n = len(slices)
    ncol = int(np.ceil(np.sqrt(n)))
    nrow = int(np.ceil(n / ncol))
    fig = st.new_figure((4.2 * ncol + 2.6, 4.1 * nrow + 1.4))
    axes = fig.subplots(nrow, ncol, squeeze=False)

    norm = Normalize(*st.span(P[cc_col]))
    xlim = st.span(P["x_mm"])
    zlim = st.span(P["z_mm"])
    fs_tile = max(11, st.FS_ANNOT - int(1.5 * max(0, ncol - 2)))

    for k, ax in enumerate(axes.ravel()):
        if k >= n:
            ax.set_visible(False)
            continue
        Q = P[P["slice_idx"] == slices[k]]
        if hide_outliers and "outlier" in Q.columns:
            Q = Q[~Q["outlier"].astype(bool)]
        ok = np.isfinite(Q[cc_col])

        ax.scatter(Q["x_mm"][~ok], Q["z_mm"][~ok], s=26, c=st.C_MISS,
                   marker="x", linewidths=1.0, zorder=2)
        ax.scatter(Q["x_mm"][ok], Q["z_mm"][ok], s=marker_size, c=Q[cc_col][ok],
                   cmap=st.CMAP_SEQ, norm=norm, edgecolors="none", zorder=3)
        _draw_vessel(ax, Q)

        ax.set_aspect("equal")
        ax.set_xlim(*xlim)
        ax.set_ylim(zlim[1], zlim[0])          # z grows downwards
        st.style_axes(ax)
        ax.tick_params(labelsize=fs_tile)
        ax.set_title(_tile_title(Q, slices[k], cc_col),
                     fontsize=fs_tile, pad=8, linespacing=1.3)
        if k + ncol >= n:                      # nothing below this tile
            ax.set_xlabel("x [mm]", fontsize=fs_tile + 2, labelpad=6)
        if k % ncol == 0:
            ax.set_ylabel("z [mm]", fontsize=fs_tile + 2, labelpad=6)

    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=st.CMAP_SEQ),
                      ax=axes.ravel().tolist(), fraction=0.035, pad=0.02)
    cb.set_label(f"{cc_col} peak CC", fontsize=st.FS_LABEL, labelpad=12)
    cb.ax.tick_params(labelsize=st.FS_TICK)
    st.set_suptitle(fig, f"{title}  ROI through-plane peak CC by slice".strip())
    return fig


def _draw_vessel(ax, Q):
    """Outline the lumen at the outermost sampled radius."""
    if "r_mm" not in Q.columns or not np.isfinite(Q["r_mm"]).any():
        return
    R = float(np.nanmax(Q["r_mm"]))
    th = np.linspace(0, 2 * np.pi, 240)
    ax.plot(R * np.cos(th), R * np.sin(th), "-", color="#333333", lw=1.2,
            zorder=1)


def _tile_title(Q, slice_value, cc_col):
    """Two lines -- who the tile is, then how it did. One line is too wide to
    survive next to its neighbour once the grid has more than one column."""
    cc = Q[cc_col]
    head = f"Slice {slice_value:g}"
    finite = cc[np.isfinite(cc)]
    n_ok = finite.size
    valid = 100 * n_ok / max(len(cc), 1)
    if n_ok:
        q25, q75 = np.percentile(finite, [25, 75])
        stats = (f"valid={n_ok}/{len(cc)} ({valid:.1f}%)\n"
                 f"CC med={np.median(finite):.2f} [IQR {q25:.2f}-{q75:.2f}]")
    else:
        stats = f"valid=0/{len(cc)} (0.0%), CC unavailable"
    return f"{head}\n{stats}"
