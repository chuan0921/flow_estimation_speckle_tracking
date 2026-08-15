"""Figure made from points.csv: the lumen cross-section vector field.

Ports viz.plot_vector_field. Marker colour is the through-plane component Vy,
arrows are the in-plane (Vx,Vz) pair. Each panel scales its arrows to its own
peak in-plane speed -- printed under the axes -- because the error field is
typically an order of magnitude below the flow itself.
"""

from __future__ import annotations

import numpy as np

from . import style as st


def plot_vector_field(T, title="", truth=False, hide_outliers=True,
                      marker_size=90, arrow_scale=1.0):
    # One realisation only: repeats sit on identical coordinates, so their
    # arrows and markers would just be drawn on top of each other.
    T = st.representative_seed(T)
    m = (np.isfinite(T["Vy_mms"]) & np.isfinite(T["Vx_mms"])
         & np.isfinite(T["Vz_mms"]))
    if hide_outliers and "outlier" in T.columns:
        m &= ~T["outlier"].astype(bool)
    if truth:
        m &= np.isfinite(T["VyTrue_mms"])
    if not m.any():
        raise ValueError("nothing to plot: no finite rows in the points table")

    Q = T[m]
    x = Q["x_mm"].to_numpy(float)
    z = Q["z_mm"].to_numpy(float)
    vx, vz, vy = (Q["Vx_mms"].to_numpy(float), Q["Vz_mms"].to_numpy(float),
                  Q["Vy_mms"].to_numpy(float))

    if truth:
        tx, tz, ty = (Q["VxTrue_mms"].to_numpy(float),
                      Q["VzTrue_mms"].to_numpy(float),
                      Q["VyTrue_mms"].to_numpy(float))
        panels = [(vx, vz, vy, "estimate", False),
                  (tx, tz, ty, "CFD truth", False),
                  (vx - tx, vz - tz, vy - ty, "error (est - CFD)", True)]
        shared = st.span(np.concatenate([vy, ty]))
    else:
        panels = [(vx, vz, vy, "estimate", False)]
        shared = st.span(vy)

    gs = _grid_spacing(x, z)
    fig = st.new_figure((6.6 * len(panels) + 1.0, 7.6))
    axes = fig.subplots(1, len(panels), squeeze=False)[0]

    for ax, (u, v, c, name, diverge) in zip(axes, panels):
        if diverge:
            L = max(float(np.nanmax(np.abs(c))), 1e-12)
            clim, cmap = (-L, L), st.CMAP_DIV
        else:
            clim, cmap = shared, st.CMAP_SEQ

        speed = np.hypot(u, v)
        vmax = float(np.nanmax(speed))
        length = arrow_scale * gs / max(vmax, 1e-12)   # mm of arrow per mm/s

        sc = ax.scatter(x, z, s=marker_size, c=c, cmap=cmap, vmin=clim[0],
                        vmax=clim[1], edgecolors="none", zorder=2)
        ax.quiver(x, z, length * u, length * v, angles="xy",
                  scale_units="xy", scale=1, color=st.C_LINE, width=0.004,
                  headwidth=4, headlength=5, zorder=3)

        ax.set_aspect("equal")
        ax.invert_yaxis()                              # z grows downwards
        st.style_axes(ax, grid=False)
        st.label_axes(ax, f"x [mm]   (arrow = {vmax:.2g} mm/s)", None)
        ax.set_title(name, fontsize=st.FS_LABEL, fontweight="bold", pad=12)
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("$v_y$ [mm/s]", fontsize=st.FS_ANNOT, labelpad=10)
        cb.ax.tick_params(labelsize=st.FS_ANNOT)

    st.label_axes(axes[0], None, "z [mm]  (depth)")
    hidden = ", outliers hidden" if hide_outliers else ""
    st.set_suptitle(fig, f"{title}  {int(m.sum())} points{hidden}".strip())
    return fig


def _grid_spacing(x, z):
    """Median lattice pitch [mm]: the length one reference arrow stands for."""
    dx = np.median(np.diff(np.unique(np.round(x, 6))))
    dz = np.median(np.diff(np.unique(np.round(z, 6))))
    gs = np.nanmin([dx, dz])
    return float(gs) if np.isfinite(gs) and gs > 0 else 0.25
