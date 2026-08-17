"""Centre-diameter BFVP figures made from profile.csv."""

from __future__ import annotations

import numpy as np
from matplotlib.ticker import MaxNLocator

from . import style as st


def plot_profile(prof, title="", band=True):
    """Plot one slice from the left wall (-R) to the right wall (+R)."""
    prof = st.representative_seed(prof)
    if "diameter_pos_mm" not in prof:
        raise ValueError("profile.csv is radial; re-run main for diameter BFVP")
    if prof["slice_idx"].nunique() != 1:
        raise ValueError("a BFVP is slice-specific; pass exactly one slice")
    prof = prof.sort_values("diameter_pos_mm")
    x = prof["diameter_pos_mm"]

    fig = st.new_figure(st.FIG_SINGLE)
    ax = fig.add_subplot(111)

    ax.axvline(0, color=st.C_LINE, lw=1.2, ls=":", zorder=1)
    ax.scatter(x, prof["vy_raw"], s=62, color=st.C_RAW, edgecolors="white",
               linewidths=0.6, zorder=5, label="Raw ROI estimate")
    ax.plot(x, prof["vy_fill"], "--", color=st.C_RAW, lw=2.0, zorder=2,
            label="No-slip interpolation")
    ax.plot(x, prof["vy_sg"], "-", color=st.C_SG, lw=3.0, zorder=3,
            label="Estimate after 2-D SG")
    ax.plot(x, prof["vy_true"], "--", color=st.C_TRUTH, lw=2.4, zorder=4,
            label="CFD truth")

    st.style_axes(ax)
    # Same three-segment naming the vessel summary shades with.
    region = (f" - {st.segment(prof['vessel_region'].iloc[0])[2]}"
              if "vessel_region" in prof else "")
    slice_idx = prof["slice_idx"].iloc[0]
    pos = prof["slice_pos_mm"].iloc[0]
    st.label_axes(ax, "Position along centre diameter x [mm]", "$v_y$ [mm/s]",
                  f"{title} | Slice {slice_idx:g} (y = {pos:.1f} mm)\n"
                  f"{region.strip(' -')} | Diameter BFVP".strip())
    radius = np.nanmax(np.abs(x))
    ax.set_xlim(-radius, radius)
    st.legend(ax, loc="lower left")
    return fig


def plot_error_vs_radius(prof, title="", relative=False):
    """Velocity error along one centre diameter.

    The upper panel is the raw and SG-filtered bias with a +/-1 std ribbon
    across slices; the lower panel counts how many slices actually measured
    each radius, which is what makes the near-wall bias interpretable -- the
    outer radii are both the hardest to track and the thinnest sampled.
    """
    prof = st.representative_seed(prof)
    coord = "diameter_pos_mm" if "diameter_pos_mm" in prof else "r_mm"
    e_raw = prof["vy_raw"] - prof["vy_true"]
    e_sg = prof["vy_sg"] - prof["vy_true"]
    unit = "mm/s"
    if relative:
        ref = prof["vy_true"].replace(0, np.nan)
        e_raw = 100 * e_raw / ref
        e_sg = 100 * e_sg / ref
        unit = "%"

    r, braw, sraw, _ = st.collapse(prof[coord], e_raw)
    _, bsg, ssg, _ = st.collapse(prof[coord], e_sg)
    # Slices that actually measured each radius. Counted as distinct slices,
    # not as rows, so repeats cannot inflate it. (The MATLAB version averaged
    # the 0/1 flag here, so its "slices" axis topped out at 1; a count is what
    # the panel claims to show.)
    nval = _slices_per_radius(prof, r, coord)

    fig = st.new_figure(st.FIG_STACKED)
    gs = fig.add_gridspec(3, 1)
    ax_b = fig.add_subplot(gs[0:2, 0])
    ax_n = fig.add_subplot(gs[2, 0], sharex=ax_b)

    ax_b.axhline(0, color="k", lw=1.5, zorder=2)
    st.band(ax_b, r, braw, sraw, st.C_RAW, zorder=1)
    st.band(ax_b, r, bsg, ssg, st.C_SG, zorder=1)
    ax_b.plot(r, braw, "-o", color=st.C_RAW, lw=2.2, ms=6, zorder=3, label="raw")
    ax_b.plot(r, bsg, "-", color=st.C_SG, lw=3.0, zorder=4,
              label="Savitzky-Golay")
    st.style_axes(ax_b)
    st.label_axes(ax_b, None, f"$v_y$ error [{unit}]",
                  f"{title}  bias $\\pm$ 1 std across slices".strip())
    st.legend(ax_b, loc="upper left")
    ax_b.tick_params(labelbottom=False)

    step = np.median(np.diff(r)) if r.size > 1 else 0.05
    ax_n.bar(r, nval, width=0.8 * step, color=st.C_BAR, edgecolor="none")
    st.style_axes(ax_n)
    xlabel = ("Position along centre diameter x [mm]"
              if coord == "diameter_pos_mm" else "radial position r [mm]")
    st.label_axes(ax_n, xlabel, "slices")
    ax_n.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
    ax_n.set_xlim(r.min() - 0.05, r.max() + 0.05)
    return fig


def _slices_per_radius(prof, radii, coord="r_mm"):
    """How many distinct slices measured each radius, aligned to `radii`."""
    ok = prof[np.isfinite(prof["vy_raw"])]
    key = "slice_idx" if "slice_idx" in prof.columns else coord
    n = ok.groupby(np.round(ok[coord], 6))[key].nunique()
    return n.reindex(np.round(radii, 6), fill_value=0).to_numpy(float)
