"""Figures that compare conditions: the _summary set.

Ports viz.plot_nrmse_vs_speed, plot_valid_fraction, plot_est_vs_true and
plot_bland_altman. The first two read the pooled per-slice metrics
(_summary/metrics_slice_all.csv); the last two read the pooled points table
(every dataset's points.csv concatenated).
"""

from __future__ import annotations

import numpy as np

from . import style as st

# Fixed seed: the jitter in plot_valid_fraction is cosmetic, and a figure
# should not change every time it is re-rendered.
_RNG_SEED = 0


def plot_nrmse_vs_speed(msl, title=""):
    """Headline figure: accuracy against flow speed.

    Bars are the across-slice mean NRMSE; the error bar is the across-slice
    std, i.e. how reproducible the accuracy is along the vessel rather than
    how large the residual is within one slice.
    """
    names, speed = st.order_datasets(msl["dataset"])
    n = len(names)
    cols = ("nrmse_raw", "nrmse_sg")
    M = np.full((n, 2), np.nan)
    E = np.full((n, 2), np.nan)
    for k, nm in enumerate(names):
        Q = msl[msl["dataset"] == nm]
        for j, c in enumerate(cols):
            M[k, j] = np.nanmean(100 * Q[c]) if Q[c].notna().any() else np.nan
            E[k, j] = np.nanstd(100 * Q[c], ddof=1) if Q[c].notna().sum() > 1 else 0.0

    fig = st.new_figure(st.FIG_SINGLE)
    ax = fig.add_subplot(111)
    x = np.arange(n, dtype=float)
    for j, (off, colour, lbl) in enumerate(
            [(-0.18, st.C_RAW, "raw profile"), (0.18, st.C_SG, "Savitzky-Golay")]):
        ax.bar(x + off, M[:, j], width=0.32, color=colour, edgecolor="none",
               label=lbl, zorder=2)
        ax.errorbar(x + off, M[:, j], E[:, j], fmt="none", ecolor="k",
                    elinewidth=1.6, capsize=8, capthick=1.6, zorder=3)

    st.style_axes(ax)
    ax.set_xticks(x, st.speed_labels(names, speed))
    st.label_axes(ax, "peak velocity",
                  "NRMSE [%]   (mean $\\pm$ std across slices)",
                  f"{title}  through-plane accuracy vs flow speed".strip())
    st.legend(ax, loc="upper left")

    dead = ~np.isfinite(M).any(axis=1)
    for k in np.flatnonzero(dead):
        ax.text(x[k], 0, "no valid points", ha="center", va="bottom",
                style="italic", color="0.45", fontsize=st.FS_ANNOT)
    return fig


def plot_valid_fraction(msl, title=""):
    """How much of the lumen the estimator could answer, per condition.

    Companion to plot_nrmse_vs_speed and has to be read with it: NRMSE is
    computed only over points that produced an estimate, so a condition can
    post a flattering NRMSE precisely because the hard points dropped out.
    """
    names, speed = st.order_datasets(msl["dataset"])
    n = len(names)
    rng = np.random.default_rng(_RNG_SEED)

    fig = st.new_figure(st.FIG_SINGLE)
    ax = fig.add_subplot(111)
    for k, nm in enumerate(names):
        v = 100 * msl.loc[msl["dataset"] == nm, "valid_frac"].to_numpy(float)
        jitter = (rng.random(v.size) - 0.5) * 0.28
        ax.scatter(k + jitter, v, s=44, color=st.C_POINT2, alpha=0.65,
                   edgecolors="none", zorder=2)
        m = np.nanmean(v)
        s = np.nanstd(v, ddof=1) if np.isfinite(v).sum() > 1 else 0.0
        ax.errorbar(k, m, s, color=st.C_RAW, lw=2.4, capsize=12, capthick=2.0,
                    marker="o", ms=11, zorder=3)
        ax.text(k, 102, f"{m:.0f}%", ha="center", fontweight="bold",
                color="0.25", fontsize=st.FS_ANNOT)

    st.style_axes(ax)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(0, 110)
    ax.set_xticks(np.arange(n), st.speed_labels(names, speed))
    ax.set_yticks(np.arange(0, 101, 20))
    st.label_axes(ax, "peak velocity", "points with a valid $v_y$ [%]",
                  f"{title}  estimator yield per slice".strip())
    return fig


def plot_est_vs_true(P, title="", hide_outliers=True):
    """Estimated against true v_y, all conditions on one axis.

    Departure from the diagonal is the bias, vertical scatter is the
    precision, and the velocity ceiling shows up directly: points stop
    appearing above the speed the lag scan can still resolve.
    """
    P = _paired(P, hide_outliers)
    names, speed = st.order_datasets(P["dataset"])
    cols = st.dataset_colors(len(names))

    fig = st.new_figure(st.FIG_SQUARE)
    ax = fig.add_subplot(111)
    lo = min(0.0, float(np.nanmin([P["Vy_mms"].min(), P["VyTrue_mms"].min()]))) * 1.05
    hi = float(np.nanmax([P["Vy_mms"].max(), P["VyTrue_mms"].max()])) * 1.05

    for k, nm in enumerate(names):
        Q = P[P["dataset"] == nm]
        r = _corr_safe(Q["VyTrue_mms"].to_numpy(float), Q["Vy_mms"].to_numpy(float))
        ax.scatter(Q["VyTrue_mms"], Q["Vy_mms"], s=26, color=cols[k], alpha=0.35,
                   edgecolors="none", zorder=2,
                   label=f"{st.speed_labels([nm], speed[k:k + 1])[0]}  "
                         f"(n={len(Q)}, r={r:.3f})")
    ax.plot([lo, hi], [lo, hi], "--k", lw=2.0, zorder=3, label="1:1")

    ax.set_aspect("equal")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    st.style_axes(ax)
    st.label_axes(ax, "true $v_y$ [mm/s]", "estimated $v_y$ [mm/s]",
                  f"{title}  estimate vs truth".strip())
    st.legend(ax, loc="upper left")
    return fig


def plot_bland_altman(P, title="", hide_outliers=True, relative=False,
                      estimate_col="Vy_mms"):
    """Agreement between estimate and truth.

    A trend in the cloud means the error scales with velocity rather than
    being a fixed offset -- the distinction an RMSE number alone cannot make.
    estimate_col picks which estimate to score (Vy_raw_mms vs Vy_sg_mms), the
    same knob main() uses to draw the raw and post-SG versions.
    """
    P = _paired(P, hide_outliers, estimate_col)
    avg = ((P[estimate_col] + P["VyTrue_mms"]) / 2).to_numpy(float)
    dif = (P[estimate_col] - P["VyTrue_mms"]).to_numpy(float)
    unit = "mm/s"
    if relative:
        ref = np.where(avg == 0, np.nan, avg)
        dif = 100 * dif / ref
        unit = "%"
    ok = np.isfinite(avg) & np.isfinite(dif)
    avg, dif, P = avg[ok], dif[ok], P[ok]

    bias = float(np.mean(dif))
    sd = float(np.std(dif, ddof=1))
    loa = (bias - 1.96 * sd, bias + 1.96 * sd)

    names, speed = st.order_datasets(P["dataset"])
    cols = st.dataset_colors(len(names))

    fig = st.new_figure(st.FIG_SINGLE)
    ax = fig.add_subplot(111)
    for k, nm in enumerate(names):
        s = (P["dataset"] == nm).to_numpy()
        ax.scatter(avg[s], dif[s], s=26, color=cols[k], alpha=0.35,
                   edgecolors="none", zorder=2,
                   label=st.speed_labels([nm], speed[k:k + 1])[0])

    xr = (float(avg.min()), float(avg.max()))
    ax.plot(xr, [bias, bias], "-", color=st.C_LINE, lw=2.4, zorder=3)
    for y in loa:
        ax.plot(xr, [y, y], "--", color=st.C_LOA, lw=2.0, zorder=3)
    # Anchored to the axes, not to the data: the rightmost point sits on the
    # frame, so a data-anchored label hangs off the edge.
    lbl = ax.get_yaxis_transform()
    ax.text(0.99, bias, f"bias {bias:.1f}", transform=lbl, ha="right",
            va="bottom", color=st.C_LINE, fontsize=st.FS_ANNOT,
            fontweight="bold")
    ax.text(0.99, loa[1], f"+1.96sd {loa[1]:.1f}", transform=lbl, ha="right",
            va="bottom", color=st.C_LOA, fontsize=st.FS_ANNOT)
    ax.text(0.99, loa[0], f"-1.96sd {loa[0]:.1f}", transform=lbl, ha="right",
            va="top", color=st.C_LOA, fontsize=st.FS_ANNOT)

    st.style_axes(ax)
    st.label_axes(ax, "mean of estimate and truth [mm/s]",
                  f"estimate - truth [{unit}]",
                  f"{title}  Bland-Altman (n = {dif.size})".strip())
    st.legend(ax, loc="lower left")
    return fig


def plot_bland_altman_slices(T, value="flow", title=""):
    """Bland-Altman with one point per slice and seed, not per grid cell.

    The per-cell version drowned in ~60k points; agreement is a per-slice
    question anyway -- does the reconstructed slice deliver the right flow
    (or mean velocity)? T needs columns dataset, est, tru; typically 41
    slices x seeds x speeds, a few hundred points, colored by speed.
    value only labels the axes: "flow" [ml/min] or "velocity" [mm/s].
    """
    unit = "ml/min" if value == "flow" else "mm/s"
    est = T["est"].to_numpy(float)
    tru = T["tru"].to_numpy(float)
    ok = np.isfinite(est) & np.isfinite(tru)
    T, est, tru = T[ok], est[ok], tru[ok]
    avg = (est + tru) / 2
    dif = est - tru

    bias = float(np.mean(dif))
    sd = float(np.std(dif, ddof=1))
    loa = (bias - 1.96 * sd, bias + 1.96 * sd)

    names, speed = st.order_datasets(T["dataset"])
    cols = st.dataset_colors(len(names))

    fig = st.new_figure(st.FIG_SINGLE)
    ax = fig.add_subplot(111)
    for k, nm in enumerate(names):
        sel = (T["dataset"] == nm).to_numpy()
        ax.scatter(avg[sel], dif[sel], s=42, color=cols[k], alpha=0.75,
                   edgecolors="none", zorder=2,
                   label=st.speed_labels([nm], speed[k:k + 1])[0])
    xr = (float(avg.min()), float(avg.max()))
    ax.plot(xr, [bias, bias], "-", color=st.C_LINE, lw=2.4, zorder=3)
    for yv in loa:
        ax.plot(xr, [yv, yv], "--", color=st.C_LOA, lw=2.0, zorder=3)
    lbl = ax.get_yaxis_transform()
    ax.text(0.99, bias, f"bias {bias:.1f}", transform=lbl, ha="right",
            va="bottom", color=st.C_LINE, fontsize=st.FS_ANNOT,
            fontweight="bold")
    ax.text(0.99, loa[1], f"+1.96sd {loa[1]:.1f}", transform=lbl, ha="right",
            va="bottom", color=st.C_LOA, fontsize=st.FS_ANNOT)
    ax.text(0.99, loa[0], f"-1.96sd {loa[0]:.1f}", transform=lbl, ha="right",
            va="top", color=st.C_LOA, fontsize=st.FS_ANNOT)
    st.style_axes(ax)
    st.label_axes(ax, f"mean of estimate and truth [{unit}]",
                  f"estimate - truth [{unit}]",
                  f"{title}  Bland-Altman, one point per slice "
                  f"(n = {dif.size})".strip())
    st.legend(ax, loc="lower left")
    return fig


def _paired(P, hide_outliers, estimate_col="Vy_mms"):
    """Rows where estimate and truth are both finite (and not flagged).

    One realisation per slice, matching what main() feeds these figures
    (representative_seed_points); repeats would otherwise weigh a slice by how
    many times it happened to be re-run.
    """
    if estimate_col not in P.columns:
        raise ValueError(f"points table has no {estimate_col} column")
    P = st.representative_seed(P)
    m = np.isfinite(P[estimate_col]) & np.isfinite(P["VyTrue_mms"])
    if hide_outliers and "outlier" in P.columns:
        m &= ~P["outlier"].astype(bool)
    out = P[m]
    if out.empty:
        raise ValueError("no finite paired points in the points table")
    return out


def _corr_safe(a, b):
    if a.size < 2 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])
