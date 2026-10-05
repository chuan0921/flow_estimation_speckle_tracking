#!/usr/bin/env python3
"""Render every figure for a results folder written by main().

    python3 make_figures.py results/sim_stenosis50_80cms   # one dataset
    python3 make_figures.py results                        # all + _summary
    python3 make_figures.py results --outdir /tmp/check    # somewhere else
    python3 make_figures.py --speeds results/sim_v0*_5seed/sim_v0*  # one metric per PNG

Figures land in <dataset>/figures/ and results/_summary/figures/ unless
--outdir says otherwise. Each figure is a pure function of a CSV, so
re-rendering never needs MATLAB or a re-run of the estimation.
"""

from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import pandas as pd

import viz
from viz import style

SUMMARY_DIR = "_summary"


def dataset_figures(ds_dir, out_dir=None, title=None, slice_idx=None,
                    roi_slices=None, no_truth=False, z0=None):
    """Per-dataset figures. Returns the paths written."""
    ds_dir = os.path.abspath(ds_dir)
    out_dir = out_dir or os.path.join(ds_dir, "figures")
    title = os.path.basename(ds_dir) if title is None else title
    written = []

    prof_csv = os.path.join(ds_dir, "profile.csv")
    if os.path.isfile(prof_csv):
        prof = pd.read_csv(prof_csv)
        if "diameter_pos_mm" in prof:
            for idx in sorted(prof["slice_idx"].unique()):
                one = prof[prof["slice_idx"] == idx]
                written.append(style.save(
                    viz.plot_profile(one, title=title),
                    os.path.join(out_dir, "bfvp", f"slice_{int(idx):03d}.png")))
        else:
            print(f"skip: radial profile.csv in {ds_dir}; re-run main for BFVP")

    msl_csv = os.path.join(ds_dir, "metrics_slice.csv")
    if os.path.isfile(msl_csv):
        written.append(style.save(
            viz.plot_vessel_summary(pd.read_csv(msl_csv),
                                    regions=_read(ds_dir, "vessel_regions.csv"),
                                    title=title),
            os.path.join(out_dir, "vessel_summary.png")))

    pts_csv = os.path.join(ds_dir, "points.csv")
    if os.path.isfile(pts_csv):
        P = pd.read_csv(pts_csv)
        written.append(style.save(viz.plot_cc_peak_map(P, title=title),
                                  os.path.join(out_dir, "cc_peak_map.png")))
        # One slice carries the vector field: the middle one by default, the
        # same choice main() made when it still drew this figure itself.
        slices = sorted(P["slice_idx"].unique()) if "slice_idx" in P else [None]
        pick = slice_idx if slice_idx is not None else slices[len(slices) // 2]
        Tm = P[P["slice_idx"] == pick] if pick is not None else P
        if Tm["Vy_mms"].notna().any():
            written.append(style.save(
                viz.plot_vector_field(Tm, title=f"{title} slice {pick:g}",
                                      truth=True),
                os.path.join(out_dir, "vector_field.png")))
        written += roi_peak_figures(P, out_dir, title, roi_slices, no_truth, z0)
    if not written:
        print(f"skip: no CSV to plot in {ds_dir}")
    return written


def roi_peak_figures(P, out_dir, title, want=None, no_truth=False, z0=None):
    """Full ROI peak-lag map plus the three 3x3 zooms, per requested slice.

    Defaults to one slice per flow segment rather than all of them: the full map
    carries per-tile text and is a large canvas, and three slices already answer
    "does the lag behave differently upstream, in the throat and downstream".
    Pass want to name the slices instead, for a specific one worth reading.
    """
    written = []
    color = "velocity" if no_truth else "error"
    for idx in (want if want else _representative_slices(P)):
        d = viz.roi_peak_frame(P, idx)
        lim = viz.velocity_limit(d) if no_truth else viz.lag_error_limit(d)
        stem = f"slice_{int(idx):03d}"
        head = f"{title} slice {idx:g}"
        written.append(style.save(
            viz.plot_roi_peak_map(d, title=head, lim=lim, color=color, z0=z0),
            os.path.join(out_dir, "roi_peak", f"{stem}_full.png")))
        for name, target in viz.ZOOMS:
            written.append(style.save(
                viz.plot_roi_peak_3x3(d, target, name, title=head, lim=lim,
                                      color=color, z0=z0),
                os.path.join(out_dir, "roi_peak", f"{stem}_3x3_{name}.png")))
    return written


def _representative_slices(P):
    """One slice per flow segment: the one nearest that segment's midpoint."""
    if "vessel_region" not in P.columns:
        slices = sorted(P["slice_idx"].unique())
        return [slices[len(slices) // 2]]
    d = P[["slice_idx", "slice_pos_mm", "vessel_region"]].drop_duplicates()
    d = d.assign(seg=[style.segment(v)[0] for v in d["vessel_region"]])
    picked = []
    for seg, g in d.groupby("seg", sort=False):
        mid = g["slice_pos_mm"].median()
        row = g.iloc[(g["slice_pos_mm"] - mid).abs().to_numpy().argmin()]
        picked.append((seg, row["slice_idx"], row["slice_pos_mm"]))
    picked.sort(key=lambda t: t[2])
    print("roi_peak: " + ", ".join(
        f"{seg} -> slice {int(i)} (y={y:+.1f} mm)" for seg, i, y in picked))
    return [i for _, i, _ in picked]


def speed_figures(dataset_dirs, out_dir, title=""):
    """One PNG per metric, every speed condition drawn in it.

    The transpose of vessel_summary: that figure stacks a run's metrics, this
    one lines a metric up across runs. Written to speed_<metric>.png so a
    single metric can be opened without the others around it.
    """
    P = viz.collect_speeds(dataset_dirs)
    names, _ = style.order_datasets(P["dataset"])
    print(f"speeds: {', '.join(str(n) for n in names)}")
    return [style.save(viz.plot_speed_metric(P, key, title=title),
                       os.path.join(out_dir, f"speed_{key}.png"))
            for key in viz.metric_keys()]


def summary_figures(results_dir, out_dir=None, title=""):
    """Across-dataset figures. Returns the paths written."""
    results_dir = os.path.abspath(results_dir)
    out_dir = out_dir or os.path.join(results_dir, SUMMARY_DIR, "figures")
    written = []

    msl_csv = os.path.join(results_dir, SUMMARY_DIR, "metrics_slice_all.csv")
    if os.path.isfile(msl_csv):
        msl = pd.read_csv(msl_csv)
        regions = _read(os.path.join(results_dir, SUMMARY_DIR),
                        "vessel_regions_all.csv")
        for nm in pd.unique(msl["dataset"]):
            reg = None if regions is None else regions[regions["dataset"] == nm]
            written.append(style.save(
                viz.plot_vessel_summary(msl[msl["dataset"] == nm], regions=reg,
                                        title=str(nm)),
                os.path.join(out_dir, f"{nm}_vessel_summary.png")))
        written.append(style.save(viz.plot_valid_fraction(msl, title=title),
                                  os.path.join(out_dir, "valid_fraction.png")))
        # Only meaningful with more than one condition to line up.
        if msl["dataset"].nunique() > 1:
            written.append(style.save(viz.plot_nrmse_vs_speed(msl, title=title),
                                      os.path.join(out_dir, "nrmse_vs_speed.png")))
    else:
        print(f"skip: no {SUMMARY_DIR}/metrics_slice_all.csv in {results_dir}")

    # main() never writes the pooled points table, so rebuild it from the
    # per-dataset files -- same rows, just never persisted in one place.
    parts = [pd.read_csv(f) for f in
             sorted(glob.glob(os.path.join(results_dir, "*", "points.csv")))]
    if parts:
        P = pd.concat(parts, ignore_index=True)
        written.append(style.save(viz.plot_est_vs_true(P, title=title),
                                  os.path.join(out_dir, "est_vs_true.png")))
    else:
        print(f"skip: no <dataset>/points.csv under {results_dir}")

    # Bland-Altman at the per-slice scale: the per-cell clouds (60k points)
    # were unreadable, and agreement is a per-slice question -- does each
    # reconstructed slice deliver the right flow and mean velocity.
    fparts = [pd.read_csv(f) for f in
              sorted(glob.glob(os.path.join(results_dir, "*", "field.csv")))]
    if fparts:
        F = pd.concat(fparts, ignore_index=True)
        F = F[F["sample_type"] == "grid"]
        ok = np.isfinite(F["Vy_sg_mms"]) & np.isfinite(F["VyTrue_mms"]) \
            & np.isfinite(F["VyPlaneTrue_mms"])
        F = F[ok]
        F["fe"] = F["Vy_sg_mms"] * F["cell_area_mm2"]
        F["ft"] = F["VyPlaneTrue_mms"] * F["cell_area_mm2"]
        g = F.groupby(["dataset", "slice_idx", "seed_base"])
        agg = g.agg(fe=("fe", "sum"), ft=("ft", "sum"),
                    ve=("Vy_sg_mms", "mean"),
                    vt=("VyTrue_mms", "mean")).reset_index()
        flow = agg.rename(columns={"fe": "est", "ft": "tru"})
        flow[["est", "tru"]] *= 60.0 / 1000.0        # mm^3/s -> ml/min
        written.append(style.save(
            viz.plot_bland_altman_slices(flow, value="flow",
                                         title="Flow per slice"),
            os.path.join(out_dir, "bland_altman_flow.png")))
        vel = agg.rename(columns={"ve": "est", "vt": "tru"})
        written.append(style.save(
            viz.plot_bland_altman_slices(vel, value="velocity",
                                         title="Mean velocity per slice"),
            os.path.join(out_dir, "bland_altman_velocity.png")))
    else:
        print(f"skip: no <dataset>/field.csv under {results_dir}")
    return written


def _read(dir_path, name):
    """Optional CSV: the table, or None when the pipeline did not write it."""
    path = os.path.join(dir_path, name)
    return pd.read_csv(path) if os.path.isfile(path) else None


def is_dataset_dir(path):
    return any(os.path.isfile(os.path.join(path, f))
               for f in ("points.csv", "profile.csv", "metrics_slice.csv"))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", nargs="+",
                    help="results/<dataset>, or the results root; with "
                         "--speeds, one dataset folder per condition")
    ap.add_argument("--speeds", action="store_true",
                    help="combine the given dataset folders into one PNG per "
                         "metric instead of per-dataset figures")
    ap.add_argument("--outdir", default=None,
                    help="default: <dataset>/figures and _summary/figures")
    ap.add_argument("--title", default=None, help="default: the folder name")
    ap.add_argument("--slice", type=float, default=None,
                    help="slice_idx for the vector field (default: the middle one)")
    ap.add_argument("--roi-slice", type=float, action="append", default=None,
                    metavar="IDX",
                    help="slice_idx for the ROI peak-lag figures; repeatable. "
                         "Default: one slice per flow segment")
    ap.add_argument("--z0", type=float, default=None,
                    help="absolute depth of the vessel centre [mm]; when "
                         "given, ROI-peak z axes show image depth instead "
                         "of vessel-centred z")
    ap.add_argument("--no-truth", action="store_true",
                    help="experimental data: no ground truth exists, so the "
                         "ROI peak tiles are coloured by estimated velocity "
                         "and the true-lag line is dropped")
    args = ap.parse_args()

    written = []
    if args.speeds:
        out = args.outdir or os.path.join(
            os.path.dirname(os.path.abspath(args.path[0])), "speed_figures")
        for p in speed_figures([os.path.abspath(d) for d in args.path], out,
                               args.title or ""):
            print(p)
        return

    if len(args.path) > 1:
        ap.error("pass one path, or use --speeds with several dataset folders")
    path = os.path.abspath(args.path[0])
    if is_dataset_dir(path):
        written += dataset_figures(path, args.outdir, args.title, args.slice,
                                   args.roi_slice, args.no_truth, args.z0)
    else:
        for ds in sorted(d for d in glob.glob(os.path.join(path, "*"))
                         if os.path.isdir(d) and is_dataset_dir(d)):
            sub = None if args.outdir is None else os.path.join(
                args.outdir, os.path.basename(ds))
            written += dataset_figures(ds, sub, args.title, args.slice,
                                       args.roi_slice, args.no_truth, args.z0)
        sub = None if args.outdir is None else os.path.join(args.outdir, SUMMARY_DIR)
        written += summary_figures(path, sub, args.title or "")

    for p in written:
        print(p)


if __name__ == "__main__":
    main()
