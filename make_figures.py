#!/usr/bin/env python3
"""Render every figure for a results folder written by main().

    python3 make_figures.py results/sim_stenosis50_80cms   # one dataset
    python3 make_figures.py results                        # all + _summary
    python3 make_figures.py results --outdir /tmp/check    # somewhere else

Figures land in <dataset>/figures/ and results/_summary/figures/ unless
--outdir says otherwise. Each figure is a pure function of a CSV, so
re-rendering never needs MATLAB or a re-run of the estimation.
"""

from __future__ import annotations

import argparse
import glob
import os

import pandas as pd

import viz
from viz import style

SUMMARY_DIR = "_summary"


def dataset_figures(ds_dir, out_dir=None, title=None, slice_idx=None):
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
    if not written:
        print(f"skip: no CSV to plot in {ds_dir}")
    return written


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
        for col, ba_title, png in _bland_altman_variants(P):
            written.append(style.save(
                viz.plot_bland_altman(P, title=ba_title, estimate_col=col),
                os.path.join(out_dir, png)))
    else:
        print(f"skip: no <dataset>/points.csv under {results_dir}")
    return written


def _bland_altman_variants(P):
    """(column, title, filename) per estimate worth a Bland-Altman panel.

    main() draws the raw and post-SG estimates separately; tables written
    before those columns existed get the single Vy_mms version.
    """
    want = [("Vy_raw_mms", "Raw", "bland_altman_raw.png"),
            ("Vy_sg_mms", "After SG (2-D)", "bland_altman_sg.png")]
    have = [w for w in want if w[0] in P.columns]
    return have or [("Vy_mms", "", "bland_altman.png")]


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
    ap.add_argument("path", help="results/<dataset>, or the results root")
    ap.add_argument("--outdir", default=None,
                    help="default: <dataset>/figures and _summary/figures")
    ap.add_argument("--title", default=None, help="default: the folder name")
    ap.add_argument("--slice", type=float, default=None,
                    help="slice_idx for the vector field (default: the middle one)")
    args = ap.parse_args()

    path = os.path.abspath(args.path)
    written = []
    if is_dataset_dir(path):
        written += dataset_figures(path, args.outdir, args.title, args.slice)
    else:
        for ds in sorted(d for d in glob.glob(os.path.join(path, "*"))
                         if os.path.isdir(d) and is_dataset_dir(d)):
            sub = None if args.outdir is None else os.path.join(
                args.outdir, os.path.basename(ds))
            written += dataset_figures(ds, sub, args.title, args.slice)
        sub = None if args.outdir is None else os.path.join(args.outdir, SUMMARY_DIR)
        written += summary_figures(path, sub, args.title or "")

    for p in written:
        print(p)


if __name__ == "__main__":
    main()
