"""Figures for the flow-estimation pipeline, rendered from the CSVs.

main() (MATLAB) computes and writes results/<dataset>/*.csv; everything in
here reads those CSVs and writes PNGs. No figure reaches back into MATLAB, so
plots can be restyled and re-rendered without re-running the estimation.

    from viz import plot_profile, style
    fig = plot_profile(pd.read_csv('results/<ds>/profile.csv'), title='<ds>')
    style.save(fig, 'out.png')

Per dataset            from                              figure
  plot_profile           profile.csv                       profile.png
  plot_error_vs_radius   profile.csv                       error_vs_radius.png
  plot_vessel_summary    metrics_slice.csv + vessel_regions.csv
                                                           vessel_summary.png
  plot_cc_peak_map       points.csv                        cc_peak_map.png
  plot_vector_field      points.csv (one slice)            vector_field.png
Across datasets        from                              figure
  plot_vessel_summary    _summary/metrics_slice_all.csv    <name>_vessel_summary.png
  plot_valid_fraction    _summary/metrics_slice_all.csv    valid_fraction.png
  plot_nrmse_vs_speed    _summary/metrics_slice_all.csv    nrmse_vs_speed.png
  plot_est_vs_true       pooled points.csv                 est_vs_true.png
  plot_bland_altman      pooled points.csv                 bland_altman_{raw,sg}.png

Figures that show one point per ROI keep a single realisation per slice
(style.representative_seed), matching representative_seed_points in main.m.
"""

from . import style
from .cc_peak_map import plot_cc_peak_map
from .profile import plot_error_vs_radius, plot_profile
from .summary import (plot_bland_altman, plot_est_vs_true,
                      plot_nrmse_vs_speed, plot_valid_fraction)
from .vector_field import plot_vector_field
from .vessel_summary import plot_vessel_summary

__all__ = [
    "style",
    "plot_profile", "plot_error_vs_radius", "plot_vessel_summary",
    "plot_cc_peak_map", "plot_vector_field",
    "plot_nrmse_vs_speed", "plot_valid_fraction", "plot_est_vs_true",
    "plot_bland_altman",
]
