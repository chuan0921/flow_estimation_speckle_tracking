function o = default_opts(o)
%DEFAULT_OPTS Fill src.estimate_slice() options with defaults.
def.grid_step = 0.25e-3;
% One speckle grain: lateral FWHM is 0.507 mm, and the transit peak's
% usable height is the fraction of the window arriving in sync, so the
% window must not span more velocity than the survival tolerance allows
% (relative spread <= 2*sigma_y/d_row). Halving from 1.0 mm moved the
% measurable boundary from r/R 0.65 to 0.80 and cut the straight-vessel
% full-lumen error from 20% to 7%; below one grain there is no speckle
% left to correlate.
def.roi_half = 0.25e-3;
% An ROI is sampled when the echo energy outside the lumen is at most this
% fraction of the energy inside it (src.geometry). It replaces the old fixed
% 0.9*R cut, so the lattice reaches the wall instead of stopping short of it.
def.roi_out_energy_max = 0.25;
def.lag_min = 20;
def.lag_max = 4000;
def.n_lag_scan = 240;
% Empty means "derive from the acquisition" in src.estimate_slice: the row1
% self-correlation is bounded by the elevation beam (2*sigma_y), not by the
% row spacing, so the ladder is built in units of the beam crossing time.
% Set it explicitly only to override that for an experiment.
def.inplane_lags = [];
% How far a scatterer may travel through the beam, in beam widths, before
% its row1 speckle is written off. 1.0 leaves residual correlation ~exp(-1)
% = 0.37, still well clear of cc_abs_min; measured breakdown on this dataset
% sits between 1.24 and 2.16 beam widths.
def.inplane_survival = 1.0;
% 'auto' follows the phantom: analytic straight pipe is off, CFD is on.
% Logical true/false provides an explicit override for ablation tests.
def.compensate_xz = 'auto';
def.inplane_win = 3;
% Frame pairs averaged per block match, spread evenly over the ensemble. The
% sub-pixel error sets a floor of 0.5 px * d_row/(2*sigma_y) = 1.1 px on the
% row2 shift and averaging is the only lever on it, but the gain saturates
% once the pairs stop overlapping: measured residual at 40 / 85 / 170 / 400
% pairs is 0.84 / 0.62 / 0.57 / 0.58 px on the axis, so 85 buys the whole
% improvement for 5% more runtime where 400 costs 2.4x. Worth having and no
% more than that -- it gains 12 ROIs of 767, all of them mid-radius, and
% leaves the near-wall residual untouched at 3.1 px.
def.inplane_pairs = 85;
def.cc_abs_min = 0.15;
def.cc_rel_min = 0.35;
def.cc_min = 0.5;           % accept only peaks that are actually locked. On
                            % the 5-seed run: 0.2 keeps 92% of ROIs at point
                            % RMSE 103 mm/s, 0.6 keeps 32% at 16 mm/s but
                            % empties the outer half of the lumen and drives
                            % the flow 13-27% low. 0.5 sits between them.
def.fft_min_lags = 24;    % shift groups with >= this many lags use FFT
def.transit_cc = 'envelope';  % 'coherent' correlates complex IQ instead of
                            % the envelope in every transit scan. On
                            % experimental data the inter-row envelope
                            % correlation collapses in the fast core (the
                            % clutter's brightness common mode swamps the
                            % ~0.25 coherent remnant), while the complex
                            % correlation still locks: phase-misaligned
                            % clutter self-cancels in the coherent sum.
                            % Simulation stays on 'envelope' (equivalent
                            % there, and all tuning was done on it).
def.drift_win = 3;
def.fine_on = false;      % fine drift stage (two-grain window)
def.fine_roi_half_x = 0.6e-3;   % half-extent: window holds two lateral grains
def.fine_roi_half_z = 0.15e-3;  % one axial grain
def.fine_win = 2;         % fine search half-width around the anchor [px]
% The applied XZ shift is vx(row1)*L*dt, which grows without bound in L even
% though the physical displacement cannot: if lag L is the true transit, the
% shift is (vx/vy)*d_row regardless of L. In-plane flow above this fraction
% of the through-plane speed does not occur in these vessels, so shifts are
% clamped at xz_shift_max_ratio * d_row. Without the clamp the long-lag tail
% of a scan walks the row2 window up to 76 px away -- four window widths --
% onto unrelated speckle, and the reverse fallback manufactured 36 confident
% wrong answers per three slices out of exactly that.
def.xz_shift_max_ratio = 0.25;
def.vy_fine_on = true;    % local speckle-ROI vy rescan, also used at zero XZ shift
def.vy_fine_roi_half_x = 0.3e-3;  % one lateral grain half-extent
def.vy_fine_roi_half_z = 0.15e-3; % one axial grain half-extent
% Stage B rescans vy with the speckle window; this is how many lateral
% neighbour sub-windows join a curve-level weighted average around it.
% Zero. The neighbours sit 0.65 mm away at a different radius, and averaging
% their CC curves drags the peak exactly the way the stage A window mixture
% did -- it was this pipeline's largest remaining bias once the window
% shrank: switching 1 -> 0 took the straight-vessel full-lumen NRMSE from
% 5.4% to 2.7% (80 cm/s) and 7.3% to 3.2% (20 cm/s), flow from +2.9% to
% +1.3%, at identical coverage, and left the stenosis unchanged. The
% remaining ring bias then matches the physical window-mixture model, so
% the averaging was the missing mechanism, not a noise suppressor.
def.vy_spatial_half = 0;
def.scan_reverse = 'auto';  % auto follows phantom, then rev_y_mm
def.rev_y_mm = 3.5;         % auto: scan reverse when slice_pos >= this
def.rev_vy_max = 0.15;      % max plausible reverse speed [m/s]
def.z_sig_rev = 2.2;        % significance gates for the reverse peak
def.z_broad = 1.5;
def.z_half = 1.0;
def.nb_half = 5;
def.prom_frac = 0.3;        % peak prominence: required drop within prom_w
def.prom_w = 12;            % samples to each side; edge-truncated = reject.
                            % Also one of the reverse scan's significance
                            % gates, so it cannot be tuned for the forward
                            % tail alone: widening it to 20 during the
                            % long-lag experiments silently edge-rejected a
                            % quarter of the reverse lag grid.
def.prom_tail_min = 1200;   % forward: only gate peaks beyond this lag
def.prom_abs_min = 0;       % minimum absolute peak height above the curve
                            % median. The relative prominence test cannot
                            % reject a flat curve with micro-ripples (its
                            % 20%-of-nothing is satisfied by noise); static
                            % wall-texture plateaus in experimental data
                            % peak 0.0002-0.004 above their own median
                            % while genuine transit peaks sit 0.02-0.3
                            % above. 0 = off (simulation default).
def.wall_margin_mm = 0; % Wall clearance [mm]: exclude ROIs whose centre is
                        % within this distance of the local wall, so no
                        % window straddles it (0.25 = the ROI half-width).
                        % 0 = off (historical lattice-to-the-wall).
def.sharp_band = [];    % Peak-shape band [lo hi] for the transit peak's
                        % dimensionless sharpness curv*lag^2/cc. Empty =
                        % gate off (historical behaviour). A genuine peak
                        % sits at (d_row/sigma_y)^2 ~= 20 by beam geometry;
                        % [10 45] rejects frozen-speckle plateaus and
                        % micro-ripple spikes that clear the cc gate with a
                        % shapeless curve. The lower edge hugs the honest
                        % floor (straight-pipe p05 ~= 10.4): the 8-10 alley
                        % holds mixture-widened wall liars, and clearing it
                        % improved all seven datasets (straight 5-seed
                        % finals another -4..-9%, sten-20 upstream
                        % 11.3->9.5%) for 1-3pp of yield.
% Tissue localization (src.localize_scan). Margin in lumen radii kept
% clear of the mask (blood moves, PSF smears ~0.5 mm past the wall),
% spatial search range in scan steps, and the frame the snapshot is
% taken from. The ref is ALWAYS row2 -- the sign convention lives in
% localize_scan's help text and must not be reparameterised.
def.loc_margin = 1.3;
def.loc_range = 8;
def.loc_frame = 101;
def.loc_search_mm = 0;  % In-plane block-search half-width [mm]. 0 = plain
                        % same-position NCC (correct while the probe only
                        % translates -- every simulation). >0 = classic
                        % block matching for experimental scans where
                        % probe pressure or an uneven surface shifts the
                        % tissue in-plane between the paired planes: each
                        % tissue block keeps its best ZNCC within +/- this
                        % many mm in x and z per slice offset.
def.loc_block_mm = 2;   % Block edge [mm] in block mode.
def.loc_blocks = 3;     % Tissue blocks; slice-offset peak is the median
                        % across blocks and the majority must agree.
def.nmt_thresh = 2;
def.nmt_eps_vy = 5e-3;
def.nmt_eps_ip = 1.5e-3;
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(o, fn{k}) || isempty(o.(fn{k}))
        o.(fn{k}) = def.(fn{k});
    end
end
end
