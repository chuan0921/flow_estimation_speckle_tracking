function o = default_opts(o)
%DEFAULT_OPTS Fill src.estimate_slice() options with defaults.
def.grid_step = 0.25e-3;
def.roi_half = 0.5e-3;
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
def.drift_win = 3;
def.fine_on = false;      % fine drift stage (two-grain window)
def.fine_roi_half_x = 0.6e-3;   % half-extent: window holds two lateral grains
def.fine_roi_half_z = 0.15e-3;  % one axial grain
def.fine_win = 2;         % fine search half-width around the anchor [px]
def.vy_fine_on = true;    % local speckle-ROI vy rescan, also used at zero XZ shift
def.vy_fine_roi_half_x = 0.3e-3;  % one lateral grain half-extent
def.vy_fine_roi_half_z = 0.15e-3; % one axial grain half-extent
def.vy_spatial_half = 1;  % weighted CC filter over centre +/- one speckle ROI
def.scan_reverse = 'auto';  % 'auto' | 'on' | 'off'
def.rev_y_mm = 3.5;         % auto: scan reverse when slice_pos >= this
def.rev_vy_max = 0.15;      % max plausible reverse speed [m/s]
def.z_sig_rev = 2.2;        % significance gates for the reverse peak
def.z_broad = 1.5;
def.z_half = 1.0;
def.nb_half = 5;
def.prom_frac = 0.3;        % peak prominence: required drop within prom_w
def.prom_w = 12;            % samples to each side; edge-truncated = reject
def.prom_tail_min = 1200;   % forward: only gate peaks beyond this lag
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
