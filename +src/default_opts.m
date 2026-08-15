function o = default_opts(o)
%DEFAULT_OPTS Fill src.estimate_slice() options with defaults.
def.grid_step = 0.25e-3;
def.roi_half = 0.5e-3;
def.lag_min = 20;
def.lag_max = 4000;
def.n_lag_scan = 240;
def.inplane_lags = [15 30 60 120 240];
def.inplane_win = 3;
def.inplane_pairs = 40;
def.cc_abs_min = 0.15;
def.cc_rel_min = 0.35;
def.cc_min = 0.2;
def.fft_min_lags = 24;    % shift groups with >= this many lags use FFT
def.drift_win = 3;
def.fine_on = false;      % fine drift stage (two-grain window)
def.fine_roi_half_x = 0.6e-3;   % half-extent: window holds two lateral grains
def.fine_roi_half_z = 0.15e-3;  % one axial grain
def.fine_win = 2;         % fine search half-width around the anchor [px]
def.vy_fine_on = true;    % after xz compensation, rescan vy with speckle ROIs
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
