function T = estimate_slice(data_file, phantom, opts)
%ESTIMATE_SLICE Three-component flow estimation for one slice.
% The per-slice estimator called by main() for every slice of a run. Core
% primitives live alongside it in +src, plotting in +viz:
%     T = src.estimate_slice(data_file);                       % analytic truth
%     T = src.estimate_slice(data_file, src.phantom_cfd(cfd)); % CFD truth
%     viz.plot_vector_field(T);
%
% The second argument supplies the vessel geometry and the ground truth,
% and takes either of two forms:
%   []  or omitted   pure-simulation data (psf_sim): a straight pipe with a
%                    parabolic profile is built from the slice file itself
%                    via src.phantom_parabolic
%   a phantom        the interface returned by src.phantom_cfd or
%                    src.phantom_parabolic
% A raw CFD struct (one carrying grid_Uy) is still accepted and wrapped
% with src.phantom_cfd automatically.
%
% Per grid point:
%   0a  rough vy from a coarse uncompensated transit scan -> in-plane
%       lag cap (survival-bias guard 0.25 * 2*sigma_y / vy)
%   0b  in-plane vx/vz: chained multi-lag row1 self-correlation
%       (src.in_plane)
%   A   transit-lag scan with lag-dependent integer xz shift of the row2
%       patch (scan_transit_lags), cc_min peak (pick_peak)
%   B   coarse drift at the transit peak (drift_match); optional
%       speckle-size fine stage (refine_drift, opts.fine_on); then a
%       recompensated rescan
% Points the forward scan cannot answer fall back to a significance-gated
% negative-lag (reverse flow) scan, controlled by opts.scan_reverse.
% nmt outlier flag on the grid (flag only, estimates untouched).
%
% Requires sigma_y, k0, lambda in the slice file.
%
% See also MAIN, SRC.PHANTOM_PARABOLIC, SRC.PHANTOM_CFD.

% No path bootstrap here: reaching src.estimate_slice already means the
% project root resolves, so every other src.* call resolves the same way.

if nargin < 2
    phantom = [];
end
if nargin < 3
    opts = struct();
end
opts = src.default_opts(opts);
tp = src.through_plane();   % handles to the through-plane primitives

d = load(data_file, 'iq_row1', 'iq_row2', 'x', 'z', 'dx', 'dz', 'Nt', ...
    'dt', 'xc', 'zc', 'R', 'd_row', 'y_row1', 'y_row2', 'slice_pos', ...
    'sigma_y', 'k0', 'lambda', 'vmax', 'flow_direction');
ph = resolve_phantom(phantom, d);
IQ1 = d.iq_row1;                       % complex single, for vz phase
IQ2 = d.iq_row2;
E1 = abs(IQ1);                         % single envelopes throughout
E2 = abs(IQ2);
Nt = d.Nt;
Nz = numel(d.z);
Nx = numel(d.x);
dt = d.dt;
d_row = d.d_row;

rz = round(opts.roi_half / d.dz);
rx = round(opts.roi_half / d.dx);
[grid_x, grid_z, keep, GX, GZ, R_loc, ysp] = src.geometry(d, ph, opts); %#ok<ASGLU>
nP = numel(grid_x);
ix0v = round((grid_x - d.x(1)) / d.dx) + 1;
iz0v = round((grid_z - d.z(1)) / d.dz) + 1;
r_pt = hypot(grid_x - d.xc, grid_z - d.zc);

lag_max = min(opts.lag_max, Nt - 64);
lags_f = unique(round(logspace(log10(opts.lag_min), log10(lag_max), ...
    opts.n_lag_scan)));
lags_coarse = unique(round(logspace(log10(opts.lag_min), log10(lag_max), 40)));
lag_ip = opts.inplane_lags(opts.inplane_lags < Nt / 2);
zshift = zeros(size(lags_coarse));
rxf = round(opts.fine_roi_half_x / d.dx);
rzf = round(opts.fine_roi_half_z / d.dz);
rx_vy = round(opts.vy_fine_roi_half_x / d.dx);
rz_vy = round(opts.vy_fine_roi_half_z / d.dz);

% Reverse scan setup: negative log-spaced lags within the plausible
% reverse speed range, only for points where the forward scan fails.
scan_rev = strcmp(opts.scan_reverse, 'on') || ...
    (strcmp(opts.scan_reverse, 'auto') && ...
    d.slice_pos * 1e3 >= opts.rev_y_mm);
lags_r = [];
if scan_rev
    dtc_min = ceil(d_row / opts.rev_vy_max / dt);
    if dtc_min < lag_max
        lags_r = -flip(unique(round(logspace(log10(dtc_min), ...
            log10(lag_max), 80))));
    else
        scan_rev = false;
    end
end

vy_g = nan(nP, 1);
cc_g = nan(nP, 1);
vx_e = nan(nP, 1);
vz_e = nan(nP, 1);
cc_B = nan(nP, 1);
n_ip = zeros(nP, 1);
refined = false(nP, 1);
fine_ok = false(nP, 1);
is_rev = false(nP, 1);
t_stage = zeros(1, 5);         % rough / inplane / scanA / stageB / reverse

for p = 1:nP
    iz0 = iz0v(p);
    ix0 = ix0v(p);
    An = tp.normalize_patch(E1, iz0, ix0, rz, rx, Nt);

    % 0a: rough vy for the in-plane lag cap.
    t0 = tic;
    cc = tp.scan_transit_lags(An, E2, iz0, ix0, rz, rx, Nt, lags_coarse, ...
        zshift, zshift, Nz, Nx, opts.fft_min_lags);
    [vy_rough, ~] = tp.pick_peak(cc, lags_coarse, d_row, dt, 0.15);
    cap = max(lag_ip);
    if ~isnan(vy_rough)
        cap = max(0.25 * 2 * d.sigma_y / max(abs(vy_rough), 1e-3) / dt, ...
            lag_ip(1));
    end
    t_stage(1) = t_stage(1) + toc(t0);

    % 0b: chained multi-lag in-plane vx/vz.
    t0 = tic;
    [vx_e(p), vz_e(p), n_ip(p)] = src.in_plane(E1, IQ1, iz0, ix0, ...
        rz, rx, Nt, lag_ip, cap, dt, d, opts, Nz, Nx);
    t_stage(2) = t_stage(2) + toc(t0);

    vxc = vx_e(p);
    vzc = vz_e(p);
    if isnan(vxc), vxc = 0; end
    if isnan(vzc), vzc = 0; end

    % A: compensated transit scan.
    t0 = tic;
    oz = vzc * lags_f * dt / d.dz;
    ox = vxc * lags_f * dt / d.dx;
    cc = tp.scan_transit_lags(An, E2, iz0, ix0, rz, rx, Nt, lags_f, oz, ox, ...
        Nz, Nx, opts.fft_min_lags);
    [vy_g(p), cc_g(p)] = tp.pick_peak(cc, lags_f, d_row, dt, opts.cc_min, ...
        opts.prom_frac, opts.prom_w, opts.prom_tail_min);
    t_stage(3) = t_stage(3) + toc(t0);

    % Reverse fallback: significance-gated negative-lag scan for points
    % the forward scan could not answer.
    if isnan(vy_g(p))
        if scan_rev
            t0 = tic;
            oz = vzc * lags_r * dt / d.dz;
            ox = vxc * lags_r * dt / d.dx;
            [ccr, ccra, ccrb] = tp.scan_transit_lags(An, E2, iz0, ix0, ...
                rz, rx, Nt, lags_r, oz, ox, Nz, Nx, opts.fft_min_lags);
            [vy_r, cc_r] = tp.pick_significant(ccr, ccra, ccrb, lags_r, ...
                d_row, dt, opts);
            if ~isnan(vy_r)
                vy_g(p) = vy_r;
                cc_g(p) = cc_r;
                is_rev(p) = true;
            end
            t_stage(5) = t_stage(5) + toc(t0);
        end
        continue;
    end

    % B: drift at the transit peak, optional fine stage, rescan.
    t0 = tic;
    L_pk = round(d_row / (vy_g(p) * dt));
    L_pk = min(max(L_pk, lags_f(1)), lags_f(end));
    w2 = opts.drift_win + ceil(L_pk / 800);
    ctr_z = vzc * L_pk * dt / d.dz;
    ctr_x = vxc * L_pk * dt / d.dx;
    [ox_t, oz_t, pk_d] = tp.drift_match(An, E2, iz0, ix0, rz, rx, L_pk, ...
        ctr_z, ctr_x, w2, Nt, Nz, Nx);
    cc_B(p) = pk_d;
    if ~isnan(ox_t)
        if opts.fine_on
            % Fine stage: three one-grain window panel anchored at the
            % coarse drift (locally-normalized spatial correlation filter).
            [ox_f, oz_f] = tp.refine_drift(E1, E2, IQ1, IQ2, iz0, ix0, ...
                rzf, rxf, L_pk, round(oz_t), round(ox_t), ...
                opts.fine_win, Nt, Nz, Nx, d.k0, d.lambda, d.dz);
            if ~isnan(ox_f)
                ox_t = ox_f;
                oz_t = oz_f;
                fine_ok(p) = true;
            end
        end
        vx_r = ox_t * d.dx / (L_pk * dt);
        vz_r = oz_t * d.dz / (L_pk * dt);
        if fine_ok(p)
            vx_e(p) = vx_r;            % final in-plane from transit drift
            vz_e(p) = vz_r;
        end
        oz = vz_r * lags_f * dt / d.dz;
        ox = vx_r * lags_f * dt / d.dx;
        if opts.vy_fine_on
            cc = tp.scan_transit_lags_spatial(E1, E2, iz0, ix0, rz_vy, rx_vy, ...
                Nt, lags_f, oz, ox, Nz, Nx, opts.fft_min_lags, ...
                opts.vy_spatial_half);
        else
            cc = tp.scan_transit_lags(An, E2, iz0, ix0, rz, rx, Nt, lags_f, ...
                oz, ox, Nz, Nx, opts.fft_min_lags);
        end
        [vy_2, cc_2] = tp.pick_peak(cc, lags_f, d_row, dt, opts.cc_min, ...
            opts.prom_frac, opts.prom_w, opts.prom_tail_min);
        if ~isnan(vy_2)
            vy_g(p) = vy_2;
            cc_g(p) = cc_2;
            refined(p) = true;
        end
    end
    t_stage(4) = t_stage(4) + toc(t0);
end
fprintf(['v2 timing [s]: rough %.0f, inplane %.0f, scanA %.0f, ' ...
    'stageB %.0f, reverse %.0f\n'], t_stage);

% Ground truth averaged over the transit path. The phantom takes
% vessel-centred x,z and absolute y, all in metres, and returns m/s.
Xq = repmat(grid_x - d.xc, 1, 5);
Zq = repmat(grid_z - d.zc, 1, 5);
Yq = repmat(ysp, nP, 1);
[ux_path, uy_path, uz_path] = ph.velocity(Xq, Yq, Zq);
vx_true = mean(ux_path, 2);
vy_true = mean(uy_path, 2);
vz_true = mean(uz_path, 2);

% Reverse contiguity: a genuine recirculation pocket is spatially
% connected; a reverse claim with no reverse neighbor is demoted to a
% rejection.
if any(is_rev)
    Fr = false(size(GX)); Fr(keep) = is_rev;
    nb = conv2(double(Fr), ones(3), 'same') - double(Fr);
    lone = Fr & nb < 1;
    lone_k = lone(keep);
    vy_g(lone_k) = nan;
    cc_g(lone_k) = nan;
    is_rev(lone_k) = false;
end

% nmt outlier flag on the grid.
Fy = nan(size(GX)); Fy(keep) = vy_g;
Fx = nan(size(GX)); Fx(keep) = vx_e;
Fz = nan(size(GX)); Fz(keep) = vz_e;
out_g = src.nmt_outlier(Fy, opts.nmt_thresh, opts.nmt_eps_vy) | ...
    src.nmt_outlier(Fx, opts.nmt_thresh, opts.nmt_eps_ip) | ...
    src.nmt_outlier(Fz, opts.nmt_thresh, opts.nmt_eps_ip);
is_out = out_g(keep);
lag_peak_frames = d_row ./ (vy_g * dt);

% Positions are vessel-centred in x and z, matching r_mm and the phantom.
T = table((grid_x - d.xc) * 1e3, (grid_z - d.zc) * 1e3, r_pt * 1e3, ...
    vx_e * 1e3, vy_g * 1e3, vz_e * 1e3, ...
    vx_true * 1e3, vy_true * 1e3, vz_true * 1e3, ...
    cc_g, lag_peak_frames, cc_B, is_rev, is_out, n_ip, refined, fine_ok, ...
    'VariableNames', {'x_mm', 'z_mm', 'r_mm', 'Vx_mms', 'Vy_mms', ...
    'Vz_mms', 'VxTrue_mms', 'VyTrue_mms', 'VzTrue_mms', 'ccA', ...
    'lag_peak_frames', 'ccB', 'reverse', 'outlier', 'n_lags_inplane', ...
    'refined', 'fine'});
end


function ph = resolve_phantom(arg, d)
%RESOLVE_PHANTOM Pick the ground-truth source for this run.
% Empty means the dataset is a pure simulation with no CFD solution, so the
% analytic straight-pipe phantom is built from the slice file itself. A raw
% CFD struct is wrapped for backward compatibility with older call sites.
if isempty(arg)
    ph = src.phantom_parabolic(d);
    return;
end
if ~isstruct(arg) || ~isscalar(arg)
    error('src:estimate_slice:badPhantom', ...
        ['The second argument must be [], a phantom from ' ...
         'src.phantom_cfd / src.phantom_parabolic, or a raw CFD struct.']);
end
if isfield(arg, 'velocity') && isa(arg.velocity, 'function_handle')
    ph = arg;
    return;
end
if isfield(arg, 'grid_Uy')
    ph = src.phantom_cfd(arg);
    return;
end
error('src:estimate_slice:badPhantom', ...
    ['Unrecognised second argument: a struct with neither a velocity ' ...
     'handle (phantom) nor grid_Uy (CFD solution).']);
end
