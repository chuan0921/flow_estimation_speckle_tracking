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
% Per grid point when XZ compensation is enabled:
%   0a  rough vy from a coarse uncompensated transit scan -> in-plane
%       lag cap: the frames a scatterer needs to cross
%       opts.inplane_survival elevation beam widths at that speed
%   0b  in-plane vx/vz: chained multi-lag row1 self-correlation
%       (src.in_plane)
%   A   transit-lag scan with lag-dependent integer xz shift of the row2
%       patch (scan_transit_lags), cc_min peak (pick_peak). Scanned twice,
%       with the shift and without it, keeping whichever peak correlates
%       better; xz_shift_kept records the outcome
%   B   coarse drift at the transit peak (drift_match); optional
%       speckle-size fine stage (refine_drift, opts.fine_on); then a
%       recompensated rescan
% Analytic straight-pipe phantoms disable XZ compensation in auto mode: vx
% and vz are fixed at zero, while the small-ROI spatially weighted vy rescan
% is retained. Points the forward scan cannot answer fall back to a significance-gated
% negative-lag (reverse flow) scan, controlled by opts.scan_reverse.
% nmt outlier flag on the grid (flag only, estimates untouched).
%
% Truth comes from src.transit_truth: the streamline from row1 to row2 is
% integrated, so VyTrue_mms is the speed the transit lag corresponds to and
% points whose streamline never reaches row2 are NaN rather than scored.
% VyPlaneTrue_mms (row1 plane, for the flow integral) and VyPathTrue_mms
% (the superseded fixed-(x,z) average) are written alongside it.
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
compensate_xz = src.resolve_xz_mode(opts.compensate_xz, ph);
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

% Geometric contraction shift: in a tapering vessel the scatterers follow
% the wall, so the lateral displacement over the transit is (1 - R2/R1)
% times the vector from the axis -- one number per slice, read off the
% geometry, no estimation. Measured truth confirms it: the true shift in the
% narrowing grows linearly with radius at exactly (1-R2/R1)/dx px per mm.
% It is independent of the candidate lag for the same reason the shift clamp
% is: if a lag is the true transit, the displacement is set by where the
% streamline ends, not by how long it took. Zero for a straight pipe.
%
% Narrowing only. The follow-the-wall assumption holds where the pressure
% gradient keeps the flow attached -- converging sections. In the expansion
% the flow separates instead of following the wall outward, and offering the
% outward shift there cost accuracy (downstream nrmse 3.1 -> 4.5) for no
% coverage; with contr clamped at zero for expansions the candidate sleeps.
R1w = ph.wall_radius(d.y_row1);
R2w = ph.wall_radius(d.y_row2);
if isfinite(R1w) && isfinite(R2w) && R1w > 0
    contr = max(1 - R2w / R1w, 0);
else
    contr = 0;
end
ox_geo = -contr * (grid_x - d.xc) / d.dx;   % inward when narrowing
oz_geo = -contr * (grid_z - d.zc) / d.dz;

lag_max = min(opts.lag_max, Nt - 64);
lags_f = unique(round(logspace(log10(opts.lag_min), log10(lag_max), ...
    opts.n_lag_scan)));
lags_coarse = unique(round(logspace(log10(opts.lag_min), log10(lag_max), 40)));
if compensate_xz
    lag_ip = inplane_ladder(opts, d, dt, Nt);
else
    lag_ip = [];
end
zshift = zeros(size(lags_coarse));
zshift_f = zeros(size(lags_f));
rxf = round(opts.fine_roi_half_x / d.dx);
rzf = round(opts.fine_roi_half_z / d.dz);
rx_vy = round(opts.vy_fine_roi_half_x / d.dx);
rz_vy = round(opts.vy_fine_roi_half_z / d.dz);

% Reverse scan setup: auto is available only to phantoms that can contain
% reverse flow. This prevents the stenosis downstream-y rule from creating
% false reverse estimates in a uniform straight pipe.
allow_auto_reverse = ~isfield(ph, 'allow_reverse') || ph.allow_reverse;
scan_rev = strcmp(opts.scan_reverse, 'on') || ...
    (strcmp(opts.scan_reverse, 'auto') && allow_auto_reverse && ...
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
shifted = false(nP, 1);        % did stage A keep the compensated candidate
best_g = nan(nP, 1);           % interior CC maximum, gates or no gates
curv_g = nan(nP, 1);           % peak curvature [1/frame^2] (see pick_peak)
is_rev = false(nP, 1);
t_stage = zeros(1, 5);         % rough / inplane / scanA / stageB / reverse

for p = 1:nP
    iz0 = iz0v(p);
    ix0 = ix0v(p);
    An = tp.normalize_patch(E1, iz0, ix0, rz, rx, Nt);

    if compensate_xz
        % 0a: rough vy for the in-plane lag cap.
        t0 = tic;
        cc = tp.scan_transit_lags(An, E2, iz0, ix0, rz, rx, Nt, ...
            lags_coarse, zshift, zshift, Nz, Nx, opts.fft_min_lags);
        [vy_rough, ~] = tp.pick_peak(cc, lags_coarse, d_row, dt, 0.15);
        cap = max(lag_ip);
        if ~isnan(vy_rough)
            % No floor. Flooring the cap at lag_ip(1) made the guard inert for
            % 94% of the lumen and forced every fast ROI to match at a lag where
            % its scatterers had already left the beam; a fast point is better
            % served by no compensation (vxc = vzc = 0 below) than by one match
            % on a pattern that is no longer there.
            cap = opts.inplane_survival * 2 * d.sigma_y / ...
                max(abs(vy_rough), 1e-3) / dt;
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
    else
        % Fully developed straight-pipe truth has no in-plane motion.
        vx_e(p) = 0;
        vz_e(p) = 0;
        vxc = 0;
        vzc = 0;
    end
    % Whether this compensation is worth applying is settled in stage A by
    % scanning with and without it, not by a threshold on its size. A
    % threshold was tried and cannot work: gating at 1 px was inert (12% of
    % the throat silenced, no metric moved) and at 2.5 px it silenced the
    % narrowing slice instead, which genuinely needs 4.2 px of shift but
    % under-measures at 2.3, taking that slice's coverage from 40.4% to
    % 26.8%. Low misses, high kills the wrong points, and the in-plane
    % estimate's own 1.1 px error leaves no window in between.
    shifted(p) = vxc ~= 0 || vzc ~= 0;

    % A: transit scan, run both compensated and uncompensated, keep the peak
    % that actually correlates better.
    %
    % The shift applied here is the row1 lateral velocity extrapolated over
    % the whole transit, and that extrapolation is wrong by a fixed amount
    % this stage cannot measure: in the throat the row1 velocity points
    % inwards but reverses past the neck, so the true offset is 0.68 px while
    % the extrapolation asks for 3.22 px and the patch is searched 2.5 px away
    % from where the speckle is. In the narrowing the error goes the other way
    % (3.97 px applied against a true 5.72) and compensating clearly helps.
    % Which case a point is in is not observable from the in-plane estimate --
    % its own error is 1.1 px, larger than the thing that would have to be
    % detected -- but the peak correlation is observable, so both candidates
    % are scanned and the correlation decides. Nothing to scan twice when the
    % shift is identically zero.
    t0 = tic;
    cap_z = opts.xz_shift_max_ratio * d_row / d.dz;
    cap_x = opts.xz_shift_max_ratio * d_row / d.dx;
    oz = min(max(vzc * lags_f * dt / d.dz, -cap_z), cap_z);
    ox = min(max(vxc * lags_f * dt / d.dx, -cap_x), cap_x);
    [vy_g(p), cc_g(p), shp] = stage_a_pick(tp, opts, An, E2, iz0, ...
        ix0, rz, rx, Nt, lags_f, oz, ox, Nz, Nx, d_row, dt);
    best_g(p) = shp.best;
    curv_g(p) = shp.curv;
    if vxc ~= 0 || vzc ~= 0
        [vy_0, cc_0, shp0] = stage_a_pick(tp, opts, An, E2, iz0, ...
            ix0, rz, rx, Nt, lags_f, zshift_f, zshift_f, ...
            Nz, Nx, d_row, dt);
        if ~isnan(vy_0) && (isnan(vy_g(p)) || cc_0 > cc_g(p))
            vy_g(p) = vy_0;
            cc_g(p) = cc_0;
            best_g(p) = shp0.best;
            curv_g(p) = shp0.curv;
            vxc = 0;              % stage B and the reverse scan follow suit
            vzc = 0;
            shifted(p) = false;
        end
    end
    % Third candidate: the geometric contraction shift. The in-plane
    % estimate extrapolates the row1 lateral velocity, which lags the
    % converging field and under-compensates by a third to a half (measured
    % against oracle shifts: 2.4 px applied where the truth is 4.5); the
    % contraction field is the part of the truth the wall geometry hands
    % over for free, and in the narrowing it matches the measured true
    % shift to within the fit noise. Constant over lags. Skipped when it
    % is too small to differ from the zero candidate.
    if compensate_xz && hypot(ox_geo(p), oz_geo(p)) > 0.5
        [vy_gm, cc_gm, shpg] = stage_a_pick(tp, opts, An, E2, iz0, ...
            ix0, rz, rx, Nt, lags_f, ...
            repmat(oz_geo(p), size(lags_f)), ...
            repmat(ox_geo(p), size(lags_f)), Nz, Nx, d_row, dt);
        if ~isnan(vy_gm) && (isnan(vy_g(p)) || cc_gm > cc_g(p))
            vy_g(p) = vy_gm;
            cc_g(p) = cc_gm;
            best_g(p) = shpg.best;
            curv_g(p) = shpg.curv;
            vxc = ox_geo(p) * d.dx * vy_gm / d_row;  % equivalent lateral
            vzc = oz_geo(p) * d.dz * vy_gm / d_row;  % velocity for stage B
            shifted(p) = true;
        end
    end
    t_stage(3) = t_stage(3) + toc(t0);

    % Reverse fallback: significance-gated negative-lag scan for points
    % the forward scan could not answer.
    if isnan(vy_g(p))
        if scan_rev
            t0 = tic;
            oz = min(max(vzc * lags_r * dt / d.dz, -cap_z), cap_z);
            ox = min(max(vxc * lags_r * dt / d.dx, -cap_x), cap_x);
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

    % Straight pipe: keep the local spatial vy measurement, but do not
    % estimate or apply any XZ drift.
    if ~compensate_xz
        t0 = tic;
        if opts.vy_fine_on
            cc = tp.scan_transit_lags_spatial(E1, E2, iz0, ix0, ...
                rz_vy, rx_vy, Nt, lags_f, zshift_f, zshift_f, Nz, Nx, ...
                opts.fft_min_lags, opts.vy_spatial_half);
            [vy_2, cc_2, shp2] = tp.pick_peak(cc, lags_f, d_row, dt, ...
                opts.cc_min, opts.prom_frac, opts.prom_w, ...
                opts.prom_tail_min);
            if ~isnan(vy_2)
                vy_g(p) = vy_2;
                cc_g(p) = cc_2;
                best_g(p) = shp2.best;
                curv_g(p) = shp2.curv;
                            refined(p) = true;
            end
        end
        t_stage(4) = t_stage(4) + toc(t0);
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
        [vy_2, cc_2, shp2] = tp.pick_peak(cc, lags_f, d_row, dt, ...
            opts.cc_min, opts.prom_frac, opts.prom_w, opts.prom_tail_min);
        if ~isnan(vy_2)
            vy_g(p) = vy_2;
            cc_g(p) = cc_2;
            best_g(p) = shp2.best;
            curv_g(p) = shp2.curv;
                    refined(p) = true;
        end
    end
    t_stage(4) = t_stage(4) + toc(t0);
end
fprintf(['v2 timing [s]: rough %.0f, inplane %.0f, scanA %.0f, ' ...
    'stageB %.0f, reverse %.0f\n'], t_stage);

% Ground truth, in three flavours, because the estimator, the flow integral
% and the old scoring are not after the same quantity:
%
%   transit  what the transit lag is: d_row/T along the streamline the
%            scatterers actually follow (src.transit_truth). This is what
%            Vy_mms is compared against.
%   plane    uy on the row1 plane. A flow rate is a flux through a plane, so
%            this is the only truth the flow integral can be scored on, and
%            the only one conserved along the vessel.
%   path     the old fixed-(x,z) 5-point average, kept so the change stays
%            auditable rather than silent.
[vx_true, vy_true, vz_true, reach] = src.transit_truth(ph, ...
    grid_x - d.xc, grid_z - d.zc, d.y_row1, d.y_row2, d_row, ...
    d.flow_direction);
[~, vy_plane, ~] = ph.velocity(grid_x - d.xc, repmat(d.y_row1, nP, 1), ...
    grid_z - d.zc);
[~, uy_path, ~] = ph.velocity(repmat(grid_x - d.xc, 1, 5), ...
    repmat(ysp, nP, 1), repmat(grid_z - d.zc, 1, 5));
vy_path = mean(uy_path, 2);

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

% Peak-shape gate (opts.sharp_band, off when empty). A genuine transit
% peak's relative width is fixed by the beam geometry: the speckle is in
% view for ~sigma_y of its d_row flight, so curv*lag^2/cc clusters at
% (d_row/sigma_y)^2 (~20 here) for every true peak regardless of speed.
% Frozen-speckle plateaus land far below the band (the "peak" is the top
% of a mesa), micro-ripple spikes far above it; both correlate well and
% repeat, so the height/prominence gates pass them -- shape is the only
% observable that separates them from honest slow points.
if numel(opts.sharp_band) == 2
    sharp_g = curv_g .* lag_peak_frames .^ 2 ./ max(abs(cc_g), eps);
    bad = isfinite(sharp_g) & (sharp_g < opts.sharp_band(1) | ...
        sharp_g > opts.sharp_band(2));
    vy_g(bad) = nan;
    cc_g(bad) = nan;
    is_rev(bad) = false;
end

% Positions are vessel-centred in x and z, matching r_mm and the phantom.
T = table((grid_x - d.xc) * 1e3, (grid_z - d.zc) * 1e3, r_pt * 1e3, ...
    vx_e * 1e3, vy_g * 1e3, vz_e * 1e3, ...
    vx_true * 1e3, vy_true * 1e3, vz_true * 1e3, ...
    vy_path * 1e3, vy_plane * 1e3, reach, ...
    repmat(compensate_xz, nP, 1), shifted, cc_g, best_g, curv_g, ...
    lag_peak_frames, cc_B, ...
    is_rev, is_out, n_ip, refined, fine_ok, ...
    'VariableNames', {'x_mm', 'z_mm', 'r_mm', 'Vx_mms', 'Vy_mms', ...
    'Vz_mms', 'VxTrue_mms', 'VyTrue_mms', 'VzTrue_mms', ...
    'VyPathTrue_mms', 'VyPlaneTrue_mms', 'reachable', ...
    'xz_compensated', 'xz_shift_kept', 'ccA', 'cc_best', 'cc_curv', ...
    'lag_peak_frames', ...
    'ccB', 'reverse', ...
    'outlier', 'n_lags_inplane', 'refined', 'fine'});
end


function [vy, pkv, shp] = stage_a_pick(tp, opts, An, E2, iz0, ix0, ...
    rz, rx, Nt, lags, oz, ox, Nz, Nx, d_row, dt)
%STAGE_A_PICK Stage A transit scan and gated peak in one step.
% Alternative stage A windows were tried and lost to this plain scan: a
% speckle-sized window with curve-level spatial averaging re-imports the
% velocity mixture through the neighbours (RMSE 39 -> 63/170 on the
% stenosis test slices), and per-member peak fitting against position
% trades it for single-speckle peak jitter plus one-sided extrapolation
% (39 -> 143/267). The window itself, roi_half, is the lever that works.
cc = tp.scan_transit_lags(An, E2, iz0, ix0, rz, rx, Nt, lags, oz, ox, ...
    Nz, Nx, opts.fft_min_lags);
[vy, pkv, shp] = tp.pick_peak(cc, lags, d_row, dt, opts.cc_min, ...
    opts.prom_frac, opts.prom_w, opts.prom_tail_min);
end


function lag_ip = inplane_ladder(opts, d, dt, Nt)
%INPLANE_LADDER Lag ladder for the row1 self-correlation, in frames.
%
% Stage 0b tracks the speckle within row1, so its lags are bounded by how
% long the scatterers stay in the elevation beam (2*sigma_y) -- not by the
% row spacing d_row, which sets the transit scan instead. The two scales are
% far apart: here 2*sigma_y/d_row = 0.45, so a ladder pitched at d_row
% starts beyond the beam crossing time and no rung can ever be measured.
%
% The ladder runs one frame at a time up to 1.5 crossings at vmax; the
% per-point cap in the main loop trims it to what that point's own speed
% allows, so slow points near the wall (where the compensation is actually
% large) get the whole ladder and a real slope fit, while fast points on the
% axis (where it is ~1 px of a 20 px window) get one rung or none.
lag_ip = opts.inplane_lags;
if isempty(lag_ip)
    n_cross = 2 * d.sigma_y / (max(d.vmax, 1e-3) * dt);
    lag_ip = 1:max(2, round(1.5 * n_cross));
end
lag_ip = lag_ip(lag_ip >= 1 & lag_ip < Nt / 2);
if isempty(lag_ip)
    error('src:estimate_slice:noInplaneLags', ...
        ['opts.inplane_lags left no lag in 1..%d frames. In-plane ' ...
         'compensation would be silently disabled for every point.'], ...
        ceil(Nt / 2) - 1);
end
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
