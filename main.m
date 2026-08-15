function S = main(data_root, varargin)
%MAIN Run every speed condition, write CSVs and figures to results/.
%
%   S = main('/home/chihyu/psf_sim/data')
%   S = main(root, 'Slices', 1:4:41)          % subsample, fast pass
%   S = main(root, 'Datasets', {'sim_v020cms'})
%   S = main(root, 'Opts', struct('grid_step', 0.5e-3))
%   S = main(root, 'Figures', false)          % CSVs only
%
% data_root is either a folder of dataset folders (each holding manifest.csv
% and slices/) or a single dataset folder. Output lands in
% <project>/results/<dataset name>/, mirroring the input folder names, plus
% a _summary/ folder comparing the conditions.
%
% Per dataset:
%   points.csv         one row per grid point per slice (raw estimates)
%   profile.csv        radial profile per slice: raw mean, SG-filtered, truth
%   metrics_slice.csv  one row per slice
%   metrics.csv        one row: slice metrics averaged over the vessel
%   figures/
%
% Metrics. The lattice is identical on every slice, so radii line up and the
% analysis works at two levels, as requested:
%   per slice    RMSE_s and NRMSE_s = RMSE_s / mean(VyTrue_s)
%   per vessel   mean over slices, with the spread across slices as the std
% RMSE is reported on three bases so raw and filtered are comparable:
%   *_pt   over individual grid points
%   *_raw  over the per-radius profile means
%   *_sg   over the Savitzky-Golay filtered profile
% std_mms is the standard deviation of the point-level residual.
%
% Lag range. The stock defaults in src.default_opts target a much longer,
% faster acquisition. With 'AutoLags' (default true) the scan range is
% derived per dataset from the transit lag d_row/(vmax*dt), so a bare
% main(root) produces usable numbers. Anything set explicitly in
% 'Opts' wins over the derived value.
%
% See also SRC.ESTIMATE_SLICE, SRC.PHANTOM_PARABOLIC, SRC.PHANTOM_CFD.

root = fileparts(mfilename('fullpath'));
if isempty(which('src.default_opts'))
    % exist() reports 0 for package functions however the package resolves,
    % so the reachability test has to be which().
    addpath(root);
end

p = inputParser;
p.FunctionName = 'main';
p.addParameter('Datasets', {}, @(v) iscellstr(v) || isstring(v) || ischar(v));
p.addParameter('Slices', [], @(v) isempty(v) || isnumeric(v));
p.addParameter('Opts', struct(), @isstruct);
p.addParameter('Phantom', [], @(v) isempty(v) || isstruct(v));
p.addParameter('ResultsDir', fullfile(root, 'results'), @(v) ischar(v) || isstring(v));
p.addParameter('Figures', true, @(v) islogical(v) && isscalar(v));
p.addParameter('AutoLags', true, @(v) islogical(v) && isscalar(v));
p.addParameter('SGOrder', 3, @(v) isnumeric(v) && isscalar(v) && v >= 1);
p.addParameter('SGWindow', 11, @(v) isnumeric(v) && isscalar(v) && v >= 3);
p.addParameter('ProfileStep', 0.05, @(v) isnumeric(v) && isscalar(v) && v > 0);
p.parse(varargin{:});
o = p.Results;
o.ResultsDir = char(o.ResultsDir);

sets = find_datasets(data_root, o.Datasets);
if isempty(sets)
    error('main:noDatasets', ...
        'No dataset folder with a manifest.csv found under %s', data_root);
end
fprintf('main: %d dataset(s) under %s\n', numel(sets), data_root);

% Figures are produced off-screen; restore whatever the user had.
if o.Figures
    vis0 = get(0, 'DefaultFigureVisible');
    restore = onCleanup(@() set(0, 'DefaultFigureVisible', vis0));
    set(0, 'DefaultFigureVisible', 'off');
end

allPoints = {};
allProfile = {};
allSlice = {};
allMetrics = {};

for ds = 1:numel(sets)
    name = sets(ds).name;
    outDir = fullfile(o.ResultsDir, name);
    ensure_dir(outDir);
    fprintf('\n=== %s ===\n', name);

    M = readtable(fullfile(sets(ds).path, 'manifest.csv'));
    files = resolve_files(M, sets(ds).path);
    pick = select_slices(height(M), o.Slices);
    opts = dataset_opts(files{pick(1)}, o);

    pts = cell(1, numel(pick));
    t0 = tic;
    for k = 1:numel(pick)
        i = pick(k);
        Ti = src.estimate_slice(files{i}, o.Phantom, opts);
        Ti.dataset = repmat(string(name), height(Ti), 1);
        Ti.slice_idx = repmat(M.slice_idx(i), height(Ti), 1);
        Ti.slice_pos_mm = repmat(M.slice_pos_mm(i), height(Ti), 1);
        pts{k} = Ti;
        fprintf('  slice %3d (%+6.1f mm)  %3d/%3d valid\n', M.slice_idx(i), ...
            M.slice_pos_mm(i), nnz(isfinite(Ti.Vy_mms)), height(Ti));
    end
    P = movevars(vertcat(pts{:}), {'dataset', 'slice_idx', 'slice_pos_mm'}, ...
        'Before', 1);
    fprintf('  %d slices in %.0f s\n', numel(pick), toc(t0));

    prof = build_profile(P, o.SGOrder, o.SGWindow, ...
        slice_R(files{pick(1)}), o.ProfileStep);
    msl = slice_metrics(P, prof);
    met = vessel_metrics(msl, name, slice_vmax(files{pick(1)}));

    writetable(P, fullfile(outDir, 'points.csv'));
    writetable(prof, fullfile(outDir, 'profile.csv'));
    writetable(msl, fullfile(outDir, 'metrics_slice.csv'));
    writetable(met, fullfile(outDir, 'metrics.csv'));

    if o.Figures
        figDir = fullfile(outDir, 'figures');
        ensure_dir(figDir);
        save_fig(viz.plot_profile(prof, 'Title', name), ...
            fullfile(figDir, 'profile.png'));
        save_fig(viz.plot_error_vs_radius(prof, 'Title', name), ...
            fullfile(figDir, 'error_vs_radius.png'));
        save_fig(viz.plot_metrics_vs_slice(msl, 'Title', name), ...
            fullfile(figDir, 'metrics_vs_slice.png'));
        save_fig(viz.plot_cc_peak_map(P, 'Title', name), ...
            fullfile(figDir, 'cc_peak_map.png'));
        mid = pick(ceil(numel(pick) / 2));
        Tm = P(P.slice_idx == M.slice_idx(mid), :);
        if any(isfinite(Tm.Vy_mms))
            save_fig(viz.plot_vector_field(Tm, 'Truth', true), ...
                fullfile(figDir, 'vector_field.png'));
        end
    end

    allPoints{end + 1} = P;      %#ok<AGROW>
    allProfile{end + 1} = prof;  %#ok<AGROW>
    allSlice{end + 1} = msl;     %#ok<AGROW>
    allMetrics{end + 1} = met;   %#ok<AGROW>
end

S.points = vertcat(allPoints{:});
S.profile = vertcat(allProfile{:});
S.metrics_slice = vertcat(allSlice{:});
S.metrics = vertcat(allMetrics{:});

sumDir = fullfile(o.ResultsDir, '_summary');
ensure_dir(sumDir);
writetable(S.metrics, fullfile(sumDir, 'metrics_all.csv'));
writetable(S.metrics_slice, fullfile(sumDir, 'metrics_slice_all.csv'));

if o.Figures
    figDir = fullfile(sumDir, 'figures');
    ensure_dir(figDir);
    save_fig(viz.plot_nrmse_vs_speed(S.metrics_slice), ...
        fullfile(figDir, 'nrmse_vs_speed.png'));
    save_fig(viz.plot_valid_fraction(S.metrics_slice), ...
        fullfile(figDir, 'valid_fraction.png'));
    save_fig(viz.plot_est_vs_true(S.points), ...
        fullfile(figDir, 'est_vs_true.png'));
    save_fig(viz.plot_bland_altman(S.points), ...
        fullfile(figDir, 'bland_altman.png'));
end

fprintf('\nmain: results in %s\n', o.ResultsDir);
disp(S.metrics);
end


% ---------------------------------------------------------------- discovery

function sets = find_datasets(data_root, want)
%FIND_DATASETS A dataset is any folder holding a manifest.csv.
data_root = char(data_root);
sets = struct('name', {}, 'path', {});
if isfile(fullfile(data_root, 'manifest.csv'))
    [~, nm] = fileparts(strip_sep(data_root));
    sets(1) = struct('name', nm, 'path', data_root);
    return;
end
d = dir(data_root);
d = d([d.isdir] & ~startsWith({d.name}, '.'));
for k = 1:numel(d)
    pth = fullfile(data_root, d(k).name);
    if isfile(fullfile(pth, 'manifest.csv'))
        sets(end + 1) = struct('name', d(k).name, 'path', pth); %#ok<AGROW>
    end
end
if ~isempty(want)
    want = string(want);
    sets = sets(ismember(string({sets.name}), want));
end
end


function s = strip_sep(s)
while ~isempty(s) && (s(end) == filesep)
    s(end) = [];
end
end


function files = resolve_files(M, dsPath)
%RESOLVE_FILES Manifest paths are absolute from generation time; fall back to
% the slices/ folder next to the manifest when the data has been moved.
files = cellstr(string(M.file));
for k = 1:numel(files)
    if ~isfile(files{k})
        [~, nm, ext] = fileparts(files{k});
        files{k} = fullfile(dsPath, 'slices', [nm ext]);
    end
end
missing = ~cellfun(@isfile, files);
if any(missing)
    error('main:missingSlices', '%d slice file(s) not found, e.g. %s', ...
        nnz(missing), files{find(missing, 1)});
end
end


function pick = select_slices(n, want)
if isempty(want)
    pick = 1:n;
    return;
end
pick = want(want >= 1 & want <= n);
if isempty(pick)
    error('main:noSlices', 'Slices selection left nothing in 1..%d.', n);
end
end


function opts = dataset_opts(sample_file, o)
%DATASET_OPTS Derive a lag range from this dataset's own transit lag.
opts = o.Opts;
if ~o.AutoLags
    return;
end
s = load(sample_file, 'd_row', 'dt', 'vmax', 'Nt');
if ~all(isfield(s, {'d_row', 'dt', 'vmax', 'Nt'}))
    return;
end
L = s.d_row / (s.vmax * s.dt);        % transit lag at peak velocity [frames]
lo = max(2, floor(0.5 * L));
hi = min(s.Nt - 64, max(ceil(12 * L), lo + 20));
if ~isfield(opts, 'lag_min') || isempty(opts.lag_min)
    opts.lag_min = lo;
end
if ~isfield(opts, 'lag_max') || isempty(opts.lag_max)
    opts.lag_max = hi;
end
if ~isfield(opts, 'inplane_lags') || isempty(opts.inplane_lags)
    opts.inplane_lags = unique(max(2, round(L * [0.5 1 2 4])));
end
fprintf('  transit lag at vmax = %.1f frames -> lag scan %d..%d\n', ...
    L, opts.lag_min, opts.lag_max);
end


function v = slice_vmax(f)
s = load(f, 'vmax');
if isfield(s, 'vmax')
    v = s.vmax * 1e3;
else
    v = nan;
end
end


function r = slice_R(f)
%SLICE_R Vessel radius [mm]; the no-slip anchor sits here.
s = load(f, 'R');
if isfield(s, 'R')
    r = s.R * 1e3;
else
    r = nan;
end
end


% ------------------------------------------------------------------ analysis

function prof = build_profile(P, order, win, R_mm, step_mm)
%BUILD_PROFILE Radial velocity profile per slice, on a uniform radial
% lattice 0:step_mm:R_mm.
%
% Points sharing a radius are averaged first. Those per-radius means land on
% the lattice that carries the final answer, in this fixed order:
%
%   1 interpolate  gaps between measured radii (linear)
%   2 extrapolate  outward from the outermost measured radius (linear)
%   3 no-slip      v = 0 at r = R, which is what the extrapolation aims at
%   4 Savitzky-Golay along increasing radius
%   5 restore_wall reinstates no-slip, which step 4 does not preserve
%                  ->  vy_sg, the final result
%
% The lattice must be uniform before step 4: sgolayfilt assumes equal
% sample spacing, and the raw radii (hypot of a square grid) are anything
% but -- 37 distinct radii spaced 0.014..0.25 mm apart on the default grid.
%
% vy_raw keeps its old meaning and stays NaN away from a measured radius,
% so the scatter in viz.plot_profile and the raw metric are unaffected.
% vy_true is carried across the lattice the same way and pinned to 0 at the
% wall; beyond the outermost measured radius it is inferred, not measured,
% so slice_metrics scores vy_sg only where vy_raw is finite.
slices = unique(P.slice_idx, 'stable');
rl = (0 : step_mm : R_mm).';
if rl(end) < R_mm - 1e-9
    rl(end + 1) = R_mm;             % always land the no-slip anchor exactly
end
out = cell(1, numel(slices));
for k = 1:numel(slices)
    Q = P(P.slice_idx == slices(k), :);
    [r, ~, g] = unique(round(Q.r_mm, 6));
    n = accumarray(g, 1);
    vy = accumarray(g, Q.Vy_mms, [], @meannan);
    sd = accumarray(g, Q.Vy_mms, [], @stdnan);
    tr = accumarray(g, Q.VyTrue_mms, [], @meannan);
    nv = accumarray(g, double(isfinite(Q.Vy_mms)), [], @sum);

    [vy_l, src_l, r_out] = fill_radial(rl, r, vy, R_mm);
    tr_l = fill_radial(rl, r, tr, R_mm);
    [n_l, nv_l, sd_l, raw_l] = bin_to_lattice(rl, r, n, nv, sd, vy);

    % A measured radius rarely lands exactly on the lattice, so let the
    % binned raw value decide what counts as measured.
    src_l(isfinite(raw_l) & src_l ~= "noslip") = "meas";

    sg = restore_wall(rl, sg_filter(vy_l, order, win), r_out, R_mm);

    out{k} = table( ...
        repmat(Q.dataset(1), numel(rl), 1), ...
        repmat(slices(k), numel(rl), 1), ...
        repmat(Q.slice_pos_mm(1), numel(rl), 1), ...
        rl, n_l, nv_l, raw_l, sd_l, vy_l, sg, tr_l, src_l, ...
        'VariableNames', {'dataset', 'slice_idx', 'slice_pos_mm', 'r_mm', ...
        'n_points', 'n_valid', 'vy_raw', 'vy_std', 'vy_fill', 'vy_sg', ...
        'vy_true', 'src'});
end
prof = vertcat(out{:});
end


function [v, src_l, r_out] = fill_radial(rl, rm, vm, R_mm)
%FILL_RADIAL Measured per-radius means onto the lattice: interpolate inside,
% extrapolate outside, pin v = 0 at the wall. Radii whose points were all
% NaN drop out here and are filled as interior gaps.
% src_l is the per-lattice-radius provenance flag; it is not called src
% because that would shadow the +src package inside this function.
v = nan(size(rl));
src_l = repmat("none", numel(rl), 1);
ok = isfinite(vm) & isfinite(rm);
r_out = nan;
if ~any(ok)
    src_l(rl >= R_mm - 1e-9) = "noslip";
    v(rl >= R_mm - 1e-9) = 0;
    return;
end
rm = rm(ok);
vm = vm(ok);
r_out = max(rm);

in = rl <= r_out + 1e-9;
if sum(ok) == 1
    v(in) = vm;                     % nothing to interpolate between
else
    v(in) = interp1(rm, vm, min(rl(in), r_out), 'linear');
end
src_l(in) = "interp";
src_l(ismembertol(rl, rm, 1e-9, 'DataScale', 1)) = "meas";

% Outward: straight line from the outermost measured value to (R, 0).
out = ~in;
if r_out < R_mm - 1e-9
    v(out) = vm(end) * (R_mm - rl(out)) / (R_mm - r_out);
else
    v(out) = 0;
end
src_l(out) = "extrap";

wall = rl >= R_mm - 1e-9;
v(wall) = 0;
src_l(wall) = "noslip";
end


function y = restore_wall(rl, y, r_out, R_mm)
%RESTORE_WALL Put no-slip back after filtering. sgolayfilt does not honour
% boundary values -- run over the ramp it drags the wall sample off zero and
% can push it negative. Outside the measured range there is no noise to
% filter anyway, so the tail is rebuilt as the straight line from the
% filtered outermost measured value down to (R, 0).
if ~isfinite(r_out) || r_out >= R_mm - 1e-9
    y(rl >= R_mm - 1e-9) = 0;
    return;
end
[~, i_out] = min(abs(rl - r_out));
tail = rl > rl(i_out);
y(tail) = y(i_out) * (R_mm - rl(tail)) / (R_mm - rl(i_out));
y(rl >= R_mm - 1e-9) = 0;
end


function [n_l, nv_l, sd_l, raw_l] = bin_to_lattice(rl, rm, n, nv, sd, vy)
%BIN_TO_LATTICE Per-radius counts and raw means to their nearest lattice
% cell. raw_l stays NaN where nothing was measured, so vy_raw keeps meaning
% "a measurement happened here".
n_l = zeros(numel(rl), 1);
nv_l = zeros(numel(rl), 1);
sd_l = nan(numel(rl), 1);
raw_l = nan(numel(rl), 1);
if isempty(rm)
    return;
end
[~, j] = min(abs(rl(:).' - rm(:)), [], 2);       % nearest lattice cell
for i = 1:numel(rm)
    n_l(j(i)) = n_l(j(i)) + n(i);
    nv_l(j(i)) = nv_l(j(i)) + nv(i);
    if isfinite(vy(i))
        if isnan(raw_l(j(i)))
            raw_l(j(i)) = vy(i);
            sd_l(j(i)) = sd(i);
        else
            raw_l(j(i)) = mean([raw_l(j(i)), vy(i)]);
        end
    end
end
end


function y = sg_filter(x, order, win)
%SG_FILTER Savitzky-Golay along a gap-free uniform lattice. The window and
% order are clipped to what the sample count supports. x arrives filled, so
% there is nothing to skip -- the old version compacted around NaNs, which
% silently treated non-adjacent radii as neighbours.
y = x;
fin = isfinite(x);
if ~all(fin)
    y(:) = nan;
    return;
end
n = numel(x);
if n < 3
    return;
end
w = min(win, n);
if mod(w, 2) == 0
    w = w - 1;
end
k = min(order, w - 1);
if w < 3 || k < 1
    return;
end
y = sgolayfilt(x, k, w);
end


function msl = slice_metrics(P, prof)
%SLICE_METRICS One row per slice: point-level, profile-level and filtered.
slices = unique(P.slice_idx, 'stable');
name = strings(numel(slices), 1);
pos = nan(numel(slices), 1);
nTot = nan(numel(slices), 1);
nVal = nan(numel(slices), 1);
frac = nan(numel(slices), 1);
mtru = nan(numel(slices), 1);
bias = nan(numel(slices), 1);
sdev = nan(numel(slices), 1);
[rmPt, nrPt] = deal(nan(numel(slices), 1));
[rmRa, nrRa] = deal(nan(numel(slices), 1));
[rmSg, nrSg] = deal(nan(numel(slices), 1));

for k = 1:numel(slices)
    Q = P(P.slice_idx == slices(k), :);
    R = prof(prof.slice_idx == slices(k), :);
    name(k) = Q.dataset(1);
    pos(k) = Q.slice_pos_mm(1);
    nTot(k) = height(Q);
    m = isfinite(Q.Vy_mms) & isfinite(Q.VyTrue_mms);
    nVal(k) = nnz(m);
    frac(k) = nVal(k) / max(nTot(k), 1);
    if ~any(m)
        continue;
    end
    e = Q.Vy_mms(m) - Q.VyTrue_mms(m);
    mtru(k) = mean(Q.VyTrue_mms(m));
    bias(k) = mean(e);
    sdev(k) = std(e);
    [rmPt(k), nrPt(k)] = rmse_pair(e, mtru(k));

    mr = isfinite(R.vy_raw) & isfinite(R.vy_true);
    if any(mr)
        [rmRa(k), nrRa(k)] = rmse_pair(R.vy_raw(mr) - R.vy_true(mr), ...
            mean(R.vy_true(mr)));
    end
    % vy_sg is defined across the whole lattice now, including the
    % interpolated/extrapolated stretch. Score it only where something was
    % actually measured, so it stays comparable with rmse_raw.
    ms = isfinite(R.vy_sg) & isfinite(R.vy_true) & isfinite(R.vy_raw);
    if any(ms)
        [rmSg(k), nrSg(k)] = rmse_pair(R.vy_sg(ms) - R.vy_true(ms), ...
            mean(R.vy_true(ms)));
    end
end

msl = table(name, slices, pos, nTot, nVal, frac, mtru, bias, sdev, ...
    rmPt, nrPt, rmRa, nrRa, rmSg, nrSg, 'VariableNames', ...
    {'dataset', 'slice_idx', 'slice_pos_mm', 'n_points', 'n_valid', ...
    'valid_frac', 'mean_true_mms', 'bias_mms', 'std_mms', ...
    'rmse_pt_mms', 'nrmse_pt', 'rmse_raw_mms', 'nrmse_raw', ...
    'rmse_sg_mms', 'nrmse_sg'});
end


function [r, nr] = rmse_pair(err, ref)
r = sqrt(mean(err.^2));
nr = r / ref;
end


function met = vessel_metrics(msl, name, vmax_mms)
%VESSEL_METRICS Average the per-slice metrics over the vessel; the spread
% across slices becomes the error bar in the summary figures.
cols = {'valid_frac', 'bias_mms', 'std_mms', 'rmse_pt_mms', 'nrmse_pt', ...
    'rmse_raw_mms', 'nrmse_raw', 'rmse_sg_mms', 'nrmse_sg'};
met = table(string(name), vmax_mms, height(msl), sum(msl.n_valid), ...
    'VariableNames', {'dataset', 'vmax_mms', 'n_slices', 'n_valid_total'});
for k = 1:numel(cols)
    v = msl.(cols{k});
    met.([cols{k} '_mean']) = mean(v, 'omitnan');
    met.([cols{k} '_std']) = std(v, 'omitnan');
end
end


function m = meannan(v)
m = mean(v, 'omitnan');
end


function s = stdnan(v)
s = std(v, 'omitnan');
end


% -------------------------------------------------------------------- output

function ensure_dir(p)
if ~isfolder(p)
    mkdir(p);
end
end


function save_fig(fig, file)
if isempty(fig) || ~isgraphics(fig)
    return;
end
exportgraphics(fig, file, 'Resolution', 150);
close(fig);
end
