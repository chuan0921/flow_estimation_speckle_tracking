function S = main(data_root, varargin)
%MAIN Run every speed condition, write CSVs and figures to results/.
%
%   S = main('/home/chihyu/psf_sim/data')
%   S = main(root, 'Slices', 1:4:41)          % subsample, fast pass
%   S = main(root, 'Datasets', {'sim_v020cms'})
%   S = main(root, 'Seeds', [100 200 300])
%   S = main(root, 'Opts', struct('grid_step', 0.5e-3))
%   S = main(root, 'Figures', true)           % also draw the MATLAB figures
%
% Figures are off by default: main() writes CSVs, and make_figures.py renders
% every figure from them.
%
% data_root is either a folder of dataset folders (each holding manifest.csv
% and slices/) or a single dataset folder. Output lands in
% <project>/results/<dataset name>/, mirroring the input folder names, plus
% a _summary/ folder comparing the conditions.
%
% Per dataset:
%   points.csv              one row per grid point per slice and seed
%   field.csv               reconstructed full lumen plus no-slip wall
%   profile.csv             centre-diameter BFVP per slice and seed
%   metrics_seed_slice.csv  one row per independently estimated seed/slice
%   metrics_slice.csv       one row per slice, aggregated across seeds
%   metrics.csv        one row: slice metrics averaged over the vessel
%   figures/
%
% Seed files are discovered beside each manifest file. The unsuffixed file
% and any slice_XXX_sNNN.mat files are independent realisations. Each is
% estimated and scored separately; only then are metrics averaged across
% seeds. A single available seed therefore follows exactly the same path.
%
% Metrics compare the raw and 2-D Savitzky-Golay estimates point by point on
% the same valid ROI mask. NRMSE = RMSE / mean(VyTrue) for that seed/slice.
%
% Lag range. The stock defaults in src.default_opts target a much longer,
% faster acquisition. With 'AutoLags' (default true) the scan range is
% derived per dataset from the transit lag d_row/(vmax*dt), so a bare
% main(root) produces usable numbers. Anything set explicitly in
% 'Opts' wins over the derived value.
% opts.compensate_xz defaults to 'auto': off for analytic straight pipe and
% on for CFD. Set it to logical true/false for an explicit override.
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
p.addParameter('Seeds', 'auto', @valid_seeds);
p.addParameter('Opts', struct(), @isstruct);
p.addParameter('Phantom', [], @(v) isempty(v) || isstruct(v));
p.addParameter('ResultsDir', fullfile(root, 'results'), @(v) ischar(v) || isstring(v));
% Off by default: figures are rendered from the CSVs by make_figures.py, so
% a run costs nothing in plotting and a restyle costs nothing in compute.
p.addParameter('Figures', false, @(v) islogical(v) && isscalar(v));
p.addParameter('AutoLags', true, @(v) islogical(v) && isscalar(v));
p.addParameter('SGOrder', 2, @(v) isnumeric(v) && isscalar(v) && v >= 1);
p.addParameter('SGWindow', 5, @(v) isnumeric(v) && isscalar(v) && v >= 3);
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
allField = {};
allProfile = {};
allSeedSlice = {};
allSlice = {};
allMetrics = {};
allRegions = {};

for ds = 1:numel(sets)
    name = sets(ds).name;
    outDir = fullfile(o.ResultsDir, name);
    ensure_dir(outDir);
    fprintf('\n=== %s ===\n', name);

    M = read_manifest(fullfile(sets(ds).path, 'manifest.csv'));
    files = resolve_files(M, sets(ds).path);
    pick = select_slices(height(M), o.Slices);
    jobs = discover_seed_jobs(M, files, pick, o.Seeds);
    if isempty(jobs)
        error('main:noSeedFiles', ...
            'No files matched the requested seeds for dataset %s.', name);
    end
    opts = dataset_opts(char(jobs.file(1)), o);
    phantom = dataset_phantom(name, sets(ds).path, char(jobs.file(1)), o.Phantom);
    opts = src.default_opts(opts);
    opts.compensate_xz = src.resolve_xz_mode(opts.compensate_xz, phantom);
    if opts.compensate_xz
        xzLabel = 'on';
    else
        xzLabel = 'off';
    end
    fprintf('  xz compensation: %s\n', xzLabel);
    regions = build_vessel_regions(M, char(jobs.file(1)), phantom, opts);

    pts = cell(1, height(jobs));
    t0 = tic;
    snapIdx = [];
    snapPos = [];
    snaps = struct('E1', {}, 'E2', {}, 'x', {}, 'z', {}, 'xc', {}, ...
        'zc', {}, 'radius', {});
    for k = 1:height(jobs)
        if ~ismember(jobs.slice_idx(k), snapIdx)
            [Ti, sn] = src.estimate_slice(char(jobs.file(k)), phantom, opts);
            snapIdx(end+1) = jobs.slice_idx(k); %#ok<AGROW>
            snapPos(end+1) = jobs.slice_pos_mm(k); %#ok<AGROW>
            snaps(end+1) = sn; %#ok<AGROW>
        else
            Ti = src.estimate_slice(char(jobs.file(k)), phantom, opts);
        end
        Ti = apply_lr_table(Ti, opts, char(jobs.file(k)));
        ir = find(regions.slice_idx == jobs.slice_idx(k), 1);
        Ti.dataset = repmat(string(name), height(Ti), 1);
        Ti.slice_idx = repmat(jobs.slice_idx(k), height(Ti), 1);
        Ti.slice_pos_mm = repmat(jobs.slice_pos_mm(k), height(Ti), 1);
        Ti.seed_base = repmat(jobs.seed_base(k), height(Ti), 1);
        Ti.slice_seed = repmat(jobs.slice_seed(k), height(Ti), 1);
        Ti.lumen_radius_mm = repmat(regions.lumen_radius_mm(ir), height(Ti), 1);
        Ti.radius_ratio = repmat(regions.radius_ratio(ir), height(Ti), 1);
        Ti.truth_inplane_ratio = repmat(regions.truth_inplane_ratio(ir), height(Ti), 1);
        Ti.truth_asymmetry = repmat(regions.truth_asymmetry(ir), height(Ti), 1);
        Ti.truth_reverse_frac = repmat(regions.truth_reverse_frac(ir), height(Ti), 1);
        Ti.truth_disturbance = repmat(regions.truth_disturbance(ir), height(Ti), 1);
        Ti.vessel_region = repmat(regions.vessel_region(ir), height(Ti), 1);
        pts{k} = Ti;
        fprintf('  slice %3d seed %4g (%+6.1f mm)  %3d/%3d valid\n', ...
            jobs.slice_idx(k), jobs.seed_base(k), jobs.slice_pos_mm(k), ...
            nnz(isfinite(Ti.Vy_mms)), height(Ti));
    end
    P = movevars(vertcat(pts{:}), ...
        {'dataset', 'slice_idx', 'slice_pos_mm', 'seed_base', 'slice_seed', ...
        'vessel_region', 'lumen_radius_mm', 'radius_ratio'}, ...
        'Before', 1);
    [P, F] = reconstruct_lumen(P, o.SGOrder, o.SGWindow, ...
        phantom, char(jobs.file(1)), opts);
    fprintf('  %d slice/seed run(s), %d slice(s), in %.0f s\n', ...
        height(jobs), numel(unique(jobs.slice_idx)), toc(t0));

    % Tissue localization: the static background is the position
    % reference (direction + probe step per scan step). Report-only by
    % design -- the truth comparison downstream keeps using the manifest
    % positions, so localization error never leaks into the velocity
    % error budget; the 3-D assembly in viz is the consumer. Skipped
    % automatically when there is nothing outside the lumen to lock
    % onto (the clean simulations).
    if numel(snaps) >= 5
        [ord, ~] = sort_snaps(snapPos);
        [hasTissue, eRatio] = tissue_present(snaps(ord(1)), opts);
        if hasTissue
            dl = load(char(jobs.file(1)), 'd_row', 'y_row1', 'y_row2');
            Lloc = src.localize_scan(snaps(ord), dl.d_row, opts);
            Lloc.slice_idx = snapIdx(ord)';
            Lloc.slice_pos_mm = snapPos(ord)';
            writetable(Lloc, fullfile(outDir, 'localization.csv'));
            ok = isfinite(Lloc.step_mm);
            trueStep = mean(diff(snapPos(ord)));
            trueDir = sign(trueStep) * sign(dl.y_row2 - dl.y_row1);
            nDirOk = nnz(Lloc.direction(ok) == trueDir);
            fprintf(['  localization: dir %+d (%d/%d agree), step ' ...
                '%.4f +/- %.4f mm (manifest %.4f)\n'], trueDir, ...
                nDirOk, nnz(ok), mean(Lloc.step_mm(ok)), ...
                std(Lloc.step_mm(ok)), abs(trueStep));
        else
            fprintf(['  localization: skipped, tissue/lumen energy ' ...
                'ratio %.3f\n'], eRatio);
        end
    end

    prof = build_profile(F);
    mseed = seed_slice_metrics(P, F);
    msl = aggregate_seed_metrics(mseed);
    met = vessel_metrics(msl, name, slice_vmax(char(jobs.file(1))));

    writetable(P, fullfile(outDir, 'points.csv'));
    writetable(F, fullfile(outDir, 'field.csv'));
    writetable(prof, fullfile(outDir, 'profile.csv'));
    writetable(regions, fullfile(outDir, 'vessel_regions.csv'));
    writetable(mseed, fullfile(outDir, 'metrics_seed_slice.csv'));
    writetable(msl, fullfile(outDir, 'metrics_slice.csv'));
    writetable(met, fullfile(outDir, 'metrics.csv'));

    if o.Figures
        figDir = fullfile(outDir, 'figures');
        ensure_dir(figDir);
        Pfig = representative_seed_points(P);
        Ffig = representative_seed_points(F);
        bfvpDir = fullfile(figDir, 'bfvp');
        ensure_dir(bfvpDir);
        slicesFig = unique(Ffig.slice_idx, 'stable');
        for is = 1:numel(slicesFig)
            profFig = build_profile(Ffig(Ffig.slice_idx == slicesFig(is), :));
            save_fig(viz.plot_profile(profFig, 'Title', name), ...
                fullfile(bfvpDir, sprintf('slice_%03d.png', slicesFig(is))));
        end
        save_fig(viz.plot_metrics_vs_slice(msl, 'Title', name, ...
            'Regions', regions), fullfile(figDir, 'vessel_summary.png'));
        save_fig(viz.plot_cc_peak_map(Pfig, 'Title', name), ...
            fullfile(figDir, 'cc_peak_map.png'));
        mid = pick(ceil(numel(pick) / 2));
        Tm = Pfig(Pfig.slice_idx == M.slice_idx(mid), :);
        if any(isfinite(Tm.Vy_mms))
            save_fig(viz.plot_vector_field(Tm, 'Truth', true), ...
                fullfile(figDir, 'vector_field.png'));
        end
    end

    allPoints{end + 1} = P;      %#ok<AGROW>
    allField{end + 1} = F;       %#ok<AGROW>
    allProfile{end + 1} = prof;  %#ok<AGROW>
    allSeedSlice{end + 1} = mseed; %#ok<AGROW>
    allSlice{end + 1} = msl;     %#ok<AGROW>
    allMetrics{end + 1} = met;   %#ok<AGROW>
    allRegions{end + 1} = add_dataset(regions, name); %#ok<AGROW>
end

S.points = vertcat(allPoints{:});
S.field = vertcat(allField{:});
S.profile = vertcat(allProfile{:});
S.metrics_seed_slice = vertcat(allSeedSlice{:});
S.metrics_slice = vertcat(allSlice{:});
S.metrics = vertcat(allMetrics{:});
S.regions = vertcat(allRegions{:});

sumDir = fullfile(o.ResultsDir, '_summary');
ensure_dir(sumDir);
writetable(S.metrics, fullfile(sumDir, 'metrics_all.csv'));
writetable(S.metrics_seed_slice, fullfile(sumDir, 'metrics_seed_slice_all.csv'));
writetable(S.metrics_slice, fullfile(sumDir, 'metrics_slice_all.csv'));
writetable(S.regions, fullfile(sumDir, 'vessel_regions_all.csv'));

if o.Figures
    figDir = fullfile(sumDir, 'figures');
    ensure_dir(figDir);
    names = unique(S.metrics_slice.dataset, 'stable');
    for k = 1:numel(names)
        Q = S.metrics_slice(S.metrics_slice.dataset == names(k), :);
        R = S.regions(S.regions.dataset == names(k), :);
        save_fig(viz.plot_metrics_vs_slice(Q, 'Title', names(k), ...
            'Regions', R), fullfile(figDir, ...
            sprintf('%s_vessel_summary.png', char(names(k)))));
    end
    save_fig(viz.plot_valid_fraction(S.metrics_slice), ...
        fullfile(figDir, 'valid_fraction.png'));
    Pfig = representative_seed_points(S.points);
    save_fig(viz.plot_est_vs_true(Pfig), ...
        fullfile(figDir, 'est_vs_true.png'));
    save_fig(viz.plot_bland_altman(Pfig, ...
        'EstimateColumn', 'Vy_raw_mms', 'Title', 'Raw'), ...
        fullfile(figDir, 'bland_altman_raw.png'));
    save_fig(viz.plot_bland_altman(Pfig, ...
        'EstimateColumn', 'Vy_sg_mms', 'Title', 'After SG (2-D)'), ...
        fullfile(figDir, 'bland_altman_sg.png'));
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


function M = read_manifest(path)
%READ_MANIFEST Accept the comma- and tab-delimited manifests in the datasets.
fid = fopen(path, 'r');
if fid < 0
    error('main:manifestOpen', 'Cannot open manifest: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
header = fgetl(fid);
if ~ischar(header)
    error('main:manifestEmpty', 'Manifest is empty: %s', path);
end

if contains(header, sprintf('	'))
    delimiter = sprintf('	');
else
    delimiter = ',';
end
M = readtable(path, 'Delimiter', delimiter, 'VariableNamingRule', 'preserve');

required = {'slice_idx', 'slice_pos_mm', 'file', 'flow_direction'};
missing = required(~ismember(required, M.Properties.VariableNames));
if ~isempty(missing)
    error('main:manifestColumns', 'Manifest %s is missing column(s): %s', ...
        path, strjoin(missing, ', '));
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


function tf = valid_seeds(v)
tf = isempty(v) || isnumeric(v) || ...
    ((ischar(v) || (isstring(v) && isscalar(v))) && strcmpi(string(v), "auto"));
if isnumeric(v)
    tf = tf && isvector(v) && all(isfinite(v)) && all(v >= 0) && ...
        all(mod(v, 1) == 0);
end
end


function jobs = discover_seed_jobs(M, files, pick, want)
%DISCOVER_SEED_JOBS Expand each manifest slice into its available seeds.
% seed_base stored inside the MAT file is authoritative. Filename suffixes
% are only a fallback for older partial generations.
file = strings(0, 1);
manifest_row = zeros(0, 1);
slice_idx = zeros(0, 1);
slice_pos_mm = zeros(0, 1);
seed_base = zeros(0, 1);
slice_seed = zeros(0, 1);

for kk = 1:numel(pick)
    i = pick(kk);
    base_file = files{i};
    [folder, stem, ext] = fileparts(base_file);
    candidates = string(base_file);
    extra = dir(fullfile(folder, [stem '_s*' ext]));
    if ~isempty(extra)
        candidates = [candidates; string(fullfile({extra.folder}, {extra.name})).']; %#ok<AGROW>
    end

    seen = zeros(0, 1);
    for j = 1:numel(candidates)
        [sb, ss] = file_seed(char(candidates(j)), M.slice_idx(i));
        if ~seed_requested(sb, want)
            continue;
        end
        if ismember(sb, seen)
            warning('main:duplicateSeed', ...
                'Ignoring duplicate slice %d seed %g file %s.', ...
                M.slice_idx(i), sb, candidates(j));
            continue;
        end
        seen(end + 1, 1) = sb; %#ok<AGROW>
        file(end + 1, 1) = candidates(j); %#ok<AGROW>
        manifest_row(end + 1, 1) = i; %#ok<AGROW>
        slice_idx(end + 1, 1) = M.slice_idx(i); %#ok<AGROW>
        slice_pos_mm(end + 1, 1) = M.slice_pos_mm(i); %#ok<AGROW>
        seed_base(end + 1, 1) = sb; %#ok<AGROW>
        slice_seed(end + 1, 1) = ss; %#ok<AGROW>
    end
end

jobs = table(file, manifest_row, slice_idx, slice_pos_mm, seed_base, slice_seed);
if ~isempty(jobs)
    jobs = sortrows(jobs, {'manifest_row', 'seed_base'});
end
end


function [seed_base, slice_seed] = file_seed(file, slice_idx)
s = load(file, 'seed_base', 'slice_seed');
if isfield(s, 'seed_base') && isscalar(s.seed_base) && isfinite(s.seed_base)
    seed_base = double(s.seed_base);
else
    [~, stem] = fileparts(file);
    tok = regexp(stem, '_s(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        seed_base = str2double(tok{1});
    elseif isfield(s, 'slice_seed') && isscalar(s.slice_seed) && ...
            isfinite(s.slice_seed)
        seed_base = double(s.slice_seed) - slice_idx;
    else
        seed_base = 0;
        warning('main:legacySeed', ...
            'No seed metadata in %s; using seed_base=0.', file);
    end
end
if isfield(s, 'slice_seed') && isscalar(s.slice_seed) && isfinite(s.slice_seed)
    slice_seed = double(s.slice_seed);
else
    slice_seed = seed_base + slice_idx;
end
end


function tf = seed_requested(seed_base, want)
if isempty(want) || ischar(want) || isstring(want)
    tf = true;
else
    tf = ismember(seed_base, want);
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


function T = apply_lr_table(T, opts, file)
%APPLY_LR_TABLE Depth-dependent inter-row distance (experimental probes).
% opts.d_row_table is [depth_mm, lr_mm]; the through-plane velocity is
% scaled by lr(z)/d_row per ROI. v = d_row/(lag*dt) is linear in d_row,
% so this post-hoc scaling equals re-estimating with the local
% separation. Depths outside the table clamp to its end values.
% Simulation datasets never pass the table, so they are untouched.
if ~isfield(opts, 'd_row_table') || isempty(opts.d_row_table)
    return;
end
tab = opts.d_row_table;
d = load(file, 'zc', 'd_row');
z_abs = d.zc * 1e3 + T.z_mm;               % ROI absolute depth [mm]
lr = interp1(tab(:, 1), tab(:, 2), z_abs, 'linear');
lr(z_abs <= tab(1, 1)) = tab(1, 2);
lr(z_abs >= tab(end, 1)) = tab(end, 2);
T.lr_mm = lr;
T.Vy_mms = T.Vy_mms .* (lr / (d.d_row * 1e3));
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
% 12*L caps the scan at 87 frames here, a floor of 77 mm/s, and the
% wall-hugging ring below that floor does carry usable correlation (0.55 at
% its true lag of 160-280 frames, measured with oracle compensation) -- but
% every attempt to scan for it lost more than it gained. Opened to Nt-64,
% near-wall windows preferred the 400+ frame frozen tail of their own slow
% sliver over the genuine peak (36 confident wrong answers per 3 slices); a
% prominence gate wide enough to reject those plateaus rejected the genuine
% 160-280 frame peaks too (both live inside the gate's edge margin), and a
% 24*L compromise polluted the 0.67-0.85 ring instead. The slow sliver and
% the genuine peak are the same window's velocity mixture, so no lag-domain
% gate separates them; the floor stays.
hi = min(s.Nt - 64, max(ceil(12 * L), lo + 20));
if ~isfield(opts, 'lag_min') || isempty(opts.lag_min)
    opts.lag_min = lo;
end
if ~isfield(opts, 'lag_max') || isempty(opts.lag_max)
    opts.lag_max = hi;
end
% inplane_lags is deliberately not derived here. It belongs to the elevation
% beam rather than to the row spacing, so src.estimate_slice builds it from
% sigma_y, which this function does not load. Deriving it from L put every
% rung past the beam crossing time.
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


function ph = dataset_phantom(name, dsPath, sample_file, supplied)
%DATASET_PHANTOM Use CFD truth automatically for stenosis simulations.
if ~isempty(supplied)
    if isfield(supplied, 'wall_radius') && isfield(supplied, 'velocity')
        ph = supplied;
    else
        ph = src.phantom_cfd(supplied);
    end
    return;
end

tok = regexp(char(name), '^sim_stenosis(\d*)_?(\d+)cms', 'tokens', 'once');
if ~isempty(tok)
    degree = tok{1};
    if isempty(degree)
        degree = '50';
    end
    psf_root = fileparts(fileparts(dsPath));
    cfd_file = fullfile(psf_root, 'CFD_input', ...
        sprintf('%s_vessel_Vcenter_%scm_s', degree, tok{2}), ...
        'cfd_flow_grid.mat');
    if ~isfile(cfd_file)
        error('main:missingCFD', ...
            'Stenosis dataset %s requires CFD truth: %s', name, cfd_file);
    end
    fprintf('  truth: %s\n', cfd_file);
    ph = src.phantom_cfd(load(cfd_file));
else
    d = load(sample_file, 'R', 'vmax', 'flow_direction');
    ph = src.phantom_parabolic(d);
end
end


function regions = build_vessel_regions(M, sample_file, ph, opts)
%BUILD_VESSEL_REGIONS Classify the full vessel from geometry and CFD truth.
d = load(sample_file, 'R', 'd_row');
opts = src.default_opts(opts);
n = height(M);
radius = nan(n, 1);
ipRatio = nan(n, 1);
asymmetry = nan(n, 1);
reverseFrac = nan(n, 1);

kmax = floor(d.R / opts.grid_step);
gv = (-kmax:kmax) * opts.grid_step;
[GX, GZ] = meshgrid(gv, gv);

for k = 1:n
    y0 = M.slice_pos_mm(k) * 1e-3;
    ysp = linspace(y0 - d.d_row / 2, y0 + d.d_row / 2, 5);
    rpath = ph.wall_radius(ysp);
    radius(k) = min(rpath, [], 'omitnan') * 1e3;
    keep = hypot(GX, GZ) <= 0.9 * radius(k) * 1e-3;
    xq = GX(keep);
    zq = GZ(keep);
    Xq = repmat(xq, 1, numel(ysp));
    Zq = repmat(zq, 1, numel(ysp));
    Yq = repmat(ysp, numel(xq), 1);
    [ux, uy, uz] = ph.velocity(Xq, Yq, Zq);
    ux = mean(ux, 2);
    uy = mean(uy, 2);
    uz = mean(uz, 2);
    [ipRatio(k), asymmetry(k), reverseFrac(k)] = ...
        truth_disturbance_metrics(xq, zq, ux, uy, uz);
end

disturbance = hypot(ipRatio, asymmetry) + reverseFrac;
[region, threshold] = classify_regions(M.slice_pos_mm, radius, disturbance);
baseline = max(radius, [], 'omitnan');
ratio = radius / baseline;
regionId = region_ids(region);

regions = table(M.slice_idx, M.slice_pos_mm, radius, ratio, ipRatio, ...
    asymmetry, reverseFrac, disturbance, repmat(threshold, n, 1), ...
    regionId, region, 'VariableNames', {'slice_idx', 'slice_pos_mm', ...
    'lumen_radius_mm', 'radius_ratio', 'truth_inplane_ratio', ...
    'truth_asymmetry', 'truth_reverse_frac', 'truth_disturbance', ...
    'disturbance_threshold', 'region_id', 'vessel_region'});
end


function [ipRatio, asymmetry, reverseFrac] = ...
        truth_disturbance_metrics(x, z, ux, uy, uz)
ok = isfinite(ux) & isfinite(uy) & isfinite(uz);
if ~any(ok)
    [ipRatio, asymmetry, reverseFrac] = deal(nan);
    return;
end
x = x(ok);
z = z(ok);
ux = ux(ok);
uy = uy(ok);
uz = uz(ok);
den = sqrt(mean(uy.^2));
ipRatio = sqrt(mean(ux.^2 + uz.^2)) / max(den, eps);

[~, ~, g] = unique(round(hypot(x, z), 9));
radialMean = accumarray(g, uy, [], @mean);
asymmetry = sqrt(mean((uy - radialMean(g)).^2)) / max(den, eps);

peak = max(abs(uy));
if peak < eps
    reverseFrac = 0;
else
    reverseFrac = mean(uy < -0.02 * peak);
end
end


function [region, disturbanceThreshold] = classify_regions(y, radius, disturbance)
%CLASSIFY_REGIONS Geometry defines the stenosis; truth defines disturbance.
[ys, order] = sort(y);
rs = radius(order);
ds = disturbance(order);
n = numel(ys);
labels = repmat("uniform vessel", n, 1);

baseline = max(rs, [], 'omitnan');
minimum = min(rs, [], 'omitnan');
depth = baseline - minimum;
if ~isfinite(depth) || depth / baseline < 0.05
    disturbanceThreshold = nan;
    region = strings(n, 1);
    region(order) = labels;
    return;
end

geometryThreshold = baseline - 0.10 * depth;
throatThreshold = minimum + 0.15 * depth;
narrow = find(rs < geometryThreshold);
throat = find(rs <= throatThreshold);
[~, iMinimum] = min(rs);
if isempty(narrow)
    narrow = iMinimum;
end
if isempty(throat)
    throat = iMinimum;
end
onset = narrow(1);
throatStart = throat(1);
throatEnd = throat(end);

labels(:) = "pre-stenosis";
labels(onset:max(onset, throatStart - 1)) = "narrowing";
labels(throatStart:throatEnd) = "stenosis throat";
if throatEnd < n
    labels(throatEnd + 1:end) = "post-stenosis";
end

upstream = ds(1:max(onset - 1, 1));
upstream = upstream(isfinite(upstream));
if isempty(upstream)
    baseDisturbance = 0;
    robustSigma = 0;
else
    baseDisturbance = median(upstream);
    robustSigma = 1.4826 * median(abs(upstream - baseDisturbance));
end
disturbanceThreshold = max(0.08, baseDisturbance + 3 * robustSigma);

disturbed = isfinite(ds) & ds > disturbanceThreshold;
disturbanceStart = first_sustained(disturbed, throatEnd + 1, 2);
if ~isempty(disturbanceStart)
    labels(disturbanceStart:end) = "disturbed flow";
    recoveredGeometry = rs >= geometryThreshold;
    recoveredFlow = isfinite(ds) & ds <= disturbanceThreshold & recoveredGeometry;
    recoveryStart = first_sustained(recoveredFlow, disturbanceStart + 1, 3);
    if ~isempty(recoveryStart)
        labels(recoveryStart:end) = "recovery";
    end
end

region = strings(n, 1);
region(order) = labels;
end


function idx = first_sustained(mask, startAt, runLength)
idx = [];
startAt = max(1, startAt);
for k = startAt:max(startAt, numel(mask) - runLength + 1)
    if k + runLength - 1 <= numel(mask) && all(mask(k:k + runLength - 1))
        idx = k;
        return;
    end
end
end


function id = region_ids(region)
names = ["uniform vessel", "pre-stenosis", "narrowing", ...
    "stenosis throat", "post-stenosis", "disturbed flow", "recovery"];
id = nan(size(region));
for k = 1:numel(names)
    id(region == names(k)) = k - 1;
end
end


function T = add_dataset(T, name)
T.dataset = repmat(string(name), height(T), 1);
T = movevars(T, 'dataset', 'Before', 1);
end


% ------------------------------------------------------------------ analysis

function [P, field] = reconstruct_lumen(P, order, win, ph, sample_file, opts)
%RECONSTRUCT_LUMEN Full 2-D field with exact no-slip wall anchors.
% Raw estimates remain untouched. Missing interior and near-wall values are
% interpolated only inside the known lumen, SG is applied to that field, and
% synthetic wall samples are reset to exactly zero after filtering.
%
% The lumen is taken at the row1 plane, which is where src.geometry places
% the ROIs. Streamlines follow the taper, so a scatterer near the wall at
% row1 stays in the lumen all the way to row2 -- it simply arrives at a
% different radius.
d = load(sample_file, 'd_row', 'flow_direction');
opts = src.default_opts(opts);
step = opts.grid_step * 1e3;
P.Vy_raw_mms = P.Vy_mms;
P.Vy_fill_mms = nan(height(P), 1);
P.Vy_sg_mms = nan(height(P), 1);

G = findgroups(P.dataset, P.slice_idx, P.seed_base);
nG = max(G);
out = cell(1, nG);

% Pass 1: geometry, interpolation fill, profile fit and exam per slice.
% The blend is NOT applied here. Its three parameters (v0, n, lambda) come
% from ~100 measured points and a ~15-point exam per slice, and that
% sampling noise lands directly on the fill as slice-to-slice jitter: a
% degenerate plug fit returns n = 216 between neighbours at 8-10, and a
% lambda of 0.73 between neighbours at 0.98 costs the slice 4 NRMSE
% points. Adjacent slices are 0.5 mm apart, so the true profile cannot
% jump; the parameters are median-smoothed over the +/-2 neighbouring
% slices between the passes, which removes exactly those outliers and
% leaves smooth stretches untouched.
S = cell(1, nG);
for k = 1:nG
    st = struct();
    st.idx = find(G == k);
    Q = P(st.idx, :);
    y0 = Q.slice_pos_mm(1) * 1e-3;
    st.y0 = y0;
    radius = ph.wall_radius(y0 - d.d_row / 2) * 1e3;   % row1 plane
    if ~isfinite(radius)
        radius = Q.lumen_radius_mm(1);
    end
    st.radius = radius;

    ng = ceil(radius / step);
    gv = (-ng:ng) * step;
    [GX, GZ] = meshgrid(gv, gv);
    keep = hypot(GX, GZ) <= radius + 1e-9;
    st.xg = GX(keep);
    st.zg = GZ(keep);
    nGrid = numel(st.xg);

    nWall = max(48, ceil(2 * pi * radius / (step / 2)));
    theta = (0:nWall-1).' * (2 * pi / nWall);
    st.xWall = radius * cos(theta);
    st.zWall = radius * sin(theta);
    st.nWall = nWall;

    [st.roiPresent, roiRow] = ismember(round([st.xg st.zg], 9), ...
        round([Q.x_mm Q.z_mm], 9), 'rows');
    rawGrid = nan(nGrid, 1);
    rawGrid(st.roiPresent) = Q.Vy_raw_mms(roiRow(st.roiPresent));
    st.rawGrid = rawGrid;

    measured = isfinite(Q.Vy_raw_mms);
    % A Cartesian ROI can land exactly on the circular wall. Exclude it from
    % interpolation sources so the exact no-slip anchor wins instead of
    % scatteredInterpolant averaging two values at the same coordinate.
    sourceMeasured = measured & hypot(Q.x_mm, Q.z_mm) < radius - 1e-9;
    st.canFill = nnz(sourceMeasured) >= 3;
    st.fillGrid = nan(nGrid, 1);
    st.v0 = nan;
    st.n = nan;
    st.lam = nan;
    st.expanding = false;
    if st.canFill
        sourceX = [Q.x_mm(sourceMeasured); st.xWall];
        sourceZ = [Q.z_mm(sourceMeasured); st.zWall];
        sourceV = [Q.Vy_raw_mms(sourceMeasured); zeros(nWall, 1)];
        interpObj = scatteredInterpolant(sourceX, sourceZ, sourceV, ...
            'natural', 'none');
        fillGrid = interpObj(st.xg, st.zg);
        missing = ~isfinite(fillGrid);
        if any(missing)
            nearestObj = scatteredInterpolant(sourceX, sourceZ, sourceV, ...
                'nearest', 'nearest');
            fillGrid(missing) = nearestObj(st.xg(missing), st.zg(missing));
        end
        st.fillGrid = fillGrid;

        % The interpolation above runs straight lines to the zero anchors
        % on the wall, which cuts the shoulder off a blunt stenotic jet; a
        % robust power-law profile fitted to the measured points keeps the
        % shoulder but cannot represent an asymmetric field. Each slice
        % sits both on a held-out exam (the outer 20% of measured points)
        % and mixes them by inverse squared exam error, so a slice where
        % the model is wrong falls back to interpolation with no tuned
        % threshold. 41-slice stenosis run: full-lumen NRMSE 24.1% ->
        % 14.2%; straight vessels improve slightly (the parabola is the
        % n = 2 member of the family).
        %
        % The model assumes attached flow. Downstream of the throat the
        % jet eventually separates, the near-wall truth drops below any
        % blunt profile, and the exam cannot see it because the
        % separation lives in the unmeasured ring. The hard geometric cut
        % (expanding -> no model) proved too early -- the jet stays
        % attached well past the throat (sl22-31, full 18-32% -> 4-18%
        % when allowed) -- so expansion only disqualifies the model once
        % the measured profile itself has relaxed: the fitted exponent
        % stays blunt (n = 6.4-9.7) while attached and drops to 4.6-5.0
        % when separation appears. 5.5 splits that gap; calibrated on
        % this vessel at one seed, the constant to watch in multi-seed
        % runs.
        fd = d.flow_direction;
        rev = (isnumeric(fd) && ~isempty(fd) && fd(1) < 0) || ...
            (~isnumeric(fd) && strcmpi(strtrim(char(fd)), 'row2_to_row1'));
        if rev
            yUp = y0 + d.d_row / 2;
            yDn = y0 - d.d_row / 2;
        else
            yUp = y0 - d.d_row / 2;
            yDn = y0 + d.d_row / 2;
        end
        Rup = ph.wall_radius(yUp);
        Rdn = ph.wall_radius(yDn);
        st.expanding = isfinite(Rup) && isfinite(Rdn) && ...
            Rdn > Rup * (1 + 1e-6);

        fitOk = sourceMeasured & isfinite(Q.ccA);
        if nnz(fitOk) >= 20 && any(~isfinite(rawGrid))
            xf = Q.x_mm(fitOk);
            zf = Q.z_mm(fitOk);
            rf = hypot(xf, zf);
            vf = Q.Vy_raw_mms(fitOk);
            wf = Q.ccA(fitOk);
            qFull = fit_powerlaw(rf, vf, wf, radius, [max(vf) 3]);
            st.v0 = qFull(1);
            st.n = abs(qFull(2));
            rs = sort(rf);
            cut = interp1(linspace(0, 1, numel(rs)), rs, 0.8);
            tr = rf <= cut;
            te = ~tr;
            if nnz(te) >= 3 && nnz(tr) >= 5
                qTr = fit_powerlaw(rf(tr), vf(tr), wf(tr), radius, ...
                    [max(vf) 3]);
                predA = powerlaw_v(rf(te), qTr, radius);
                trX = [xf(tr); st.xWall];
                trZ = [zf(tr); st.zWall];
                trV = [vf(tr); zeros(nWall, 1)];
                itpB = scatteredInterpolant(trX, trZ, trV, 'linear', 'none');
                xTe = xf(te);
                zTe = zf(te);
                predB = itpB(xTe, zTe);
                nb = ~isfinite(predB);
                if any(nb)
                    itpN = scatteredInterpolant(trX, trZ, trV, ...
                        'nearest', 'nearest');
                    predB(nb) = itpN(xTe(nb), zTe(nb));
                end
                eA = sqrt(mean((predA - vf(te)) .^ 2));
                eB = sqrt(mean((predB - vf(te)) .^ 2));
                st.lam = eB^2 / (eA^2 + eB^2);
            end
        end
    end
    S{k} = st;
end
S = [S{:}];

% Median over the running +/-2 slices of the same dataset and seed.
key = strings(1, nG);
pos = zeros(1, nG);
for k = 1:nG
    Q1 = P(S(k).idx(1), :);
    key(k) = string(Q1.dataset) + "#" + string(Q1.seed_base);
    pos(k) = Q1.slice_pos_mm;
end
v0s = [S.v0];
ns = [S.n];
lams = [S.lam];
for u = unique(key)
    sel = find(key == u);
    [~, ord] = sort(pos(sel));
    sel = sel(ord);
    v0s(sel) = med_neighbours(v0s(sel));
    ns(sel) = med_neighbours(ns(sel));
    lams(sel) = med_neighbours(lams(sel));
end

% Pass 2: blend with the smoothed parameters, output smoothing, truths.
for k = 1:nG
    st = S(k);
    idx = st.idx;
    Q = P(idx, :);
    xg = st.xg;
    zg = st.zg;
    nGrid = numel(xg);
    nWall = st.nWall;
    rawGrid = st.rawGrid;
    fillGrid = st.fillGrid;
    holes = ~isfinite(rawGrid);
    separated = st.expanding && ns(k) <= 5.5;
    if st.canFill && ~separated && isfinite(ns(k)) && ...
            isfinite(lams(k)) && lams(k) > 0 && any(holes)
        vModel = powerlaw_v(hypot(xg(holes), zg(holes)), ...
            [v0s(k) ns(k)], st.radius);
        fillGrid(holes) = lams(k) * vModel + (1 - lams(k)) * fillGrid(holes);
    end
    if ~st.canFill
        sgGrid = nan(nGrid, 1);
    else
        % SG is exact for a quadratic field and destructive for a jet
        % cliff, so the smoothed exponent decides per slice instead of
        % the vessel type: straight pipes fit n = 2.0-2.2 (parabola,
        % SG is free variance reduction, raw 7.0 -> 2.7% at window 11),
        % every stenosis slice fits n = 3.9-10 (blunt profile, SG blurs
        % the cliff: 10.5 -> 10.8-11.2% measured, +2.3 points at the
        % throat). 3.0 splits the observed gap; like the separation
        % gate it is calibrated on these vessels at one seed. A slice
        % with no profile fit falls back to the tapered-vessel proxy.
        if isfinite(ns(k))
            doSG = ns(k) < 3;
        else
            doSG = ~opts.compensate_xz;
        end
        if doSG
            sgGrid = local_poly_filter([xg; st.xWall], [zg; st.zWall], ...
                [fillGrid; zeros(nWall, 1)], xg, zg, order, win, step);
        else
            sgGrid = fillGrid;
        end
    end

    % Two truths on the reconstruction grid. The field is compared against
    % the transit truth, because that is what every sample in it estimates;
    % the flow is compared against the row1 plane, because a flow rate is a
    % flux through a plane and only the plane version is conserved along the
    % vessel (465 +/- 4 ml/min over all 41 slices, against 452-473 for the
    % old fixed-column average).
    y0 = st.y0;
    y1 = y0 - d.d_row / 2;
    [~, truthTransit] = src.transit_truth(ph, xg * 1e-3, zg * 1e-3, ...
        y1, y0 + d.d_row / 2, d.d_row, d.flow_direction);
    truthGrid = truthTransit * 1e3;
    [~, truthPlane, ~] = ph.velocity(xg * 1e-3, repmat(y1, nGrid, 1), ...
        zg * 1e-3);
    truthPlane = truthPlane * 1e3;

    [matched, gridRow] = ismember(round([Q.x_mm Q.z_mm], 9), ...
        round([xg zg], 9), 'rows');
    P.Vy_fill_mms(idx(matched)) = fillGrid(gridRow(matched));
    P.Vy_sg_mms(idx(matched)) = sgGrid(gridRow(matched));

    sampleType = [repmat("grid", nGrid, 1); repmat("wall", nWall, 1)];
    roiPresentAll = [st.roiPresent; false(nWall, 1)];
    measuredAll = [isfinite(rawGrid); false(nWall, 1)];
    area = [repmat(step^2, nGrid, 1); zeros(nWall, 1)];
    rawAll = [rawGrid; nan(nWall, 1)];
    fillAll = [fillGrid; zeros(nWall, 1)];
    sgAll = [sgGrid; zeros(nWall, 1)];
    truthAll = [truthGrid; zeros(nWall, 1)];
    planeAll = [truthPlane; zeros(nWall, 1)];
    xAll = [xg; st.xWall];
    zAll = [zg; st.zWall];
    nAll = nGrid + nWall;

    out{k} = table( ...
        repmat(Q.dataset(1), nAll, 1), ...
        repmat(Q.slice_idx(1), nAll, 1), ...
        repmat(Q.slice_pos_mm(1), nAll, 1), ...
        repmat(Q.seed_base(1), nAll, 1), ...
        repmat(Q.slice_seed(1), nAll, 1), ...
        repmat(Q.vessel_region(1), nAll, 1), ...
        repmat(st.radius, nAll, 1), xAll, zAll, hypot(xAll, zAll), ...
        sampleType, roiPresentAll, measuredAll, area, rawAll, fillAll, ...
        sgAll, truthAll, planeAll, 'VariableNames', {'dataset', 'slice_idx', ...
        'slice_pos_mm', 'seed_base', 'slice_seed', 'vessel_region', ...
        'lumen_radius_mm', 'x_mm', 'z_mm', 'r_mm', 'sample_type', ...
        'roi_present', 'is_measured', 'cell_area_mm2', 'Vy_raw_mms', ...
        'Vy_fill_mms', 'Vy_sg_mms', 'VyTrue_mms', 'VyPlaneTrue_mms'});
end
field = vertcat(out{:});
P = movevars(P, {'Vy_raw_mms', 'Vy_fill_mms', 'Vy_sg_mms'}, ...
    'After', 'Vy_mms');
end


function y = local_poly_filter(x, z, v, xq, zq, order, win, step)
%LOCAL_POLY_FILTER Savitzky-Golay equivalent on a masked 2-D neighbourhood.
y = nan(size(xq));
finite = isfinite(x) & isfinite(z) & isfinite(v);
w = floor(win);
if mod(w, 2) == 0
    w = w - 1;
end
halfWidth = max(1, (w - 1) / 2) * step;
tol = 1e-6 * max(step, 1);

for i = 1:numel(xq)
    nb = finite & abs(x - xq(i)) <= halfWidth + tol & ...
        abs(z - zq(i)) <= halfWidth + tol;
    xn = (x(nb) - xq(i)) / step;
    zn = (z(nb) - zq(i)) / step;
    vn = v(nb);
    fitOrder = floor(order);
    while fitOrder >= 0
        A = polynomial_terms(xn, zn, fitOrder);
        if size(A, 1) >= size(A, 2) && rank(A) == size(A, 2)
            beta = A \ vn;
            y(i) = beta(1);
            break;
        end
        fitOrder = fitOrder - 1;
    end
end
end


function A = polynomial_terms(x, z, order)
A = ones(numel(x), 1);
for degree = 1:order
    for px = degree:-1:0
        pz = degree - px;
        A(:, end + 1) = x.^px .* z.^pz; %#ok<AGROW>
    end
end
end


function prof = build_profile(field)
%BUILD_PROFILE Velocity along the vessel-centre x-diameter (z = 0).
% The signed coordinate preserves the two sides of an asymmetric or
% disturbed profile. Quantitative metrics use the full 2-D field, not this
% diagnostic cut.
G = findgroups(field.dataset, field.slice_idx, field.seed_base);
out = cell(1, max(G));
for k = 1:max(G)
    Q = field(G == k, :);
    radius = Q.lumen_radius_mm(1);
    gridRows = Q.sample_type == "grid";
    zLevels = abs(Q.z_mm(gridRows));
    z0 = min(zLevels, [], 'omitnan');
    D = Q(gridRows & abs(abs(Q.z_mm) - z0) <= 1e-9, :);
    D = sortrows(D, 'x_mm');

    % Exact no-slip endpoints replace any coincident Cartesian boundary row.
    interior = abs(D.x_mm) < radius - 1e-9;
    D = D(interior, :);
    x = [-radius; D.x_mm; radius];
    n = [0; ones(height(D), 1); 0];
    raw = [nan; D.Vy_raw_mms; nan];
    nv = double(isfinite(raw));
    sd = nan(size(x));
    fillv = [0; D.Vy_fill_mms; 0];
    sg = [0; D.Vy_sg_mms; 0];
    truth = [0; D.VyTrue_mms; 0];
    src_l = ["no-slip wall"; repmat("centre diameter", height(D), 1); ...
        "no-slip wall"];

    out{k} = table( ...
        repmat(Q.dataset(1), numel(x), 1), ...
        repmat(Q.slice_idx(1), numel(x), 1), ...
        repmat(Q.slice_pos_mm(1), numel(x), 1), ...
        repmat(Q.seed_base(1), numel(x), 1), ...
        repmat(Q.slice_seed(1), numel(x), 1), ...
        repmat(Q.vessel_region(1), numel(x), 1), ...
        repmat(radius, numel(x), 1), x, n, nv, raw, sd, fillv, sg, truth, src_l, ...
        'VariableNames', {'dataset', 'slice_idx', 'slice_pos_mm', ...
        'seed_base', 'slice_seed', 'vessel_region', 'lumen_radius_mm', ...
        'diameter_pos_mm', 'n_points', 'n_valid', 'vy_raw', 'vy_std', ...
        'vy_fill', 'vy_sg', 'vy_true', 'src'});
end
prof = vertcat(out{:});
end


function mseed = seed_slice_metrics(P, field)
%SEED_SLICE_METRICS Score each independent seed/slice before aggregation.
G = findgroups(P.dataset, P.slice_idx, P.seed_base);
n = max(G);
name = strings(n, 1);
[slice, pos, seedBase, sliceSeed, nTot, nVal, frac, mtru] = ...
    deal(nan(n, 1));
[lumenRadius, radiusRatio, truthIp, truthAsym, truthReverse, truthDisturbance] = ...
    deal(nan(n, 1));
region = strings(n, 1);
[biasRaw, stdRaw, rmRaw, nrRaw, biasSg, stdSg, rmSg, nrSg] = ...
    deal(nan(n, 1));
[nFull, rmFillFull, nrFillFull, rmSgFull, nrSgFull, flowFill, flowSg, ...
    flowTruth, flowTransit, flowErrFill, flowErrSg] = deal(nan(n, 1));

for k = 1:n
    Q = P(G == k, :);
    name(k) = Q.dataset(1);
    slice(k) = Q.slice_idx(1);
    pos(k) = Q.slice_pos_mm(1);
    seedBase(k) = Q.seed_base(1);
    sliceSeed(k) = Q.slice_seed(1);
    lumenRadius(k) = Q.lumen_radius_mm(1);
    radiusRatio(k) = Q.radius_ratio(1);
    truthIp(k) = Q.truth_inplane_ratio(1);
    truthAsym(k) = Q.truth_asymmetry(1);
    truthReverse(k) = Q.truth_reverse_frac(1);
    truthDisturbance(k) = Q.truth_disturbance(1);
    region(k) = Q.vessel_region(1);
    nTot(k) = height(Q);
    R = field(field.dataset == Q.dataset(1) & ...
        field.slice_idx == Q.slice_idx(1) & ...
        field.seed_base == Q.seed_base(1) & field.sample_type == "grid", :);
    mf = isfinite(R.Vy_fill_mms) & isfinite(R.Vy_sg_mms) & ...
        isfinite(R.VyTrue_mms);
    nFull(k) = nnz(mf);
    if any(mf)
        truthFull = R.VyTrue_mms(mf);
        refFull = mean(truthFull);
        [rmFillFull(k), nrFillFull(k)] = rmse_pair( ...
            R.Vy_fill_mms(mf) - truthFull, refFull);
        [rmSgFull(k), nrSgFull(k)] = rmse_pair( ...
            R.Vy_sg_mms(mf) - truthFull, refFull);
    end
    % The flow integral runs on its own mask and its own truth. A cell whose
    % streamline never reaches row2 has no transit truth, but it still
    % carries flux through the row1 plane, so excluding it would understate
    % the true flow rather than the estimate. flow_transit_true_ml_min
    % integrates the transit truth over the same cells, which separates the
    % method's own bias -- the estimator reports a transit average where the
    % flux wants a plane value -- from the estimation error.
    mq = isfinite(R.Vy_fill_mms) & isfinite(R.Vy_sg_mms) & ...
        isfinite(R.VyPlaneTrue_mms);
    if any(mq)
        area = R.cell_area_mm2(mq);
        flowFill(k) = sum(R.Vy_fill_mms(mq) .* area) * 0.06;
        flowSg(k) = sum(R.Vy_sg_mms(mq) .* area) * 0.06;
        flowTruth(k) = sum(R.VyPlaneTrue_mms(mq) .* area) * 0.06;
        tt = R.VyTrue_mms(mq);
        flowTransit(k) = sum(tt(isfinite(tt)) .* area(isfinite(tt))) * 0.06;
        if abs(flowTruth(k)) >= eps
            flowErrFill(k) = 100 * (flowFill(k) - flowTruth(k)) / ...
                abs(flowTruth(k));
            flowErrSg(k) = 100 * (flowSg(k) - flowTruth(k)) / ...
                abs(flowTruth(k));
        end
    end

    m = isfinite(Q.Vy_raw_mms) & isfinite(Q.Vy_sg_mms) & ...
        isfinite(Q.VyTrue_mms);
    nVal(k) = nnz(m);
    frac(k) = nVal(k) / max(nTot(k), 1);
    if ~any(m)
        continue;
    end

    truth = Q.VyTrue_mms(m);
    eRaw = Q.Vy_raw_mms(m) - truth;
    eSg = Q.Vy_sg_mms(m) - truth;
    mtru(k) = mean(truth);
    biasRaw(k) = mean(eRaw);
    stdRaw(k) = std(eRaw);
    [rmRaw(k), nrRaw(k)] = rmse_pair(eRaw, mtru(k));
    biasSg(k) = mean(eSg);
    stdSg(k) = std(eSg);
    [rmSg(k), nrSg(k)] = rmse_pair(eSg, mtru(k));
end

mseed = table(name, slice, pos, seedBase, sliceSeed, region, lumenRadius, ...
    radiusRatio, truthIp, truthAsym, truthReverse, truthDisturbance, ...
    nTot, nVal, frac, mtru, biasRaw, stdRaw, rmRaw, nrRaw, ...
    biasSg, stdSg, rmSg, nrSg, nFull, rmFillFull, nrFillFull, rmSgFull, ...
    nrSgFull, flowFill, flowSg, flowTruth, flowTransit, flowErrFill, ...
    flowErrSg, ...
    'VariableNames', {'dataset', 'slice_idx', 'slice_pos_mm', ...
    'seed_base', 'slice_seed', 'vessel_region', 'lumen_radius_mm', ...
    'radius_ratio', 'truth_inplane_ratio', 'truth_asymmetry', ...
    'truth_reverse_frac', 'truth_disturbance', 'n_points', 'n_valid', ...
    'valid_frac', 'mean_true_mms', 'bias_raw_mms', 'std_raw_mms', ...
    'rmse_raw_mms', 'nrmse_raw', 'bias_sg_mms', 'std_sg_mms', ...
    'rmse_sg_mms', 'nrmse_sg', 'n_full_points', ...
    'rmse_fill_full_mms', 'nrmse_fill_full', 'rmse_sg_full_mms', ...
    'nrmse_sg_full', 'flow_fill_ml_min', 'flow_sg_ml_min', ...
    'flow_true_ml_min', 'flow_transit_true_ml_min', ...
    'flow_error_fill_pct', 'flow_error_sg_pct'});
end


function msl = aggregate_seed_metrics(mseed)
%AGGREGATE_SEED_METRICS One row per physical slice, mean +/- std over seeds.
G = findgroups(mseed.dataset, mseed.slice_idx);
n = max(G);
name = strings(n, 1);
[slice, pos, nSeeds, nTot, nVal, frac, fracSeedStd, mtru] = ...
    deal(nan(n, 1));
[lumenRadius, radiusRatio, truthIp, truthAsym, truthReverse, truthDisturbance] = ...
    deal(nan(n, 1));
region = strings(n, 1);
[biasRaw, biasRawSeedStd, stdRaw, rmRaw, rmRawSeedStd, nrRaw, nrRawSeedStd] = ...
    deal(nan(n, 1));
[biasSg, biasSgSeedStd, stdSg, rmSg, rmSgSeedStd, nrSg, nrSgSeedStd] = ...
    deal(nan(n, 1));

for k = 1:n
    Q = mseed(G == k, :);
    name(k) = Q.dataset(1);
    slice(k) = Q.slice_idx(1);
    pos(k) = Q.slice_pos_mm(1);
    region(k) = Q.vessel_region(1);
    lumenRadius(k) = Q.lumen_radius_mm(1);
    radiusRatio(k) = Q.radius_ratio(1);
    truthIp(k) = Q.truth_inplane_ratio(1);
    truthAsym(k) = Q.truth_asymmetry(1);
    truthReverse(k) = Q.truth_reverse_frac(1);
    truthDisturbance(k) = Q.truth_disturbance(1);
    nSeeds(k) = height(Q);
    nTot(k) = mean(Q.n_points, 'omitnan');
    nVal(k) = mean(Q.n_valid, 'omitnan');
    frac(k) = mean(Q.valid_frac, 'omitnan');
    fracSeedStd(k) = repeat_std(Q.valid_frac);
    mtru(k) = mean(Q.mean_true_mms, 'omitnan');

    biasRaw(k) = mean(Q.bias_raw_mms, 'omitnan');
    biasRawSeedStd(k) = repeat_std(Q.bias_raw_mms);
    stdRaw(k) = mean(Q.std_raw_mms, 'omitnan');
    rmRaw(k) = mean(Q.rmse_raw_mms, 'omitnan');
    rmRawSeedStd(k) = repeat_std(Q.rmse_raw_mms);
    nrRaw(k) = mean(Q.nrmse_raw, 'omitnan');
    nrRawSeedStd(k) = repeat_std(Q.nrmse_raw);

    biasSg(k) = mean(Q.bias_sg_mms, 'omitnan');
    biasSgSeedStd(k) = repeat_std(Q.bias_sg_mms);
    stdSg(k) = mean(Q.std_sg_mms, 'omitnan');
    rmSg(k) = mean(Q.rmse_sg_mms, 'omitnan');
    rmSgSeedStd(k) = repeat_std(Q.rmse_sg_mms);
    nrSg(k) = mean(Q.nrmse_sg, 'omitnan');
    nrSgSeedStd(k) = repeat_std(Q.nrmse_sg);
end

msl = table(name, slice, pos, region, lumenRadius, radiusRatio, truthIp, ...
    truthAsym, truthReverse, truthDisturbance, nSeeds, nTot, nVal, frac, ...
    fracSeedStd, mtru, biasRaw, biasRawSeedStd, stdRaw, rmRaw, ...
    rmRawSeedStd, nrRaw, ...
    nrRawSeedStd, biasSg, biasSgSeedStd, stdSg, rmSg, rmSgSeedStd, nrSg, ...
    nrSgSeedStd, 'VariableNames', {'dataset', 'slice_idx', 'slice_pos_mm', ...
    'vessel_region', 'lumen_radius_mm', 'radius_ratio', ...
    'truth_inplane_ratio', 'truth_asymmetry', 'truth_reverse_frac', ...
    'truth_disturbance', 'n_seeds', 'n_points', 'n_valid', 'valid_frac', ...
    'valid_frac_seed_std', 'mean_true_mms', 'bias_raw_mms', ...
    'bias_raw_seed_std_mms', ...
    'std_raw_mms', 'rmse_raw_mms', 'rmse_raw_seed_std_mms', ...
    'nrmse_raw', 'nrmse_raw_seed_std', 'bias_sg_mms', ...
    'bias_sg_seed_std_mms', 'std_sg_mms', 'rmse_sg_mms', ...
    'rmse_sg_seed_std_mms', 'nrmse_sg', 'nrmse_sg_seed_std'});

fullCols = {'n_full_points', 'rmse_fill_full_mms', 'nrmse_fill_full', ...
    'rmse_sg_full_mms', 'nrmse_sg_full', 'flow_fill_ml_min', ...
    'flow_sg_ml_min', 'flow_true_ml_min', 'flow_transit_true_ml_min', ...
    'flow_error_fill_pct', 'flow_error_sg_pct'};
for j = 1:numel(fullCols)
    value = nan(n, 1);
    seedStd = nan(n, 1);
    for k = 1:n
        v = mseed.(fullCols{j})(G == k);
        value(k) = mean(v, 'omitnan');
        seedStd(k) = repeat_std(v);
    end
    msl.(fullCols{j}) = value;
    msl.([fullCols{j} '_seed_std']) = seedStd;
end

% Compatibility aliases for the current plotting functions. The dedicated
% raw/SG columns above remain authoritative.
msl.bias_mms = msl.bias_sg_mms;
msl.std_mms = msl.std_sg_mms;
msl.rmse_pt_mms = msl.rmse_raw_mms;
msl.nrmse_pt = msl.nrmse_raw;
end


function [r, nr] = rmse_pair(err, ref)
r = sqrt(mean(err.^2));
if ~isfinite(ref) || abs(ref) < eps
    nr = nan;
else
    nr = r / abs(ref);
end
end


function met = vessel_metrics(msl, name, vmax_mms)
%VESSEL_METRICS Descriptive dataset-level values; slices remain primary.
cols = {'valid_frac', 'bias_raw_mms', 'std_raw_mms', 'rmse_raw_mms', ...
    'nrmse_raw', 'bias_sg_mms', 'std_sg_mms', 'rmse_sg_mms', 'nrmse_sg', ...
    'nrmse_fill_full', 'nrmse_sg_full', 'flow_error_fill_pct', ...
    'flow_error_sg_pct'};
met = table(string(name), vmax_mms, height(msl), sum(msl.n_seeds), ...
    sum(msl.n_valid .* msl.n_seeds), 'VariableNames', ...
    {'dataset', 'vmax_mms', 'n_slices', 'n_seed_runs', 'n_valid_total'});
for k = 1:numel(cols)
    v = msl.(cols{k});
    met.([cols{k} '_mean']) = mean(v, 'omitnan');
    met.([cols{k} '_slice_std']) = repeat_std(v);
end
end


function s = repeat_std(v)
v = v(isfinite(v));
if numel(v) < 2
    s = nan;
else
    s = std(v);
end
end


function Q = representative_seed_points(P)
%REPRESENTATIVE_SEED_POINTS Keep one realisation per slice for diagnostics.
G = findgroups(P.dataset, P.slice_idx);
keep = false(height(P), 1);
for k = 1:max(G)
    idx = find(G == k);
    seeds = P.seed_base(idx);
    chosen = min(seeds);
    keep(idx(seeds == chosen)) = true;
end
Q = P(keep, :);
end


% -------------------------------------------------------------------- output

function s = med_neighbours(v)
%MED_NEIGHBOURS Running median over the +/-2 neighbouring elements,
% finite values only. NaN entries stay NaN: a slice that could not fit
% its own profile does not inherit one from its neighbours.
s = v;
for i = 1:numel(v)
    if ~isfinite(v(i))
        continue;
    end
    w = v(max(1, i-2):min(numel(v), i+2));
    s(i) = median(w(isfinite(w)));
end
end

function [ord, pos] = sort_snaps(snapPos)
%SORT_SNAPS Acquisition order of the snapshots along the scan.
[pos, ord] = sort(snapPos);
end


function [ok, ratio] = tissue_present(snap, opts)
%TISSUE_PRESENT Is there enough echo outside the lumen to localize on?
% Clean simulations seed scatterers in the lumen only; their outside is
% numerically empty and the correlation would lock onto nothing.
[X, Z] = meshgrid(snap.x, snap.z);
rr = hypot(X - snap.xc, Z - snap.zc);
vIn = mean(snap.E1(rr < 0.8 * snap.radius));
vOut = mean(snap.E1(rr > opts.loc_margin * snap.radius));
ratio = double(vOut) / max(double(vIn), eps);
ok = ratio > 0.05;
end


function v = powerlaw_v(r, q, R)
%POWERLAW_V Blunt velocity profile v0 * (1 - (r/R)^n).
v = q(1) * (1 - min(max(r / R, 0), 1) .^ abs(q(2)));
end

function q = fit_powerlaw(r, v, w, R, q0)
%FIT_POWERLAW Robust weighted fit of the blunt profile, soft-L1 loss.
% The loss saturates for outliers past f_scale (mm/s), so near-wall junk
% points cannot drag the exponent the way plain least squares lets them.
fs = 50;
obj = @(q) sum(2 * fs^2 * ...
    (sqrt(1 + ((w .* (powerlaw_v(r, q, R) - v)) / fs) .^ 2) - 1));
q = fminsearch(obj, q0, ...
    optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 4000));
end

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
