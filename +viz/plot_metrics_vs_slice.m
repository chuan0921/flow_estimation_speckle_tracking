function [fig, ax] = plot_metrics_vs_slice(msl, varargin)
%PLOT_METRICS_VS_SLICE Region-aware error summary along the vessel.
%
%   viz.plot_metrics_vs_slice(msl)
%   viz.plot_metrics_vs_slice(msl, 'Regions', regions, 'Title', name)
%
% Raw and 2-D SG metrics are plotted per physical slice. Error bars are the
% across-seed standard deviation and therefore appear only where n_seeds>1.
% Regions supplies the full-vessel classification, even when msl contains a
% subset of slices.

p = inputParser;
p.FunctionName = 'viz.plot_metrics_vs_slice';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.addParameter('Regions', table(), @istable);
p.parse(varargin{:});
o = p.Results;

[y, order] = sort(msl.slice_pos_mm);
msl = msl(order, :);
regions = normalise_regions(msl, o.Regions);

blue = [0.16 0.40 0.68];
red = [0.82 0.24 0.15];
gray = [0.35 0.38 0.42];
green = [0.12 0.52 0.34];
gold = [0.86 0.52 0.08];

fig = figure('Color', 'w', 'Name', 'vessel summary', ...
    'Position', [80 30 1080 1180]);
tl = tiledlayout(fig, 6, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% NRMSE -------------------------------------------------------------------
ax(1) = nexttile(tl);
hold(ax(1), 'on');
hRaw = plot(ax(1), y, 100 * msl.nrmse_raw, '-o', 'Color', blue, ...
    'LineWidth', 1.5, 'MarkerSize', 4.5, 'MarkerFaceColor', blue);
hSg = plot(ax(1), y, 100 * msl.nrmse_sg, '-s', 'Color', red, ...
    'LineWidth', 1.7, 'MarkerSize', 4.5, 'MarkerFaceColor', red);
seed_errorbars(ax(1), y, 100 * msl.nrmse_raw, ...
    100 * column_or_nan(msl, 'nrmse_raw_seed_std'), msl.n_seeds, blue);
seed_errorbars(ax(1), y, 100 * msl.nrmse_sg, ...
    100 * column_or_nan(msl, 'nrmse_sg_seed_std'), msl.n_seeds, red);
hold(ax(1), 'off');
ylabel(ax(1), 'NRMSE [%]');
legend(ax(1), [hRaw hSg], {'Raw', 'After SG (2-D)'}, ...
    'Location', 'northwest', 'Box', 'off', 'Orientation', 'horizontal');
title(ax(1), strtrim(sprintf('%s  vessel-region error summary', ...
    char(o.Title))), 'Interpreter', 'none');

% Full-lumen NRMSE --------------------------------------------------------
ax(2) = nexttile(tl);
hold(ax(2), 'on');
hFill = plot(ax(2), y, 100 * msl.nrmse_fill_full, '-o', 'Color', gold, ...
    'LineWidth', 1.5, 'MarkerSize', 4.5, 'MarkerFaceColor', gold);
hSgFull = plot(ax(2), y, 100 * msl.nrmse_sg_full, '-s', 'Color', red, ...
    'LineWidth', 1.7, 'MarkerSize', 4.5, 'MarkerFaceColor', red);
seed_errorbars(ax(2), y, 100 * msl.nrmse_fill_full, ...
    100 * column_or_nan(msl, 'nrmse_fill_full_seed_std'), ...
    msl.n_seeds, gold);
seed_errorbars(ax(2), y, 100 * msl.nrmse_sg_full, ...
    100 * column_or_nan(msl, 'nrmse_sg_full_seed_std'), ...
    msl.n_seeds, red);
hold(ax(2), 'off');
ylabel(ax(2), 'Full-lumen NRMSE [%]');
legend(ax(2), [hFill hSgFull], {'No-slip fill', 'After SG (2-D)'}, ...
    'Location', 'northwest', 'Box', 'off', 'Orientation', 'horizontal');

% Flow error ---------------------------------------------------------------
ax(3) = nexttile(tl);
hold(ax(3), 'on');
yline(ax(3), 0, 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
plot(ax(3), y, msl.flow_error_fill_pct, '-o', 'Color', gold, ...
    'LineWidth', 1.5, 'MarkerSize', 4.5, 'MarkerFaceColor', gold);
plot(ax(3), y, msl.flow_error_sg_pct, '-s', 'Color', red, ...
    'LineWidth', 1.7, 'MarkerSize', 4.5, 'MarkerFaceColor', red);
seed_errorbars(ax(3), y, msl.flow_error_fill_pct, ...
    column_or_nan(msl, 'flow_error_fill_pct_seed_std'), msl.n_seeds, gold);
seed_errorbars(ax(3), y, msl.flow_error_sg_pct, ...
    column_or_nan(msl, 'flow_error_sg_pct_seed_std'), msl.n_seeds, red);
hold(ax(3), 'off');
ylabel(ax(3), 'Flow error [%]');

% Bias --------------------------------------------------------------------
ax(4) = nexttile(tl);
hold(ax(4), 'on');
yline(ax(4), 0, 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
plot(ax(4), y, msl.bias_raw_mms, '-o', 'Color', blue, ...
    'LineWidth', 1.5, 'MarkerSize', 4.5, 'MarkerFaceColor', blue);
plot(ax(4), y, msl.bias_sg_mms, '-s', 'Color', red, ...
    'LineWidth', 1.7, 'MarkerSize', 4.5, 'MarkerFaceColor', red);
seed_errorbars(ax(4), y, msl.bias_raw_mms, ...
    column_or_nan(msl, 'bias_raw_seed_std_mms'), msl.n_seeds, blue);
seed_errorbars(ax(4), y, msl.bias_sg_mms, ...
    column_or_nan(msl, 'bias_sg_seed_std_mms'), msl.n_seeds, red);
hold(ax(4), 'off');
ylabel(ax(4), 'Bias [mm/s]');

% Estimator yield ----------------------------------------------------------
ax(5) = nexttile(tl);
hold(ax(5), 'on');
plot(ax(5), y, 100 * msl.valid_frac, '-o', 'Color', gray, ...
    'LineWidth', 1.5, 'MarkerSize', 4.5, 'MarkerFaceColor', gray);
seed_errorbars(ax(5), y, 100 * msl.valid_frac, ...
    100 * column_or_nan(msl, 'valid_frac_seed_std'), ...
    msl.n_seeds, gray);
multi = find(msl.n_seeds > 1);
for k = multi.'
    text(ax(5), y(k), 103, sprintf('n=%d', msl.n_seeds(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9, 'Color', [0.2 0.2 0.2]);
end
hold(ax(5), 'off');
ylim(ax(5), [0 112]);
ylabel(ax(5), 'Valid ROI [%]');

% Geometry and truth disturbance -----------------------------------------
ax(6) = nexttile(tl);
hold(ax(6), 'on');
yyaxis(ax(6), 'left');
hRadius = plot(ax(6), regions.slice_pos_mm, regions.lumen_radius_mm, ...
    '-o', 'Color', green, 'LineWidth', 1.6, 'MarkerSize', 3.5, ...
    'MarkerFaceColor', green);
ylabel(ax(6), 'Lumen radius [mm]');
yyaxis(ax(6), 'right');
hDisturbance = plot(ax(6), regions.slice_pos_mm, ...
    regions.truth_disturbance, '-', 'Color', [0.52 0.24 0.58], ...
    'LineWidth', 1.7);
threshold = column_or_nan(regions, 'disturbance_threshold');
hThreshold = plot(ax(6), regions.slice_pos_mm, threshold, '--', ...
    'Color', [0.45 0.45 0.45], 'LineWidth', 1.1);
ylabel(ax(6), 'Truth disturbance score');
xlabel(ax(6), 'Slice position [mm]');
legend(ax(6), [hRadius hDisturbance hThreshold], ...
    {'Lumen radius', 'Truth disturbance', 'Disturbance threshold'}, ...
    'Location', 'northwest', 'Box', 'off', 'Orientation', 'horizontal');
hold(ax(6), 'off');

% Shared vessel-region context -------------------------------------------
for k = 1:6
    grid(ax(k), 'on');
    box(ax(k), 'on');
    set(ax(k), 'FontSize', 11, 'Layer', 'top');
    if k < 6
        set(ax(k), 'XTickLabel', []);
    end
    if k == 6
        yyaxis(ax(k), 'left');
    end
    shade_regions(ax(k), regions);
end
label_regions(ax(1), regions);
linkaxes(ax, 'x');
end


function seed_errorbars(ax, x, value, spread, nSeeds, color)
ok = nSeeds > 1 & isfinite(value) & isfinite(spread);
if any(ok)
    errorbar(ax, x(ok), value(ok), spread(ok), 'LineStyle', 'none', ...
        'Color', color, 'LineWidth', 1.2, 'CapSize', 7, ...
        'HandleVisibility', 'off');
end
end


function v = column_or_nan(T, name)
if ismember(name, T.Properties.VariableNames)
    v = T.(name);
else
    v = nan(height(T), 1);
end
end


function regions = normalise_regions(msl, regions)
need = {'slice_pos_mm', 'lumen_radius_mm', 'truth_disturbance', ...
    'vessel_region'};
if isempty(regions)
    if all(ismember(need, msl.Properties.VariableNames))
        regions = unique(msl(:, need), 'rows', 'stable');
    else
        regions = table(msl.slice_pos_mm, nan(height(msl), 1), ...
            nan(height(msl), 1), repmat("uniform vessel", height(msl), 1), ...
            'VariableNames', need);
    end
end
[~, order] = sort(regions.slice_pos_mm);
regions = regions(order, :);
if ~ismember('disturbance_threshold', regions.Properties.VariableNames)
    regions.disturbance_threshold = nan(height(regions), 1);
end
end


function shade_regions(ax, regions)
[edges, starts, stops] = region_runs(regions);
yl = ylim(ax);
holdState = ishold(ax);
hold(ax, 'on');
for k = 1:numel(starts)
    color = region_color(regions.vessel_region(starts(k)));
    patch(ax, [edges(starts(k)) edges(stops(k) + 1) ...
        edges(stops(k) + 1) edges(starts(k))], ...
        [yl(1) yl(1) yl(2) yl(2)], color, ...
        'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
if ~holdState
    hold(ax, 'off');
end
xlim(ax, [edges(1) edges(end)]);
end


function label_regions(ax, regions)
[edges, starts, stops] = region_runs(regions);
yl = ylim(ax);
dy = diff(yl);
for k = 1:numel(starts)
    left = edges(starts(k));
    right = edges(stops(k) + 1);
    text(ax, (left + right) / 2, yl(2) - 0.025 * dy, ...
        display_region(regions.vessel_region(starts(k))), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.18 0.18 0.18], ...
        'Interpreter', 'none', 'Clipping', 'on');
end
end


function [edges, starts, stops] = region_runs(regions)
y = regions.slice_pos_mm;
if isscalar(y)
    edges = [y - 0.5; y + 0.5];
else
    mid = (y(1:end-1) + y(2:end)) / 2;
    edges = [y(1) - (mid(1) - y(1)); mid; ...
        y(end) + (y(end) - mid(end))];
end
change = [true; regions.vessel_region(2:end) ~= ...
    regions.vessel_region(1:end-1)];
starts = find(change);
stops = [starts(2:end) - 1; height(regions)];
end


function color = region_color(name)
switch string(name)
    case "pre-stenosis"
        color = [0.70 0.74 0.78];
    case "narrowing"
        color = [0.95 0.72 0.20];
    case "stenosis throat"
        color = [0.86 0.25 0.18];
    case "post-stenosis"
        color = [0.30 0.65 0.42];
    case "disturbed flow"
        color = [0.48 0.32 0.62];
    case "recovery"
        color = [0.20 0.55 0.74];
    otherwise
        color = [0.72 0.74 0.76];
end
end


function label = display_region(name)
switch string(name)
    case "pre-stenosis"
        label = 'Pre-stenosis';
    case "narrowing"
        label = 'Narrowing';
    case "stenosis throat"
        label = 'Stenosis throat';
    case "post-stenosis"
        label = 'Post-stenosis';
    case "disturbed flow"
        label = 'Disturbed flow';
    case "recovery"
        label = 'Recovery';
    otherwise
        label = 'Uniform vessel';
end
end
