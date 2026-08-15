function [fig, ax] = plot_nrmse_vs_speed(msl, varargin)
%PLOT_NRMSE_VS_SPEED Headline figure: accuracy against flow speed.
%
%   viz.plot_nrmse_vs_speed(msl)
%
% msl is the pooled per-slice metrics table (metrics_slice_all.csv), one row
% per slice per dataset. Bars are the across-slice mean NRMSE for the raw and
% Savitzky-Golay profiles; error bars are the across-slice standard
% deviation, i.e. how reproducible the accuracy is along the vessel rather
% than how large the residual is within one slice.
%
% Datasets are ordered by peak velocity, taken from the trailing number in
% the folder name (sim_v080cms -> 80 cm/s) and falling back to alphabetical.
%
% Name-value options
%   Title    ''    figure title prefix
%
% Returns the figure and axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_nrmse_vs_speed';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.parse(varargin{:});
o = p.Results;

[names, speed] = order_datasets(msl.dataset);
n = numel(names);
M = nan(n, 2);
E = nan(n, 2);
for k = 1:n
    Q = msl(msl.dataset == names(k), :);
    M(k, :) = [mean(100 * Q.nrmse_raw, 'omitnan'), ...
        mean(100 * Q.nrmse_sg, 'omitnan')];
    E(k, :) = [std(100 * Q.nrmse_raw, 'omitnan'), ...
        std(100 * Q.nrmse_sg, 'omitnan')];
end

fig = figure('Color', 'w', 'Name', 'NRMSE vs speed', ...
    'Position', [100 100 720 480]);
ax = axes(fig); %#ok<LAXES>
x = (1:n).';
off = [-0.18 0.18];
w = 0.32;
cols = [0.20 0.45 0.75; 0.85 0.33 0.10];
b = gobjects(1, 2);
hold(ax, 'on');
for j = 1:2
    xj = x + off(j);
    b(j) = bar(ax, xj, M(:, j), w, 'FaceColor', cols(j, :), ...
        'EdgeColor', 'none');
    errorbar(ax, xj, M(:, j), E(:, j), 'k', 'LineStyle', 'none', ...
        'LineWidth', 1, 'CapSize', 8);
end
hold(ax, 'off');

grid(ax, 'on');
box(ax, 'on');
set(ax, 'XTick', 1:n, 'XTickLabel', label_for(names, speed));
xlabel(ax, 'peak velocity');
ylabel(ax, 'NRMSE [%]   (mean +/- std across slices)');
legend(ax, b, {'raw profile', 'Savitzky-Golay'}, 'Location', 'northwest', ...
    'Box', 'off');
title(ax, strtrim(sprintf('%s  through-plane accuracy vs flow speed', ...
    char(o.Title))), 'Interpreter', 'none');

% A condition with no valid slices would otherwise draw as a silent gap.
dead = all(~isfinite(M), 2);
if any(dead)
    text(ax, find(dead), zeros(nnz(dead), 1), 'no valid points', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontAngle', 'italic', 'Color', [0.5 0.5 0.5]);
end
end


function [names, speed] = order_datasets(col)
names = unique(col, 'stable');
speed = nan(size(names));
for k = 1:numel(names)
    t = regexp(char(names(k)), '(\d+)\s*cms?$', 'tokens', 'once');
    if ~isempty(t)
        speed(k) = str2double(t{1});
    end
end
if all(isfinite(speed))
    [speed, k] = sort(speed);
    names = names(k);
else
    [names, k] = sort(names);
    speed = speed(k);
end
end


function L = label_for(names, speed)
L = cell(1, numel(names));
for k = 1:numel(names)
    if isfinite(speed(k))
        L{k} = sprintf('%g cm/s', speed(k));
    else
        L{k} = char(names(k));
    end
end
end
