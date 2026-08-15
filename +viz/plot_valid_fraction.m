function [fig, ax] = plot_valid_fraction(msl, varargin)
%PLOT_VALID_FRACTION How much of the lumen the estimator could answer.
%
%   viz.plot_valid_fraction(msl)
%
% msl is the pooled per-slice metrics table (metrics_slice_all.csv). Each
% dataset gets a box-like summary across its slices: the across-slice mean
% with a +/-1 std bar, over the individual slice values.
%
% This is the companion to plot_nrmse_vs_speed and has to be read with it.
% NRMSE is computed only over points that produced an estimate, so a
% condition can post a flattering NRMSE precisely because the hard points
% dropped out. This figure is where that shows up.
%
% Name-value options
%   Title    ''    figure title prefix
%
% Returns the figure and axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_valid_fraction';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.parse(varargin{:});
o = p.Results;

[names, speed] = order_datasets(msl.dataset);
n = numel(names);

fig = figure('Color', 'w', 'Name', 'valid fraction', ...
    'Position', [100 100 720 480]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

for k = 1:n
    v = 100 * msl.valid_frac(msl.dataset == names(k));
    jitter = (rand(size(v)) - 0.5) * 0.28;
    scatter(ax, k + jitter, v, 22, [0.62 0.66 0.72], 'filled', ...
        'MarkerFaceAlpha', 0.6);
    m = mean(v, 'omitnan');
    s = std(v, 'omitnan');
    errorbar(ax, k, m, s, 'Color', [0.20 0.45 0.75], 'LineWidth', 1.6, ...
        'CapSize', 12, 'Marker', 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', [0.20 0.45 0.75]);
    text(ax, k, 102, sprintf('%.0f%%', m), 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'Color', [0.25 0.25 0.25]);
end

hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [0.4, n + 0.6]);
ylim(ax, [0 110]);
set(ax, 'XTick', 1:n, 'XTickLabel', label_for(names, speed), 'YTick', 0:20:100);
xlabel(ax, 'peak velocity');
ylabel(ax, 'points with a valid v_y [%]');
title(ax, strtrim(sprintf('%s  estimator yield per slice', char(o.Title))), 'Interpreter', 'none');
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
