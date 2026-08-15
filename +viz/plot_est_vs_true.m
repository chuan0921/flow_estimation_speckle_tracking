function [fig, ax] = plot_est_vs_true(P, varargin)
%PLOT_EST_VS_TRUE Estimated against true v_y, all conditions on one axis.
%
%   viz.plot_est_vs_true(P)
%   viz.plot_est_vs_true(P, 'HideOutliers', false)
%
% P is the pooled points table (points.csv concatenated across datasets).
% One marker per grid point per slice, coloured by dataset, against the 1:1
% line. Departure from the diagonal is the bias; vertical scatter is the
% precision. Because every condition is on the same axis, the plot also shows
% the velocity ceiling directly: points stop appearing above the speed the
% lag scan can still resolve.
%
% Name-value options
%   Title         ''     figure title prefix
%   HideOutliers  true   drop rows flagged by the normalized median test
%
% Returns the figure and axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_est_vs_true';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.addParameter('HideOutliers', true, @(v) islogical(v) && isscalar(v));
p.parse(varargin{:});
o = p.Results;

m = isfinite(P.Vy_mms) & isfinite(P.VyTrue_mms);
if o.HideOutliers && ismember('outlier', P.Properties.VariableNames)
    m = m & ~P.outlier;
end
P = P(m, :);
if isempty(P)
    error('viz:plot_est_vs_true:noPoints', 'No finite paired points in P.');
end

[names, speed] = order_datasets(P.dataset);
C = lines(max(numel(names), 3));

fig = figure('Color', 'w', 'Name', 'estimate vs truth', ...
    'Position', [100 100 640 600]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

lim = [0, max([P.Vy_mms; P.VyTrue_mms]) * 1.05];
lim(1) = min([0; P.Vy_mms; P.VyTrue_mms]) * 1.05;
h = gobjects(1, numel(names));
lbl = cell(1, numel(names));
for k = 1:numel(names)
    Q = P(P.dataset == names(k), :);
    h(k) = scatter(ax, Q.VyTrue_mms, Q.Vy_mms, 14, C(k, :), 'filled', ...
        'MarkerFaceAlpha', 0.35);
    r = corr_safe(Q.VyTrue_mms, Q.Vy_mms);
    lbl{k} = sprintf('%s  (n=%d, r=%.3f)', label_for(names(k), speed(k)), ...
        height(Q), r);
end
hd = plot(ax, lim, lim, '--k', 'LineWidth', 1.4);
hold(ax, 'off');

axis(ax, 'equal');
xlim(ax, lim);
ylim(ax, lim);
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, 'true v_y [mm/s]');
ylabel(ax, 'estimated v_y [mm/s]');
legend(ax, [h hd], [lbl, {'1:1'}], 'Location', 'northwest', 'Box', 'off');
title(ax, strtrim(sprintf('%s  estimate vs truth', char(o.Title))), 'Interpreter', 'none');
end


function r = corr_safe(a, b)
if numel(a) < 2 || std(a) == 0 || std(b) == 0
    r = nan;
    return;
end
c = corrcoef(a, b);
r = c(1, 2);
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


function L = label_for(name, speed)
if isfinite(speed)
    L = sprintf('%g cm/s', speed);
else
    L = char(name);
end
end
