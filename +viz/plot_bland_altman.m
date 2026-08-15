function [fig, ax] = plot_bland_altman(P, varargin)
%PLOT_BLAND_ALTMAN Agreement between estimate and truth.
%
%   viz.plot_bland_altman(P)
%   viz.plot_bland_altman(P, 'Relative', true)
%
% P is the pooled points table (points.csv concatenated across datasets).
% Difference against mean, with the bias and the +/-1.96 std limits of
% agreement drawn per dataset-independent pooling. A trend in the cloud means
% the error scales with velocity rather than being a fixed offset, which is
% the distinction the RMSE number alone cannot make.
%
% Name-value options
%   Title         ''     figure title prefix
%   HideOutliers  true   drop rows flagged by the normalized median test
%   Relative      false  plot the difference as a percentage of the mean
%
% Returns the figure and axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_bland_altman';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.addParameter('HideOutliers', true, @(v) islogical(v) && isscalar(v));
p.addParameter('Relative', false, @(v) islogical(v) && isscalar(v));
p.parse(varargin{:});
o = p.Results;

m = isfinite(P.Vy_mms) & isfinite(P.VyTrue_mms);
if o.HideOutliers && ismember('outlier', P.Properties.VariableNames)
    m = m & ~P.outlier;
end
P = P(m, :);
if isempty(P)
    error('viz:plot_bland_altman:noPoints', 'No finite paired points in P.');
end

avg = (P.Vy_mms + P.VyTrue_mms) / 2;
dif = P.Vy_mms - P.VyTrue_mms;
unit = 'mm/s';
if o.Relative
    ref = avg;
    ref(ref == 0) = nan;
    dif = 100 * dif ./ ref;
    unit = '%';
end
ok = isfinite(avg) & isfinite(dif);
avg = avg(ok);
dif = dif(ok);
P = P(ok, :);

bias = mean(dif);
sd = std(dif);
loa = bias + [-1.96, 1.96] * sd;

[names, speed] = order_datasets(P.dataset);
C = lines(max(numel(names), 3));

fig = figure('Color', 'w', 'Name', 'Bland-Altman', ...
    'Position', [100 100 720 520]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

h = gobjects(1, numel(names));
lbl = cell(1, numel(names));
for k = 1:numel(names)
    s = P.dataset == names(k);
    h(k) = scatter(ax, avg(s), dif(s), 14, C(k, :), 'filled', ...
        'MarkerFaceAlpha', 0.35);
    lbl{k} = label_for(names(k), speed(k));
end

xr = [min(avg), max(avg)];
plot(ax, xr, [bias bias], '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.6);
plot(ax, xr, [loa(1) loa(1)], '--', 'Color', [0.75 0.2 0.2], 'LineWidth', 1.3);
plot(ax, xr, [loa(2) loa(2)], '--', 'Color', [0.75 0.2 0.2], 'LineWidth', 1.3);
text(ax, xr(2), bias, sprintf(' bias %.1f', bias), 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right', 'Color', [0.15 0.15 0.15]);
text(ax, xr(2), loa(2), sprintf(' +1.96sd %.1f', loa(2)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
    'Color', [0.75 0.2 0.2]);
text(ax, xr(2), loa(1), sprintf(' -1.96sd %.1f', loa(1)), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
    'Color', [0.75 0.2 0.2]);
hold(ax, 'off');

grid(ax, 'on');
box(ax, 'on');
xlabel(ax, sprintf('mean of estimate and truth [mm/s]'));
ylabel(ax, sprintf('estimate - truth [%s]', unit));
legend(ax, h, lbl, 'Location', 'best', 'Box', 'off');
title(ax, strtrim(sprintf('%s  Bland-Altman (n = %d)', char(o.Title), numel(dif))), 'Interpreter', 'none');
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
