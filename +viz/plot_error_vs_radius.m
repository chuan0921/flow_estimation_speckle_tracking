function [fig, ax] = plot_error_vs_radius(prof, varargin)
%PLOT_ERROR_VS_RADIUS Bias and spread of the vy error against radius.
%
%   viz.plot_error_vs_radius(prof)
%   viz.plot_error_vs_radius(prof, 'Title', 'sim_v080cms', 'Relative', true)
%
% prof is the profile table written by main (profile.csv). The upper
% panel shows the raw and SG-filtered bias with a +/-1 std band taken across
% slices; the lower panel shows how many slices contributed at each radius,
% which is what makes the near-wall bias interpretable -- the outer radii are
% both the hardest to track and the thinnest sampled.
%
% Name-value options
%   Title      ''     figure title prefix
%   Relative   false  plot the error as a percentage of the local truth
%
% Returns the figure and the two axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_error_vs_radius';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.addParameter('Relative', false, @(v) islogical(v) && isscalar(v));
p.parse(varargin{:});
o = p.Results;

e_raw = prof.vy_raw - prof.vy_true;
e_sg = prof.vy_sg - prof.vy_true;
unit = 'mm/s';
if o.Relative
    ref = prof.vy_true;
    ref(ref == 0) = nan;
    e_raw = 100 * e_raw ./ ref;
    e_sg = 100 * e_sg ./ ref;
    unit = '%';
end

[r, braw, sraw] = collapse(prof.r_mm, e_raw);
[~, bsg, ssg] = collapse(prof.r_mm, e_sg);
[~, nval] = collapse(prof.r_mm, double(isfinite(prof.vy_raw)));

fig = figure('Color', 'w', 'Name', 'error vs radius', ...
    'Position', [100 100 720 560]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax(1) = nexttile(tl, 1, [2 1]);
hold(ax(1), 'on');
yline(ax(1), 0, 'k-', 'LineWidth', 1);
band(ax(1), r, braw, sraw, [0.20 0.45 0.75]);
band(ax(1), r, bsg, ssg, [0.85 0.33 0.10]);
h1 = plot(ax(1), r, braw, '-o', 'Color', [0.20 0.45 0.75], 'LineWidth', 1.4, ...
    'MarkerSize', 3.5, 'MarkerFaceColor', [0.20 0.45 0.75]);
h2 = plot(ax(1), r, bsg, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
hold(ax(1), 'off');
grid(ax(1), 'on');
box(ax(1), 'on');
ylabel(ax(1), sprintf('v_y error [%s]', unit));
legend(ax(1), [h1 h2], {'raw', 'Savitzky-Golay'}, 'Location', 'northwest', ...
    'Box', 'off');
title(ax(1), strtrim(sprintf('%s  bias +/- 1 std across slices', char(o.Title))), 'Interpreter', 'none');
set(ax(1), 'XTickLabel', []);

ax(2) = nexttile(tl);
bar(ax(2), r, nval, 1, 'FaceColor', [0.6 0.64 0.7], 'EdgeColor', 'none');
grid(ax(2), 'on');
box(ax(2), 'on');
xlabel(ax(2), 'radial position r [mm]');
ylabel(ax(2), 'slices');
linkaxes(ax, 'x');
xlim(ax(2), [min(r) - 0.05, max(r) + 0.05]);
end


function band(ax, r, m, s, c)
ok = isfinite(m) & isfinite(s);
if ~any(ok)
    return;
end
fill(ax, [r(ok); flipud(r(ok))], [m(ok) - s(ok); flipud(m(ok) + s(ok))], c, ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
end


function [r, m, s] = collapse(rv, v)
[r, ~, g] = unique(round(rv, 6));
m = accumarray(g, v, [], @(x) mean(x, 'omitnan'));
s = accumarray(g, v, [], @(x) std(x, 'omitnan'));
m(accumarray(g, double(isfinite(v)), [], @sum) == 0) = nan;
end
