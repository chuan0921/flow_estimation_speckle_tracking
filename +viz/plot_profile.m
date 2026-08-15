function [fig, ax] = plot_profile(prof, varargin)
%PLOT_PROFILE Radial blood-flow velocity profile: raw, SG-filtered, truth.
%
%   viz.plot_profile(prof)
%   viz.plot_profile(prof, 'Title', 'sim_v080cms')
%
% prof is the profile table written by main (profile.csv). Every slice
% shares the same radial lattice, so the slices are overlaid as faint points
% and summarised by the across-slice mean with a +/-1 std band. The analytic
% or CFD truth is drawn on top as a dashed reference.
%
% Name-value options
%   Title    ''      figure title prefix
%   Band     true    shade the +/-1 std across-slice spread
%
% Returns the figure and axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_profile';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.addParameter('Band', true, @(v) islogical(v) && isscalar(v));
p.parse(varargin{:});
o = p.Results;

[r, raw_m, raw_s] = collapse(prof, 'vy_raw');
[~, sg_m] = collapse(prof, 'vy_sg');
[~, tru_m] = collapse(prof, 'vy_true');

fig = figure('Color', 'w', 'Name', 'radial velocity profile', ...
    'Position', [100 100 720 480]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

h = gobjects(1, 0);
lbl = {};

% Individual slices, so the reader can see the spread the band summarises.
sc = scatter(ax, prof.r_mm, prof.vy_raw, 8, [0.7 0.75 0.82], 'filled');
sc.MarkerFaceAlpha = 0.45;
h(end + 1) = sc;
lbl{end + 1} = 'per-slice raw';

if o.Band && any(isfinite(raw_s))
    ok = isfinite(raw_m) & isfinite(raw_s);
    fill(ax, [r(ok); flipud(r(ok))], ...
        [raw_m(ok) - raw_s(ok); flipud(raw_m(ok) + raw_s(ok))], ...
        [0.20 0.45 0.75], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

h(end + 1) = plot(ax, r, raw_m, '-o', 'Color', [0.20 0.45 0.75], ...
    'LineWidth', 1.4, 'MarkerSize', 3.5, 'MarkerFaceColor', [0.20 0.45 0.75]);
lbl{end + 1} = 'raw (mean over slices)';

h(end + 1) = plot(ax, r, sg_m, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
lbl{end + 1} = 'Savitzky-Golay';

h(end + 1) = plot(ax, r, tru_m, '--k', 'LineWidth', 1.6);
lbl{end + 1} = 'truth';

hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, 'radial position r [mm]');
ylabel(ax, 'v_y [mm/s]');
legend(ax, h, lbl, 'Location', 'southwest', 'Box', 'off');
title(ax, strtrim(sprintf('%s  blood flow velocity profile (%d slices)', ...
    char(o.Title), numel(unique(prof.slice_idx)))), 'Interpreter', 'none');
end


function [r, m, s] = collapse(prof, col)
%COLLAPSE Across-slice mean and std at each radius.
[r, ~, g] = unique(round(prof.r_mm, 6));
v = prof.(col);
m = accumarray(g, v, [], @(x) mean(x, 'omitnan'));
s = accumarray(g, v, [], @(x) std(x, 'omitnan'));
m(accumarray(g, double(isfinite(v)), [], @sum) == 0) = nan;
end
