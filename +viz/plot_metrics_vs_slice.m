function [fig, ax] = plot_metrics_vs_slice(msl, varargin)
%PLOT_METRICS_VS_SLICE Per-slice error metrics along the vessel.
%
%   viz.plot_metrics_vs_slice(msl)
%   viz.plot_metrics_vs_slice(msl, 'Title', 'sim_v080cms')
%
% msl is the per-slice metrics table written by main
% (metrics_slice.csv). Three stacked panels against elevational position:
% NRMSE (raw and SG-filtered), bias with a +/-1 std residual band, and the
% fraction of grid points the estimator answered. Reading them together
% separates "the estimate is wrong here" from "there was nothing to
% estimate here".
%
% Name-value options
%   Title    ''    figure title prefix
%
% Returns the figure and the three axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_metrics_vs_slice';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.parse(varargin{:});
o = p.Results;

[y, k] = sort(msl.slice_pos_mm);
msl = msl(k, :);
blue = [0.20 0.45 0.75];
red = [0.85 0.33 0.10];

fig = figure('Color', 'w', 'Name', 'metrics vs slice', ...
    'Position', [100 100 720 640]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax(1) = nexttile(tl);
hold(ax(1), 'on');
h1 = plot(ax(1), y, 100 * msl.nrmse_raw, '-o', 'Color', blue, ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'MarkerFaceColor', blue);
h2 = plot(ax(1), y, 100 * msl.nrmse_sg, '-s', 'Color', red, ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'MarkerFaceColor', red);
hold(ax(1), 'off');
grid(ax(1), 'on');
box(ax(1), 'on');
ylabel(ax(1), 'NRMSE [%]');
legend(ax(1), [h1 h2], {'raw', 'Savitzky-Golay'}, 'Location', 'best', ...
    'Box', 'off');
title(ax(1), strtrim(sprintf('%s  per-slice metrics', char(o.Title))), 'Interpreter', 'none');
set(ax(1), 'XTickLabel', []);

ax(2) = nexttile(tl);
hold(ax(2), 'on');
yline(ax(2), 0, 'k-', 'LineWidth', 1);
ok = isfinite(msl.bias_mms) & isfinite(msl.std_mms);
if any(ok)
    fill(ax(2), [y(ok); flipud(y(ok))], ...
        [msl.bias_mms(ok) - msl.std_mms(ok); ...
         flipud(msl.bias_mms(ok) + msl.std_mms(ok))], blue, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
plot(ax(2), y, msl.bias_mms, '-o', 'Color', blue, 'LineWidth', 1.4, ...
    'MarkerSize', 4, 'MarkerFaceColor', blue);
hold(ax(2), 'off');
grid(ax(2), 'on');
box(ax(2), 'on');
ylabel(ax(2), 'bias +/- 1 std [mm/s]');
set(ax(2), 'XTickLabel', []);

ax(3) = nexttile(tl);
bar(ax(3), y, 100 * msl.valid_frac, 0.8, 'FaceColor', [0.6 0.64 0.7], ...
    'EdgeColor', 'none');
grid(ax(3), 'on');
box(ax(3), 'on');
ylim(ax(3), [0 100]);
xlabel(ax(3), 'slice position [mm]');
ylabel(ax(3), 'valid points [%]');
linkaxes(ax, 'x');
end
