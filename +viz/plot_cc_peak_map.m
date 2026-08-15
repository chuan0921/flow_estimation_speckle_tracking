function [fig, ax] = plot_cc_peak_map(P, varargin)
%PLOT_CC_PEAK_MAP Per-slice ROI map of through-plane CC peaks.
%
%   viz.plot_cc_peak_map(points)
%   viz.plot_cc_peak_map(points, 'Title', 'sim_stenosis50_80cms')
%
% points is the table returned by main, or a single-slice table from
% src.estimate_slice. Each tile is one slice. ROI locations are plotted in the x-z plane;
% marker colour is the selected peak correlation, and grey x markers are
% ROI locations where no through-plane peak passed the gate.
%
% Name-value options
%   Title          ''       figure title prefix
%   CCColumn       'ccA'    table column used for peak correlation colour
%   LagColumn      'lag_peak_frames'
%   HideOutliers   false    hide rows flagged by the NMT outlier test
%   MarkerSize     28       scatter marker area [pt^2]
%   MaxSlices      Inf      plot at most this many slices, evenly sampled
%
% Returns the figure and axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_cc_peak_map';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.addParameter('CCColumn', 'ccA', @(v) ischar(v) || isstring(v));
p.addParameter('LagColumn', 'lag_peak_frames', @(v) ischar(v) || isstring(v));
p.addParameter('HideOutliers', false, @(v) islogical(v) && isscalar(v));
p.addParameter('MarkerSize', 28, @(v) isnumeric(v) && isscalar(v) && v > 0);
p.addParameter('MaxSlices', Inf, @(v) isnumeric(v) && isscalar(v) && v > 0);
p.parse(varargin{:});
o = p.Results;
ccCol = char(o.CCColumn);
lagCol = char(o.LagColumn);

need = {'x_mm', 'z_mm', ccCol};
for k = 1:numel(need)
    if ~hasvar(P, need{k})
        error('viz:plot_cc_peak_map:missingColumn', ...
            'P must contain a %s column.', need{k});
    end
end

if hasvar(P, 'slice_idx')
    slices = unique(P.slice_idx, 'stable');
else
    slices = 1;
    P.slice_idx = ones(height(P), 1);
end
if isfinite(o.MaxSlices) && numel(slices) > o.MaxSlices
    idx = unique(round(linspace(1, numel(slices), o.MaxSlices)));
    slices = slices(idx);
end

nS = numel(slices);
nCol = ceil(sqrt(nS));
nRow = ceil(nS / nCol);
figW = min(1800, max(760, 180 * nCol));
figH = min(1200, max(520, 170 * nRow));
fig = figure('Color', 'w', 'Name', 'through-plane CC peak map', ...
    'Position', [80 80 figW figH]);
tl = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
ax = gobjects(1, nS);

ccAll = P.(ccCol);
cl = span(ccAll(isfinite(ccAll)));
xlimAll = span(P.x_mm(isfinite(P.x_mm)));
zlimAll = span(P.z_mm(isfinite(P.z_mm)));

for k = 1:nS
    ax(k) = nexttile(tl);
    Q = P(P.slice_idx == slices(k), :);
    keep = true(height(Q), 1);
    if o.HideOutliers && hasvar(Q, 'outlier')
        keep = keep & ~Q.outlier;
    end
    cc = Q.(ccCol);
    ok = keep & isfinite(cc);
    bad = keep & ~isfinite(cc);

    hold(ax(k), 'on');
    if any(bad)
        scatter(ax(k), Q.x_mm(bad), Q.z_mm(bad), 12, [0.72 0.72 0.72], ...
            'x', 'LineWidth', 0.6);
    end
    if any(ok)
        scatter(ax(k), Q.x_mm(ok), Q.z_mm(ok), o.MarkerSize, cc(ok), ...
            'filled', 'MarkerEdgeColor', 'none');
    end
    draw_vessel(ax(k), Q);
    hold(ax(k), 'off');

    axis(ax(k), 'equal');
    xlim(ax(k), xlimAll);
    ylim(ax(k), zlimAll);
    set(ax(k), 'YDir', 'reverse', 'Box', 'on', 'Layer', 'top', 'CLim', cl);
    title(ax(k), tile_title(Q, slices(k), ccCol, lagCol), ...
        'Interpreter', 'none', 'FontSize', 8);
    if k > (nRow - 1) * nCol
        xlabel(ax(k), 'x [mm]');
    end
    if mod(k - 1, nCol) == 0
        ylabel(ax(k), 'z [mm]');
    end
end

colormap(fig, parula(256));
lastAx = ax(find(isgraphics(ax), 1, 'last'));
cb = colorbar(lastAx);
cb.Label.String = sprintf('%s peak CC', ccCol);
title(tl, strtrim(sprintf('%s  ROI through-plane peak CC by slice', ...
    char(o.Title))), 'Interpreter', 'none');
end


function tf = hasvar(T, name)
tf = any(strcmp(T.Properties.VariableNames, name));
end


function cl = span(v)
v = v(isfinite(v));
if isempty(v)
    cl = [0 1];
    return;
end
lo = min(v);
hi = max(v);
if hi <= lo
    lo = lo - 0.5;
    hi = hi + 0.5;
end
cl = [lo hi];
end


function draw_vessel(ax, Q)
if ~hasvar(Q, 'r_mm') || ~any(isfinite(Q.r_mm))
    return;
end
R = max(Q.r_mm(isfinite(Q.r_mm)));
th = linspace(0, 2 * pi, 160);
plot(ax, R * cos(th), R * sin(th), '-', 'Color', [0.20 0.20 0.20], ...
    'LineWidth', 0.6);
end


function s = tile_title(Q, sliceValue, ccCol, lagCol)
pos = nan;
if hasvar(Q, 'slice_pos_mm') && height(Q) > 0
    pos = Q.slice_pos_mm(1);
end
cc = Q.(ccCol);
nOk = nnz(isfinite(cc));
ccMed = median(cc, 'omitnan');
lagMed = nan;
if hasvar(Q, lagCol)
    lagMed = median(Q.(lagCol), 'omitnan');
end
if isfinite(pos)
    head = sprintf('slice %g, y=%+.1f mm', sliceValue, pos);
else
    head = sprintf('slice %g', sliceValue);
end
if isfinite(lagMed)
    s = sprintf('%s | n=%d, cc=%.2f, lag=%.0f', head, nOk, ccMed, lagMed);
else
    s = sprintf('%s | n=%d, cc=%.2f', head, nOk, ccMed);
end
end
