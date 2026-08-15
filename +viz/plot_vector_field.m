function [fig, ax] = plot_vector_field(T, varargin)
%PLOT_VECTOR_FIELD Lumen cross-section: in-plane vectors on a vy colour map.
%
%   viz.plot_vector_field(T)
%   viz.plot_vector_field(T, 'Truth', true)      % estimate | CFD | error
%   viz.plot_vector_field(T, 'HideOutliers', false)
%
% T is the table returned by src.estimate_slice(). Marker colour is the through-plane
% component Vy; arrows are the in-plane (Vx,Vz) pair. Each panel scales its
% arrows to its own peak in-plane speed, printed under the axes, because the
% error field is typically an order of magnitude below the flow itself.
%
% Name-value options
%   Truth         false  add a CFD panel and an (estimate - CFD) panel
%   HideOutliers  true   drop rows flagged by the normalized median test
%   MarkerSize    60     scatter marker area [pt^2]
%   ArrowScale    1      multiplies the reference arrow length
%   Title         ''     prefix for the layout title
%
% Returns the figure handle and the axes handles.

p = inputParser;
p.FunctionName = 'viz.plot_vector_field';
p.addParameter('Truth', false, @(v) islogical(v) && isscalar(v));
p.addParameter('HideOutliers', true, @(v) islogical(v) && isscalar(v));
p.addParameter('MarkerSize', 60, @(v) isnumeric(v) && isscalar(v) && v > 0);
p.addParameter('ArrowScale', 1, @(v) isnumeric(v) && isscalar(v) && v > 0);
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.parse(varargin{:});
o = p.Results;

m = isfinite(T.Vy_mms) & isfinite(T.Vx_mms) & isfinite(T.Vz_mms);
if o.HideOutliers
    m = m & ~T.outlier;
end
if o.Truth
    m = m & isfinite(T.VyTrue_mms);
end
if ~any(m)
    error('viz:plot_vector_field:noPoints', ...
        'Nothing to plot: no finite%s rows in T.', ...
        subst(o.HideOutliers, ', non-outlier', ''));
end

x = T.x_mm(m);
z = T.z_mm(m);
vx = T.Vx_mms(m);
vz = T.Vz_mms(m);
vy = T.Vy_mms(m);

if o.Truth
    tx = T.VxTrue_mms(m);
    tz = T.VzTrue_mms(m);
    ty = T.VyTrue_mms(m);
    P = struct( ...
        'u', {vx, tx, vx - tx}, ...
        'v', {vz, tz, vz - tz}, ...
        'c', {vy, ty, vy - ty}, ...
        'name', {'estimate', 'CFD truth', 'error (est - CFD)'}, ...
        'diverge', {false, false, true});
    shared = [vy; ty];
else
    P = struct('u', vx, 'v', vz, 'c', vy, 'name', 'estimate', ...
        'diverge', false);
    shared = vy;
end

% Panels showing a velocity share one colour scale; an error panel gets its
% own symmetric one so the sign of the bias stays readable.
cl_shared = span(shared);
gs = grid_spacing(x, z);

fig = figure('Color', 'w', 'Name', 'flow_estimation vector field');
tl = tiledlayout(fig, 1, numel(P), 'TileSpacing', 'compact', ...
    'Padding', 'compact');
ax = gobjects(1, numel(P));

for k = 1:numel(P)
    ax(k) = nexttile(tl);
    if P(k).diverge
        L = max(abs(P(k).c));
        cl = [-max(L, eps), max(L, eps)];
    else
        cl = cl_shared;
    end

    speed = hypot(P(k).u, P(k).v);
    vmax = max(speed);
    len = o.ArrowScale * gs / max(vmax, eps);    % mm of arrow per mm/s

    hold(ax(k), 'on');
    scatter(ax(k), x, z, o.MarkerSize, P(k).c, 'filled', ...
        'MarkerEdgeColor', 'none');
    quiver(ax(k), x, z, len * P(k).u, len * P(k).v, 0, ...
        'Color', [0.15 0.15 0.15], 'LineWidth', 0.6, 'MaxHeadSize', 0.4);
    hold(ax(k), 'off');

    axis(ax(k), 'equal', 'tight');
    set(ax(k), 'CLim', cl, 'Box', 'on', 'Layer', 'top', 'YDir', 'reverse');
    colormap(ax(k), colours(P(k).diverge));
    title(ax(k), P(k).name);
    xlabel(ax(k), sprintf('x [mm]   (arrow = %.2g mm/s)', vmax));
    if k == 1
        ylabel(ax(k), 'z [mm]  (depth)');
    end
    cb = colorbar(ax(k));
    cb.Label.String = 'v_y [mm/s]';
end

title(tl, strtrim(sprintf('%s  %d points%s', char(o.Title), nnz(m), ...
    subst(o.HideOutliers, ', outliers hidden', ''))), 'Interpreter', 'none');
end


function s = subst(tf, a, b)
if tf, s = a; else, s = b; end
end


function cl = span(v)
%SPAN Finite data range, widened to a non-degenerate interval.
lo = min(v);
hi = max(v);
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    lo = lo - 1;
    hi = hi + 1;
end
cl = [lo, hi];
end


function gs = grid_spacing(x, z)
%GRID_SPACING Median lattice pitch [mm], used as the reference arrow length.
dx = median(diff(unique(round(x, 6))));
dz = median(diff(unique(round(z, 6))));
gs = min([dx, dz], [], 'omitnan');
if ~isfinite(gs) || gs <= 0
    gs = 0.25;      % default opts.grid_step, in mm
end
end


function C = colours(diverge)
%COLOURS Sequential map for velocities, blue-white-red for signed error.
if ~diverge
    C = parula(256);
    return;
end
n = 128;
t = linspace(0, 1, n)';
lower = [t, t, ones(n, 1)];                 % blue  -> white
upper = [ones(n, 1), flipud(t), flipud(t)]; % white -> red
C = [lower; upper];
end
