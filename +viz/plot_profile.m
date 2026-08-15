function [fig, ax] = plot_profile(prof, varargin)
%PLOT_PROFILE Blood-flow velocity along one vessel-centre diameter.
%
%   viz.plot_profile(prof)
%   viz.plot_profile(prof, 'Title', 'sim_stenosis50_80cms')
%
% prof is one slice from profile.csv. Its signed diameter_pos_mm coordinate
% runs from the left wall (-R), through the centre (0), to the right wall
% (+R). The wall endpoints are exact no-slip constraints.

p = inputParser;
p.FunctionName = 'viz.plot_profile';
p.addParameter('Title', '', @(v) ischar(v) || isstring(v));
p.parse(varargin{:});
o = p.Results;

if ~ismember('diameter_pos_mm', prof.Properties.VariableNames)
    error('viz:plot_profile:legacyRadialProfile', ...
        'profile.csv is radial. Re-run main to create a diameter BFVP.');
end
if numel(unique(prof.slice_idx)) ~= 1
    error('viz:plot_profile:oneSlice', ...
        'A BFVP is slice-specific; pass exactly one slice.');
end
if ismember('seed_base', prof.Properties.VariableNames)
    prof = prof(prof.seed_base == min(prof.seed_base), :);
end
prof = sortrows(prof, 'diameter_pos_mm');
x = prof.diameter_pos_mm;

fig = figure('Color', 'w', 'Name', 'diameter BFVP', ...
    'Position', [100 100 900 620]);
ax = axes(fig);
hold(ax, 'on');
xline(ax, 0, ':', 'Vessel centre', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom');

hRaw = scatter(ax, x, prof.vy_raw, 58, [0.20 0.45 0.75], 'filled', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
hFill = plot(ax, x, prof.vy_fill, '--', 'Color', [0.20 0.45 0.75], ...
    'LineWidth', 1.8);
hSg = plot(ax, x, prof.vy_sg, '-', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 2.6);
hTruth = plot(ax, x, prof.vy_true, '--k', 'LineWidth', 2.2);

radius = max(abs(x), [], 'omitnan');
xlim(ax, [-radius radius]);
hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontSize', 14);
xlabel(ax, 'Position along centre diameter x [mm]', 'FontSize', 18);
ylabel(ax, 'Through-plane velocity v_y [mm/s]', 'FontSize', 18);
legend(ax, [hRaw hFill hSg hTruth], ...
    {'Raw ROI estimate', 'No-slip interpolation', ...
    'Estimate after 2-D SG', 'CFD truth'}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', 13);

region = '';
if ismember('vessel_region', prof.Properties.VariableNames)
    region = char(prof.vessel_region(1));
end
if isempty(region)
    secondLine = 'Diameter BFVP';
else
    secondLine = sprintf('%s | Diameter BFVP', region);
end
ttl = sprintf('%s | Slice %d (y = %.1f mm)\n%s', ...
    char(o.Title), prof.slice_idx(1), prof.slice_pos_mm(1), secondLine);
title(ax, strtrim(ttl), 'Interpreter', 'none', 'FontSize', 18);
end
