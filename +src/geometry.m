function [grid_x, grid_z, keep, GX, GZ, R_loc, ysp] = geometry(d, ph, opts)
% GEOMETRY  Lumen sampling grid from the local wall radius.
%
%   [grid_x, grid_z, keep, GX, GZ, R_loc, ysp] = src.geometry(d, ph, opts)
%
% ph is a phantom interface (src.phantom_cfd or src.phantom_parabolic); only
% ph.wall_radius is used here.
%
% R_loc is the narrowest wall radius along the row1-row2 transit path, so a
% scatterer that stays inside the mask is inside the lumen for the whole
% transit. Points are kept inside 0.9*R_loc; the outer annulus is dropped
% because near-wall speckle is too weak to track.
%
% ysp is returned for the caller's path-average truth lookup.

ysp = linspace(d.y_row1, d.y_row2, 5);
r_path = ph.wall_radius(ysp);
if ~any(isfinite(r_path))
    error('src:geometry:noWallRadius', ...
        ['The phantom returned no finite wall radius over the transit ' ...
         'path y = %.3g..%.3g mm. Check that the wall table covers it.'], ...
        ysp(1) * 1e3, ysp(end) * 1e3);
end
R_loc = min(r_path, [], 'omitnan');

% Anchor the lattice on the vessel axis, not on -R: with -R:step:R the
% sequence only passes through 0 when R is an exact multiple of step, so a
% general R leaves no sample on the axis and shifts the whole lattice
% off-centre by up to half a step. Nothing is lost at the rim because keep
% drops everything outside 0.9*R_loc anyway.
k = floor(d.R / opts.grid_step);
gv = (-k:k) * opts.grid_step;
[GX, GZ] = meshgrid(gv, gv);
keep = (GX.^2 + GZ.^2) <= (0.9 * R_loc)^2;
grid_x = GX(keep) + d.xc;
grid_z = GZ(keep) + d.zc;
end
