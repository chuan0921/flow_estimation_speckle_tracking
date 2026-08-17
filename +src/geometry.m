function [grid_x, grid_z, keep, GX, GZ, R_loc, ysp] = geometry(d, ph, opts)
% GEOMETRY  Lumen sampling grid, gated by how much of each ROI is vessel.
%
%   [grid_x, grid_z, keep, GX, GZ, R_loc, ysp] = src.geometry(d, ph, opts)
%
% ph is a phantom interface (src.phantom_cfd or src.phantom_parabolic); only
% ph.wall_radius is used here.
%
% An ROI is kept when its centre is inside the lumen and the echo energy
% falling outside the lumen is at most opts.roi_out_energy_max of the energy
% inside it. Both tests are on the row1 image against the row1 wall radius:
% row1 is where the ROI is defined, and the same scatterers reach row2 at a
% different position because the streamlines follow the taper -- locating
% them there is the compensated search's job, not a sampling constraint.
% The lattice therefore runs to the wall instead of stopping at a fixed
% fraction of the radius, and the window may overhang the wall even though
% the centre may not.
%
% R_loc is the narrowest wall radius along the row1-row2 transit path and is
% returned for callers that need the conservative radius.
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
R1 = ph.wall_radius(d.y_row1);
if ~isfinite(R1), R1 = R_loc; end

% Anchor the lattice on the vessel axis, not on -R: with -R:step:R the
% sequence only passes through 0 when R is an exact multiple of step, so a
% general R leaves no sample on the axis and shifts the whole lattice
% off-centre by up to half a step. The lattice reaches one ROI half-width
% past the widest wall; the energy gate below trims it.
k = ceil(R1 / opts.grid_step);
gv = (-k:k) * opts.grid_step;
[GX, GZ] = meshgrid(gv, gv);

keep = roi_energy_gate(d, opts, GX, GZ, R1);
grid_x = GX(keep) + d.xc;
grid_z = GZ(keep) + d.zc;
end


function keep = roi_energy_gate(d, opts, GX, GZ, R1)
%ROI_ENERGY_GATE Keep ROIs centred in the lumen whose out-of-lumen energy
% stays under the limit. Energy is the time-averaged envelope power, so the
% test reflects what the correlation actually sees rather than the geometric
% overlap alone.
E1 = mean(abs(d.iq_row1) .^ 2, 3);
[PX, PZ] = meshgrid(d.x, d.z);
in1 = hypot(PX - d.xc, PZ - d.zc) <= R1;

rz = round(opts.roi_half / d.dz);
rx = round(opts.roi_half / d.dx);
Nz = numel(d.z);
Nx = numel(d.x);
limit = opts.roi_out_energy_max;

keep = hypot(GX, GZ) <= R1 + 1e-12;    % the centre stays inside the wall
for p = find(keep).'
    ix = round((GX(p) + d.xc - d.x(1)) / d.dx) + 1;
    iz = round((GZ(p) + d.zc - d.z(1)) / d.dz) + 1;
    if ix - rx < 1 || ix + rx > Nx || iz - rz < 1 || iz + rz > Nz
        keep(p) = false;               % patch would run off the image
        continue;
    end
    keep(p) = passes(E1(iz - rz:iz + rz, ix - rx:ix + rx), ...
        in1(iz - rz:iz + rz, ix - rx:ix + rx), limit);
end
end


function tf = passes(Epatch, inside, limit)
e_in = sum(Epatch(inside));
if ~(e_in > 0)
    tf = false;
    return;
end
tf = sum(Epatch(~inside)) <= limit * e_in;
end
