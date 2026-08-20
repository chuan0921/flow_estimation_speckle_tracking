function ph = phantom_cfd(cfd)
%PHANTOM_CFD Ground truth interpolated from a CFD solution.
%
%   ph = src.phantom_cfd(cfd)
%
% For datasets backed by a real CFD run, where the lumen may be stenosed
% (wall radius varies along y) and all three velocity components are
% non-trivial. Truth comes from trilinear interpolation of the solution.
%
% Required fields of cfd (grids in MILLIMETRES, velocities in m/s):
%   x_grid, y_grid, z_grid    meshgrid arrays over the solution domain,
%                             vessel-centred in x and z
%   grid_Ux, grid_Uy, grid_Uz velocity components on that grid
%   y_wall_unique             elevational stations of the wall table [mm]
%   r_wall_unique             lumen radius at those stations [mm]
%
% Returns the phantom interface consumed by src.estimate_slice() and src.geometry:
%   ph.name                    label for reporting
%   ph.wall_radius(y)          y [m] -> lumen radius [m], same size as y
%   ph.velocity(x, y, z)       vessel-centred x,z and absolute y, all [m]
%                              -> [ux, uy, uz] in m/s, same size as inputs
%
% Queries outside the solution domain return 0 velocity (interp3 extrap
% value), matching the previous in-line behaviour; queries outside the wall
% table return NaN radius, which src.geometry reports rather than ignores.
%
% See also SRC.PHANTOM_PARABOLIC.

need = {'x_grid', 'y_grid', 'z_grid', 'grid_Ux', 'grid_Uy', 'grid_Uz', ...
    'y_wall_unique', 'r_wall_unique'};
missing = need(~isfield(cfd, need));
if ~isempty(missing)
    error('src:phantom_cfd:missingFields', ...
        'cfd struct is missing: %s', strjoin(missing, ', '));
end

ph.name = 'cfd';
ph.compensate_xz = true;
ph.allow_reverse = true;
ph.wall_radius = @wall_radius;
ph.velocity = @velocity;

    function r = wall_radius(y)
        % Table is in mm; NaN outside it so the caller can complain.
        r = interp1(cfd.y_wall_unique, cfd.r_wall_unique, y * 1e3, ...
            'linear', nan) * 1e-3;
    end

    function [ux, uy, uz] = velocity(x, y, z)
        X = x * 1e3;
        Y = y * 1e3;
        Z = z * 1e3;
        ux = interp3(cfd.x_grid, cfd.y_grid, cfd.z_grid, cfd.grid_Ux, ...
            X, Y, Z, 'linear', 0);
        uy = interp3(cfd.x_grid, cfd.y_grid, cfd.z_grid, cfd.grid_Uy, ...
            X, Y, Z, 'linear', 0);
        uz = interp3(cfd.x_grid, cfd.y_grid, cfd.z_grid, cfd.grid_Uz, ...
            X, Y, Z, 'linear', 0);
    end
end
