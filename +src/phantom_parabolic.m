function ph = phantom_parabolic(d, varargin)
%PHANTOM_PARABOLIC Analytic straight-pipe ground truth (Hagen-Poiseuille).
%
%   ph = src.phantom_parabolic(d)
%   ph = src.phantom_parabolic(d, 'Radius', 2.5e-3, 'Vmax', 0.8)
%
% For pure-simulation datasets (psf_sim) where the scatterer field is a
% straight pipe with a parabolic profile, so no external CFD solution
% exists. Truth is evaluated in closed form:
%
%     uy(r) = +/- vmax * (1 - (r/R)^2),   ux = uz = 0,   r = hypot(x,z)
%
% The sign follows d.flow_direction ('row2_to_row1' reverses it). The wall
% radius is constant along y, so the transit path never leaves the lumen.
%
% d is the slice struct loaded by src.estimate_slice(); R defaults to d.R and vmax to
% d.vmax. Both can be overridden for phantoms that do not carry them.
%
% Returns the phantom interface consumed by src.estimate_slice() and src.geometry:
%   ph.name                    label for reporting
%   ph.wall_radius(y)          y [m] -> lumen radius [m], same size as y
%   ph.velocity(x, y, z)       vessel-centred x,z and absolute y, all [m]
%                              -> [ux, uy, uz] in m/s, same size as inputs
%
% See also SRC.PHANTOM_CFD.

p = inputParser;
p.FunctionName = 'src.phantom_parabolic';
p.addParameter('Radius', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v > 0));
p.addParameter('Vmax', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
p.parse(varargin{:});
o = p.Results;

R = pick(o.Radius, d, 'R', 'src.phantom_parabolic');
vmax = pick(o.Vmax, d, 'vmax', 'src.phantom_parabolic');

sgn = 1;
if isfield(d, 'flow_direction') && ...
        strcmpi(strtrim(char(d.flow_direction)), 'row2_to_row1')
    sgn = -1;
end

ph.name = sprintf('parabolic (R = %.2f mm, vmax = %+.3g m/s)', R * 1e3, sgn * vmax);
ph.wall_radius = @(y) R * ones(size(y));
ph.velocity = @velocity;

    function [ux, uy, uz] = velocity(x, ~, z)
        % Uniform along y, so the elevational coordinate is unused.
        rr = hypot(x, z);
        uy = sgn * vmax * max(0, 1 - (rr / R).^2);
        ux = zeros(size(uy));
        uz = zeros(size(uy));
    end
end


function v = pick(override, d, field, caller)
%PICK Explicit override, else the slice field, else a clear error.
if ~isempty(override)
    v = override;
    return;
end
if ~isfield(d, field) || isempty(d.(field))
    error('src:phantom_parabolic:missingField', ...
        ['%s needs ''%s'' in the slice file, or pass it explicitly ' ...
         'as a name-value option.'], caller, field);
end
v = d.(field);
end
