function [vx, vy, vz, reach] = transit_truth(ph, x0, z0, y1, y2, d_row, nStep)
%TRANSIT_TRUTH Ground truth along the path the scatterers actually take.
%
%   [vx, vy, vz, reach] = src.transit_truth(ph, x0, z0, y1, y2, d_row)
%   [...] = src.transit_truth(..., nStep)
%
% x0, z0 are vessel-centred column vectors [m]; y1, y2 are the absolute
% elevations of the row1 and row2 planes [m]; d_row is their separation [m].
% vy is the through-plane speed the transit lag corresponds to [m/s], vx and
% vz the lateral velocities the stage A compensation should apply, and reach
% is false where no transit exists at all.
%
% What the estimator measures is a time: the lag at which the row1 patch
% reappears in row2. The matching truth is therefore d_row/T along the
% streamline, and T has to be integrated along that streamline, because the
% scatterers drift sideways while they cross. Sampling a fixed (x,z) column
% instead -- the same x and z at every y -- diverges from the real path
% wherever the flow is not parallel:
%
%   * in a shear layer the column stays at one radius while the streamline
%     moves to another, where uy is different (31% at r/R = 0.95 in the
%     throat, whose wall radius does not even change)
%   * where the vessel narrows the column eventually lies outside the wall
%     and reads the zero velocity there, while the real scatterer follows
%     the taper inwards and stays in the lumen (78.6 against 198.3 mm/s)
%   * downstream the column samples the recirculation and returns a negative
%     speed for a point whose scatterers never reach row2 at all
%
% vx and vz are the net lateral displacement over T rather than an average
% of the instantaneous values: they exist to reproduce the right row2
% offset, and that offset is a displacement. This is also why they are not
% the row1 lateral velocity that src.in_plane can measure -- in the throat
% the row1 value points inwards but reverses past the neck, so applying it
% over the whole transit overshoots the true offset five-fold.
%
% Every point is advanced in lockstep, so the cost is nStep phantom queries
% per slice, not per point. Truth does not depend on the speckle realisation
% either, so a multi-seed run pays this once per slice.
%
% Assumes flow along +y, which is how the phantoms and the transit scan are
% both oriented.
%
% See also SRC.ESTIMATE_SLICE, SRC.PHANTOM_CFD, SRC.PHANTOM_PARABOLIC.

if nargin < 7
    nStep = 200;
end
x0 = x0(:);
z0 = z0(:);
n = numel(x0);
dy = (y2 - y1) / nStep;

x = x0;
z = z0;
T = zeros(n, 1);
reach = true(n, 1);
for k = 0:nStep - 1
    [ux, uy, uz] = ph.velocity(x, repmat(y1 + k * dy, n, 1), z);
    stall = ~(uy > 0);            % stagnant, reversed, or off the solution
    reach = reach & ~stall;
    uy(stall) = nan;              % poisons T and the position from here on
    T = T + dy ./ uy;
    x = x + ux ./ uy * dy;
    z = z + uz ./ uy * dy;
end

vy = d_row ./ T;
vx = (x - x0) ./ T;
vz = (z - z0) ./ T;
vx(~reach) = nan;
vy(~reach) = nan;
vz(~reach) = nan;
end
