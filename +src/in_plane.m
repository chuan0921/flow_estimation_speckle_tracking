function [vx, vz, n_used, samples] = in_plane(E1, IQ1, iz0, ix0, ...
    rz, rx, Nt, lag_ip, cap, dt, d, o, Nz, Nx)
%IN_PLANE Chained multi-lag block matching on row1: short lags
% predict the search center of longer lags; vz gets a carrier-phase
% refinement; origin-constrained weighted slope gives the velocities.
% Local functions: block_match, robust_slope.
vx = nan;
vz = nan;
tau = []; dxm = []; dzm = []; w = [];
cc_ref = nan;
vhx = 0;
vhz = 0;
for il = 1:numel(lag_ip)
    L = lag_ip(il);
    if L > cap
        continue;
    end
    ctr_x = round(vhx * L * dt / d.dx);
    ctr_z = round(vhz * L * dt / d.dz);
    if ix0 - rx + ctr_x - o.inplane_win < 1 || ...
            ix0 + rx + ctr_x + o.inplane_win > Nx || ...
            iz0 - rz + ctr_z - o.inplane_win < 1 || ...
            iz0 + rz + ctr_z + o.inplane_win > Nz
        continue;
    end
    gate = o.cc_abs_min;
    if ~isnan(cc_ref)
        gate = max(o.cc_abs_min, o.cc_rel_min * cc_ref);
    end
    [dxp, dzp, pk, oz_pk, ox_pk, tt] = block_match(E1, ix0, iz0, rx, rz, ...
        L, ctr_x, ctr_z, o.inplane_win, Nt, gate, o.inplane_pairs);
    if isnan(dxp)
        continue;
    end
    if isnan(cc_ref)
        cc_ref = pk;
    end
    Ai = double(reshape(IQ1(iz0 - rz:iz0 + rz, ix0 - rx:ix0 + rx, tt), ...
        [], numel(tt)));
    Bi = double(reshape(IQ1(iz0 - rz + oz_pk:iz0 + rz + oz_pk, ...
        ix0 - rx + ox_pk:ix0 + rx + ox_pk, tt + L), [], numel(tt)));
    C = mean(sum(conj(Ai) .* Bi, 1) ./ ...
        (vecnorm(Ai, 2, 1) .* vecnorm(Bi, 2, 1)));
    dz_res = -angle(C) / (2 * d.k0);
    if abs(dz_res) < d.lambda / 4
        dzp = oz_pk + dz_res / d.dz;
    end
    tau(end + 1) = L * dt; %#ok<AGROW>
    dxm(end + 1) = dxp * d.dx; %#ok<AGROW>
    dzm(end + 1) = dzp * d.dz; %#ok<AGROW>
    w(end + 1) = pk; %#ok<AGROW>
    vhx = sum(w .* tau .* dxm) / sum(w .* tau.^2);
    vhz = sum(w .* tau .* dzm) / sum(w .* tau.^2);
end
n_used = numel(tau);
if n_used >= 2
    vx = robust_slope(tau, dxm, w);
    vz = robust_slope(tau, dzm, w);
elseif n_used == 1
    vx = dxm / tau;
    vz = dzm / tau;
end
samples = struct('tau', tau, 'dx', dxm, 'dz', dzm, 'w', w);
end


function [dx_px, dz_px, pk, oz_pk, ox_pk, tt, Rmap] = block_match(V, ix0, ...
    iz0, rx, rz, L, ctr_x, ctr_z, win, Nt, cc_gate, n_pairs)
%BLOCK_MATCH Ensemble ZNCC block matching around a predicted shift.
% Returns the subpixel displacement, peak value, integer peak shift, the
% reference-frame list used for the ensemble and the correlation map.
dx_px = nan; dz_px = nan; pk = nan; oz_pk = nan; ox_pk = nan;
step = max(1, floor((Nt - L) / n_pairs));
tt = 1:step:(Nt - L);
A = double(reshape(V(iz0 - rz:iz0 + rz, ix0 - rx:ix0 + rx, tt), [], numel(tt)));
A = A - mean(A, 1);
A = A ./ vecnorm(A, 2, 1);
Rmap = nan(2 * win + 1, 2 * win + 1);
for jz = -win:win
    for jx = -win:win
        oz = ctr_z + jz;
        ox = ctr_x + jx;
        B = double(reshape(V(iz0 - rz + oz:iz0 + rz + oz, ...
            ix0 - rx + ox:ix0 + rx + ox, tt + L), [], numel(tt)));
        B = B - mean(B, 1);
        B = B ./ vecnorm(B, 2, 1);
        Rmap(jz + win + 1, jx + win + 1) = mean(sum(A .* B, 1));
    end
end
[pval, imax] = max(Rmap(:));
[pz, px] = ind2sub(size(Rmap), imax);
if pval < cc_gate || pz == 1 || pz == 2 * win + 1 || px == 1 || px == 2 * win + 1
    return;
end
pk = pval;
oz_pk = ctr_z + pz - (win + 1);
ox_pk = ctr_x + px - (win + 1);
dz_px = oz_pk + src.parab(Rmap(pz - 1, px), Rmap(pz, px), Rmap(pz + 1, px));
dx_px = ox_pk + src.parab(Rmap(pz, px - 1), Rmap(pz, px), Rmap(pz, px + 1));
end


function v = robust_slope(tau, disp_m, w)
%ROBUST_SLOPE Origin-constrained weighted LS slope with one trim pass
% (drop measurements whose displacement residual exceeds half a pixel).
v = sum(w .* tau .* disp_m) / sum(w .* tau.^2);
res = abs(disp_m - v * tau);
keep = res <= 2.5e-5;
if nnz(keep) >= 2 && any(~keep)
    v = sum(w(keep) .* tau(keep) .* disp_m(keep)) / ...
        sum(w(keep) .* tau(keep).^2);
end
end
