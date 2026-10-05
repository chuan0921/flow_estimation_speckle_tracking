function fn = through_plane()
%THROUGH_PLANE Through-plane (vy) primitives, grouped in one file.
% MATLAB exposes only a file's first function, so the primitives are
% handed back as handles:
%     tp = src.through_plane();
%     cc = tp.scan_transit_lags(...);
% prominence_ok and robust_z are internal helpers (no handle).
fn.normalize_patch   = @normalize_patch;
fn.scan_transit_lags = @scan_transit_lags;
fn.scan_transit_lags_spatial = @scan_transit_lags_spatial;
fn.pick_peak         = @pick_peak;
fn.pick_significant  = @pick_significant;
fn.drift_match       = @drift_match;
fn.refine_drift      = @refine_drift;
end


function An = normalize_patch(E, iz0, ix0, rz, rx, Nt)
%NORMALIZE_PATCH Zero-mean unit-norm frames of one ROI, single precision.
An = reshape(E(iz0 - rz:iz0 + rz, ix0 - rx:ix0 + rx, :), [], Nt);
An = An - mean(An, 1);
n = vecnorm(An, 2, 1);
n(n == 0) = inf;
An = An ./ n;
end

function [cc, cca, ccb] = scan_transit_lags(An, E2, iz0, ix0, rz, rx, ...
    Nt, lags, oz_list, ox_list, Nz, Nx, fft_min_lags)
%SCAN_TRANSIT_LAGS Full-pair ensemble ZNCC per signed lag with a
% lag-dependent integer xz shift of the row2 patch. Lags sharing the same
% integer shift reuse one normalized patch. Groups with many lags use the
% FFT cross-correlation identity (numerically equivalent to the direct
% sliced dot products); small groups stay on the direct path.
% Optional outputs cca/ccb are the first/second-half pair means (for
% reproducibility gating); requesting them forces the direct path.
% Complex (coherent) input correlates conj(An).*Bn and reports the modulus
% of the coherent pair mean; real (envelope) input keeps the signed mean,
% bit-identical to the historical behaviour.
want_halves = nargout > 1;
coherent = ~isreal(An);
cc = nan(size(lags));
cca = nan(size(lags));
ccb = nan(size(lags));
bz = round(oz_list(:));
bx = round(ox_list(:));
G = unique([bz, bx], 'rows');
M = 2^nextpow2(2 * Nt);
AnF = [];
for g = 1:size(G, 1)
    sz = G(g, 1);
    sx = G(g, 2);
    if iz0 - rz + sz < 1 || iz0 + rz + sz > Nz || ...
            ix0 - rx + sx < 1 || ix0 + rx + sx > Nx
        continue;
    end
    idx = find(bz == sz & bx == sx)';
    Bn = normalize_patch(E2, iz0 + sz, ix0 + sx, rz, rx, Nt);
    if numel(idx) >= fft_min_lags && ~want_halves
        if isempty(AnF)
            AnF = conj(fft(An, M, 2));
        end
        S = sum(AnF .* fft(Bn, M, 2), 1);
        if coherent
            r = ifft(S, [], 2);            % r(L+1) = coherent pair sum at L
        else
            r = ifft(S, [], 2, 'symmetric'); % r(L+1) = pair sum at lag L
        end
        for il = idx
            L = lags(il);
            cc(il) = abs_if(coherent, r(L + 1)) / (Nt - L);
        end
    else
        for il = idx
            L = lags(il);
            if L >= 0
                c = sum(conj(An(:, 1:Nt - L)) .* Bn(:, 1 + L:Nt), 1);
            else
                c = sum(conj(An(:, 1 - L:Nt)) .* Bn(:, 1:Nt + L), 1);
            end
            cc(il) = abs_if(coherent, mean(c));
            if want_halves
                h = floor(numel(c) / 2);
                cca(il) = abs_if(coherent, mean(c(1:h)));
                ccb(il) = abs_if(coherent, mean(c(h + 1:end)));
            end
        end
    end
end
cc = double(cc);
cca = double(cca);
ccb = double(ccb);
end

function v = abs_if(coherent, v)
%ABS_IF Modulus for the coherent path, signed passthrough for envelopes.
if coherent
    v = abs(v);
end
end

function cc = scan_transit_lags_spatial(E1, E2, iz0, ix0, rz, rx, ...
    Nt, lags, oz_list, ox_list, Nz, Nx, fft_min_lags, spatial_half)
%SCAN_TRANSIT_LAGS_SPATIAL Speckle-size transit scan with spatial filtering.
% Each lateral sub-window is normalized locally, scanned independently, and
% the CC-lag functions are spatially weighted before peak picking. This keeps
% the vy measurement local after xz compensation while reducing peak hopping.
if nargin < 14 || isempty(spatial_half)
    spatial_half = 1;
end
rz = max(1, rz);
rx = max(1, rx);
step_x = max(1, 2 * rx + 1);
idx = -spatial_half:spatial_half;
subs = idx * step_x;
w = spatial_half + 1 - abs(idx);
C = nan(numel(subs), numel(lags));
for m = 1:numel(subs)
    sx = subs(m);
    if iz0 - rz < 1 || iz0 + rz > Nz || ...
            ix0 + sx - rx < 1 || ix0 + sx + rx > Nx
        continue;
    end
    An = normalize_patch(E1, iz0, ix0 + sx, rz, rx, Nt);
    C(m, :) = scan_transit_lags(An, E2, iz0, ix0 + sx, rz, rx, Nt, ...
        lags, oz_list, ox_list, Nz, Nx, fft_min_lags);
end
fin = isfinite(C);
W = repmat(w(:), 1, numel(lags));
W(~fin) = 0;
C(~fin) = 0;
wsum = sum(W, 1);
cc = sum(W .* C, 1) ./ wsum;
cc(wsum == 0) = nan;
end

function ok = prominence_ok(cc, j, prom_frac, prom_w)
%PROMINENCE_OK True transit peaks are local bumps that decay on both
% sides; frozen-speckle plateaus do not. Require the curve to drop by
% prom_frac of the peak height (above the curve median) within prom_w
% samples on each side. A peak too close to the scan edge to test is
% rejected (plateaus ride the edges).
ok = false;
n = numel(cc);
if j - prom_w < 1 || j + prom_w > n
    return;
end
fin = isfinite(cc);
med = median(cc(fin));
need = cc(j) - prom_frac * (cc(j) - med);
left = min(cc(j - prom_w:j - 1));
right = min(cc(j + 1:j + prom_w));
ok = left <= need && right <= need;
end

function [vy, pk, shape] = pick_peak(cc, lags, d_row, dt, cc_min, ...
    prom_frac, prom_w, prom_tail_min, prom_abs_min)
%PICK_PEAK Global peak with cc_min gate, interior only, local-step
% parabolic subframe interpolation (valid on non-uniform lag axes).
% Optional prominence gate for long-lag peaks only (frozen-speckle
% plateaus live in the tail; genuine broad shear peaks in the mid axis
% are spared).
%
% shape carries two peak diagnostics, filled whenever an interior maximum
% exists -- including for peaks the gates reject, so a threshold can be
% swept after the fact instead of by re-running:
%   best   the raw interior maximum (pk only reports gate-passing peaks)
%   curv   -d2(cc)/d(lag)2 at the peak [1/frame^2]. Time-delay theory puts
%          the estimator's random error at var ~ 1/(SNR * curvature): a
%          broad shear-smeared peak localizes poorly no matter how high its
%          cc value is, so this predicts jitter where cc alone cannot.
% A peak-skew diagnostic (asymmetry) was also recorded for a while to test
% whether the elevation-decorrelation envelope drags the fitted vertex
% toward short lags; measured skew was zero-centred and uncorrelated with
% the bias, so the column was dropped and the bias traced to the velocity
% mixture inside the window instead.
vy = nan;
pk = nan;
shape = struct('best', nan, 'curv', nan);
[pval, j] = max(cc);
if isnan(pval)
    return;
end
shape.best = pval;
if j > 1 && j < numel(lags)
    step = (lags(j + 1) - lags(j - 1)) / 2;
    shape.curv = -(cc(j - 1) - 2 * cc(j) + cc(j + 1)) / step^2;
end
if pval < cc_min || j == 1 || j == numel(lags)
    return;
end
if nargin >= 8 && abs(lags(j)) >= prom_tail_min && ...
        ~prominence_ok(cc, j, prom_frac, prom_w)
    return;
end
% Absolute prominence: a flat curve's micro-ripple passes the relative
% test above, but its peak barely clears the curve median.
if nargin >= 9 && prom_abs_min > 0 && ...
        pval - median(cc(isfinite(cc))) < prom_abs_min
    return;
end
den = cc(j - 1) - 2 * cc(j) + cc(j + 1);
delta = 0;
if isfinite(den) && abs(den) > eps
    delta = max(-0.5, min(0.5, 0.5 * (cc(j - 1) - cc(j + 1)) / den));
end
step = (lags(j + 1) - lags(j - 1)) / 2;
vy = d_row / ((lags(j) + delta * step) * dt);
pk = pval;
end

function [vy, pk, zpk] = pick_significant(cc, cca, ccb, lags, d_row, dt, o)
%PICK_SIGNIFICANT Significance-gated peak for weak-signal lag curves.
% Gates (validated on the reverse-flow experiments): robust z-score
% against the curve's own median/MAD background, neighborhood breadth,
% first/second-half reproducibility, interior peak only. Candidates are
% tried in descending z order; the first one passing all gates wins.
vy = nan;
pk = nan;
zpk = nan;
fin = isfinite(cc);
if nnz(fin) < 5
    return;
end
z = robust_z(cc);
za = robust_z(cca);
zb = robust_z(ccb);
idx = find(fin);
[~, order] = sort(z(idx), 'descend');
for k = order(:)'
    j = idx(k);
    if j == 1 || j == numel(lags)
        continue;
    end
    nb = max(1, j - o.nb_half):min(numel(lags), j + o.nb_half);
    if z(j) < o.z_sig_rev || mean(z(nb), 'omitnan') < o.z_broad || ...
            za(j) < o.z_half || zb(j) < o.z_half || ...
            ~prominence_ok(cc, j, o.prom_frac, o.prom_w)
        continue;
    end
    den = cc(j - 1) - 2 * cc(j) + cc(j + 1);
    delta = 0;
    if isfinite(den) && abs(den) > eps
        delta = max(-0.5, min(0.5, 0.5 * (cc(j - 1) - cc(j + 1)) / den));
    end
    step = (lags(j + 1) - lags(j - 1)) / 2;
    vy = d_row / ((lags(j) + delta * step) * dt);
    pk = cc(j);
    zpk = z(j);
    return;
end
end

function z = robust_z(cc)
fin = isfinite(cc);
med = median(cc(fin));
mad = median(abs(cc(fin) - med));
z = (cc - med) / (1.4826 * mad + eps);
end

function [ox_tot, oz_tot, pk, Rmap, cz, cx] = drift_match(An, E2, iz0, ...
    ix0, rz, rx, L, ctr_z, ctr_x, win, Nt, Nz, Nx)
%DRIFT_MATCH 2D row1-row2 drift search at the transit lag around the
% predicted shift. Returns the accumulated drift in px. Optional extra
% outputs: the correlation map and its integer anchor (for plotting).
ox_tot = nan; oz_tot = nan; pk = nan;
Rmap = [];
cz = round(ctr_z);
cx = round(ctr_x);
if iz0 - rz + cz - win < 1 || iz0 + rz + cz + win > Nz || ...
        ix0 - rx + cx - win < 1 || ix0 + rx + cx + win > Nx
    return;
end
step = max(1, floor((Nt - L) / 40));
tt = 1:step:(Nt - L);
Ap = An(:, tt);
Rmap = nan(2 * win + 1, 2 * win + 1);
for jz = -win:win
    for jx = -win:win
        B = reshape(E2(iz0 - rz + cz + jz:iz0 + rz + cz + jz, ...
            ix0 - rx + cx + jx:ix0 + rx + cx + jx, tt + L), [], numel(tt));
        B = B - mean(B, 1);
        n = vecnorm(B, 2, 1);
        n(n == 0) = inf;
        B = B ./ n;
        Rmap(jz + win + 1, jx + win + 1) = mean(sum(Ap .* B, 1));
    end
end
[pval, imax] = max(Rmap(:));
[pz, px] = ind2sub(size(Rmap), imax);
if pval < 0.1 || pz == 1 || pz == 2 * win + 1 || px == 1 || px == 2 * win + 1
    return;
end
pk = double(pval);
oz_tot = cz + pz - (win + 1) + src.parab(Rmap(pz - 1, px), Rmap(pz, px), Rmap(pz + 1, px));
ox_tot = cx + px - (win + 1) + src.parab(Rmap(pz, px - 1), Rmap(pz, px), Rmap(pz, px + 1));
end

function [ox_f, oz_f, Rmap] = refine_drift(E1, E2, IQ1, IQ2, iz0, ix0, ...
    rzf, rxf, L, cz, cx, win, Nt, Nz, Nx, k0, lambda, dz_m)
%REFINE_DRIFT Fine drift stage: a consultant panel of three one-grain
% sub-windows tiling the two-grain footprint (locally-normalized spatial
% correlation filter -- small windows see less shear, their ensemble ZNCC
% maps are averaged with equal weight). x subpixel from the panel envelope
% parabola; z subpixel from the carrier phase of the complex cross-row
% correlation averaged over the sub-windows at the integer peak (envelope
% parabola fallback when the phase residual exceeds lambda/4). NaN when a
% window runs off the image or the peak sits on the search edge. Optional
% output: the panel correlation map.
ox_f = nan; oz_f = nan; Rmap = [];
rxf1 = max(1, round(rxf / 2));
subs = [-rxf1, 0, rxf1];
if iz0 - rzf + min(cz - win, 0) < 1 || ...
        iz0 + rzf + max(cz + win, 0) > Nz || ...
        ix0 - rxf + min(cx - win, 0) < 1 || ...
        ix0 + rxf + max(cx + win, 0) > Nx
    return;
end
step = max(1, floor((Nt - L) / 40));
tt = 1:step:(Nt - L);
nS = 2 * win + 1;

% Per-sub normalized row1 reference patches.
As = cell(1, 3);
for m = 1:3
    s = subs(m);
    A = double(reshape(E1(iz0 - rzf:iz0 + rzf, ...
        ix0 + s - rxf1:ix0 + s + rxf1, tt), [], numel(tt)));
    A = A - mean(A, 1);
    n = vecnorm(A, 2, 1); n(n == 0) = inf;
    As{m} = A ./ n;
end

Rmap = zeros(nS, nS);
for jz = -win:win
    for jx = -win:win
        acc = 0;
        for m = 1:3
            s = subs(m);
            B = double(reshape(E2(iz0 - rzf + cz + jz:iz0 + rzf + cz + jz, ...
                ix0 + s - rxf1 + cx + jx:ix0 + s + rxf1 + cx + jx, ...
                tt + L), [], numel(tt)));
            B = B - mean(B, 1);
            n = vecnorm(B, 2, 1); n(n == 0) = inf;
            acc = acc + mean(sum(As{m} .* (B ./ n), 1));
        end
        Rmap(jz + win + 1, jx + win + 1) = acc / 3;
    end
end

[pv, im] = max(Rmap(:));
[pz, px] = ind2sub(size(Rmap), im);
if pv < 0.1 || pz == 1 || pz == nS || px == 1 || px == nS
    return;
end
z_int = cz + pz - (win + 1);
x_int = cx + px - (win + 1);
ox_f = x_int + src.parab(Rmap(pz, px - 1), Rmap(pz, px), Rmap(pz, px + 1));
oz_f = z_int + src.parab(Rmap(pz - 1, px), Rmap(pz, px), Rmap(pz + 1, px));

% Carrier-phase z: complex cross-row correlation averaged over sub-windows.
C = 0;
for m = 1:3
    s = subs(m);
    Ai = double(reshape(IQ1(iz0 - rzf:iz0 + rzf, ...
        ix0 + s - rxf1:ix0 + s + rxf1, tt), [], numel(tt)));
    Bi = double(reshape(IQ2(iz0 - rzf + z_int:iz0 + rzf + z_int, ...
        ix0 + s - rxf1 + x_int:ix0 + s + rxf1 + x_int, tt + L), ...
        [], numel(tt)));
    C = C + mean(sum(conj(Ai) .* Bi, 1) ./ ...
        (vecnorm(Ai, 2, 1) .* vecnorm(Bi, 2, 1)));
end
dz_res = -angle(C / 3) / (2 * k0);
if abs(dz_res) < lambda / 4
    oz_f = z_int + dz_res / dz_m;
end
end
