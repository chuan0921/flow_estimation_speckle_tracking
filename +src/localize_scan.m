function L = localize_scan(snaps, d_row, opts)
%LOCALIZE_SCAN Probe direction and step from the tissue, per slice.
%
%   L = src.localize_scan(snaps, d_row, opts)
%
% The background tissue is static, so it works as the position reference:
% the row2 plane rides d_row ahead of row1 on the probe, and the tissue
% row2 sees at scan step k reappears in row1 once the probe has moved
% d_row. The convention is fixed and must not be changed piecemeal:
% ref is ALWAYS row2, moved is ALWAYS row1. A correlation peak at
% positive slice offset means the probe walks toward the row2 side; at
% negative offset, toward the row1 side. |peak| converts the known
% hardware spacing d_row into the probe step: step = d_row / |peak|.
%
% snaps is a struct array ordered by acquisition step, one element per
% slice, with fields E1, E2 (one envelope frame each), x, z, xc, zc,
% radius (lumen radius [m] at that slice). Pixels inside
% opts.loc_margin * radius are excluded -- blood moves, tissue must not.
% Each depth row of a frame is de-meaned first, so layered-tissue
% banding (y-invariant by construction) cannot prop up the baseline.
%
% Peaks are refined with a baseline-subtracted Gaussian three-point fit:
% the correlation peak has the elevation beam's Gaussian shape, and a
% plain parabola peak-locks toward the integer grid (-0.7% bias measured,
% Gaussian leaves -0.08%). Slices whose partner lies beyond the scan end
% fall back to the mirrored pairing (row1[k] against row2[j]) -- the same
% physical coincidence read from the other end, peak sign flipped back.
%
% opts.loc_search_mm > 0 switches the correlator to block matching: the
% same-position NCC above assumes the probe only translates along y, but
% probe pressure and an uneven phantom surface shift and compress the
% tissue in-plane between observations. In block mode opts.loc_blocks
% small tissue blocks (opts.loc_block_mm edge) are each searched over
% +/- loc_search_mm in x and z; per slice offset the block keeps its
% best in-window ZNCC, so an in-plane shift raises the true peak back
% instead of killing it. The slice-offset peak is taken per block and
% the median across blocks is reported (majority of blocks must agree).
% The same search window is used at every offset, so the elevated
% max-of-noise baseline cancels in the peak comparison. 0 (default)
% keeps the plain same-position path -- correct whenever the tissue
% does not deform, and the block mode's degenerate special case.
%
% L has one row per slice: direction (+1/-1), step_mm, cc_peak, pairing.
n = numel(snaps);
useBlocks = opts.loc_search_mm > 0;
if useBlocks
    px = mean(diff(snaps(1).x));
    hb = max(4, round(opts.loc_block_mm * 1e-3 / (2 * px)));
    sr = max(1, round(opts.loc_search_mm * 1e-3 / px));
    ctrs = pick_blocks(snaps, hb, sr, px, opts);
    if isempty(ctrs)
        warning('localize_scan:noBlocks', ...
            'no tissue block clears the lumen; same-position NCC used');
        useBlocks = false;
    else
        img1 = cell(1, n);
        img2 = cell(1, n);
        for k = 1:n
            img1{k} = deband(snaps(k).E1);
            img2{k} = deband(snaps(k).E2);
        end
    end
end
if ~useBlocks
    tpl1 = cell(1, n);
    tpl2 = cell(1, n);
    for k = 1:n
        s = snaps(k);
        [X, Z] = meshgrid(s.x, s.z);
        mask = hypot(X - s.xc, Z - s.zc) > opts.loc_margin * s.radius;
        tpl1{k} = norm_masked(deband(s.E1), mask);
        tpl2{k} = norm_masked(deband(s.E2), mask);
    end
end

ms = -opts.loc_range:opts.loc_range;
direction = zeros(n, 1);
step_mm = nan(n, 1);
cc_peak = nan(n, 1);
pairing = strings(n, 1);
for k = 1:n
    if useBlocks
        [mp, pk, ok] = block_peak(img2{k}, img1, k, ms, ctrs, hb, sr);
    else
        [mp, pk, ok] = find_peak(tpl2{k}, tpl1, k, ms);
    end
    pairing(k) = "forward";
    if ~ok
        % scan edge: the forward partner does not exist; read the same
        % coincidence from the other end and flip the sign back.
        if useBlocks
            [mp, pk, ok] = block_peak(img1{k}, img2, k, ms, ctrs, hb, sr);
        else
            [mp, pk, ok] = find_peak(tpl1{k}, tpl2, k, ms);
        end
        mp = -mp;
        pairing(k) = "mirrored";
    end
    if ok
        direction(k) = sign(mp);
        step_mm(k) = d_row * 1e3 / abs(mp);
        cc_peak(k) = pk;
    else
        pairing(k) = "failed";
    end
end
L = table((1:n)', direction, step_mm, cc_peak, pairing, ...
    'VariableNames', {'scan_step', 'direction', 'step_mm', 'cc_peak', ...
    'pairing'});
end


function [mp, pk, ok] = find_peak(ref, movers, k, ms)
%FIND_PEAK Best interior correlation peak of ref against the moved stack.
% |m| <= 1 is excluded: neighbouring slices correlate through the
% elevation beam overlap regardless of geometry, and m = 0 is the same
% position by definition.
n = numel(movers);
cc = nan(size(ms));
for i = 1:numel(ms)
    j = k + ms(i);
    if j >= 1 && j <= n
        cc(i) = movers{j}' * ref;
    end
end
[mp, pk, ok] = peak_from_curve(ms, cc);
end


function [mp, pk, ok] = block_peak(ref, movers, k, ms, ctrs, hb, sr)
%BLOCK_PEAK Slice-offset peak by block matching, median across blocks.
% Each block keeps its best in-window ZNCC per slice offset, so an
% in-plane tissue shift (probe pressure, uneven surface) moves the match
% inside the window instead of destroying it. Blocks vote: the majority
% must find an interior peak, and the median peak position wins -- one
% block sitting on bland tissue cannot steer the answer.
n = numel(movers);
nb = size(ctrs, 1);
mps = nan(nb, 1);
pks = nan(nb, 1);
for b = 1:nb
    iz = ctrs(b, 1); ix = ctrs(b, 2);
    tpl = ref(iz-hb:iz+hb, ix-hb:ix+hb);
    cc = nan(size(ms));
    for i = 1:numel(ms)
        j = k + ms(i);
        if j >= 1 && j <= n
            reg = movers{j}(iz-hb-sr:iz+hb+sr, ix-hb-sr:ix+hb+sr);
            cc(i) = max_zncc(tpl, reg);
        end
    end
    [mps(b), pks(b), okb, marg] = peak_from_curve(ms, cc);
    % Max-over-shifts lifts the noise baseline (max of ~441 unrelated
    % correlations), so an interior fluctuation can fake a peak where the
    % plain path would correctly fail -- and failing is load-bearing: it
    % is what hands the scan tail to the mirrored pairing. Gate on peak
    % significance: a genuine tissue re-sighting clears its own curve's
    % baseline by ~8 robust deviations, noise by ~1.
    if ~okb || marg < 4
        mps(b) = nan;
    end
end
sel = isfinite(mps);
ok = nnz(sel) > nb / 2;
mp = median(mps(sel));
pk = median(pks(sel));
if ~ok
    mp = 0; pk = nan;
end
end


function c = max_zncc(tpl, reg)
%MAX_ZNCC Best zero-normalised cross-correlation of tpl anywhere in reg.
% filter2 + box sums give the exact ZNCC at every full-overlap shift
% without toolbox dependencies.
t = double(tpl);
t = t - mean(t(:));
nt = norm(t(:));
r = double(reg);
num = filter2(t, r, 'valid');
box = ones(size(t));
s1 = filter2(box, r, 'valid');
s2 = filter2(box, r.^2, 'valid');
np = numel(t);
den = sqrt(max(s2 - s1.^2 / np, 0)) * nt;
c = max(num(:) ./ max(den(:), eps));
end


function C = pick_blocks(snaps, hb, sr, px, opts)
%PICK_BLOCKS Tissue block centres [iz ix], spread apart, clear of blood.
% Candidates keep the whole template block outside opts.loc_margin times
% the largest lumen radius of the scan (the block must stay in tissue at
% every slice it is matched against) and the search window inside the
% frame. The search window itself may graze the margin zone: a candidate
% shift that overlaps blood merely scores lower and loses the in-window
% max. Greedy max-min spread picks opts.loc_blocks of them so the
% blocks sample different tissue, not three copies of the same corner.
s = snaps(1);
[X, Z] = meshgrid(s.x, s.z);
rr = hypot(X - s.xc, Z - s.zc);
guard = opts.loc_margin * max([snaps.radius]) + sqrt(2) * hb * px;
m = hb + sr;
cand = false(size(rr));
cand(1+m:end-m, 1+m:end-m) = true;
cand = cand & rr > guard;
[iz, ix] = find(cand);
keep = mod(iz - 1, hb) == 0 & mod(ix - 1, hb) == 0;   % coarse grid
iz = iz(keep); ix = ix(keep);
if isempty(iz)
    C = zeros(0, 2);
    return;
end
nb = min(opts.loc_blocks, numel(iz));
d0 = rr(sub2ind(size(rr), iz, ix));
[~, first] = max(d0);
pick = first;
while numel(pick) < nb
    dmin = inf(numel(iz), 1);
    for p = pick(:)'
        dmin = min(dmin, hypot(iz - iz(p), ix - ix(p)));
    end
    dmin(pick) = -inf;
    [~, nxt] = max(dmin);
    pick(end+1) = nxt; %#ok<AGROW>
end
C = [iz(pick) ix(pick)];
end


function [mp, pk, ok, marg] = peak_from_curve(ms, cc)
%PEAK_FROM_CURVE Interior peak of a slice-offset curve, Gaussian sub-pixel.
% marg is the peak's significance: its height above the off-peak
% baseline in robust (MAD) deviations of the off-peak samples. The
% plain path ignores it; block mode gates on it.
sel = abs(ms) >= 2 & isfinite(cc);
mp = 0; pk = nan; ok = false; marg = 0;
if nnz(sel) < 3
    return;
end
ccs = cc;
ccs(~sel) = -inf;
[pk, ip] = max(ccs);
% interior peak with both neighbours, away from the search boundary
if ip <= 1 || ip >= numel(ms) || ~isfinite(cc(ip-1)) || ~isfinite(cc(ip+1))
    return;
end
base = median(cc(abs(ms) <= 2 & isfinite(cc)));
y = max([cc(ip-1) cc(ip) cc(ip+1)] - base, 1e-6);
g = log(y);
mp = ms(ip) + (g(1) - g(3)) / (2 * (g(1) - 2 * g(2) + g(3)));
ok = true;
off = sel;
off(max(ip-1, 1):min(ip+1, numel(ms))) = false;
o = cc(off);
if numel(o) >= 4
    spread = 1.4826 * median(abs(o - median(o)));
    marg = (pk - median(o)) / max(spread, 1e-6);
else
    marg = inf;
end
end


function E = deband(E)
E = double(E) - mean(double(E), 2);
end


function v = norm_masked(E, mask)
v = E(mask);
v = v - mean(v);
v = v / (norm(v) + eps);
end
