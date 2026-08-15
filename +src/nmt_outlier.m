function out = nmt_outlier(F, thresh, eps_v)
%NMT_OUTLIER Normalized median test against the 8 grid neighbors.
[nr, nc] = size(F);
nb = nan(nr, nc, 8);
k = 0;
for dr = -1:1
    for dc = -1:1
        if dr == 0 && dc == 0
            continue;
        end
        k = k + 1;
        S = nan(nr, nc);
        rs = max(1, 1 + dr):min(nr, nr + dr);
        cs = max(1, 1 + dc):min(nc, nc + dc);
        S(rs - dr, cs - dc) = F(rs, cs);
        nb(:, :, k) = S;
    end
end
med = median(nb, 3, 'omitnan');
madv = median(abs(nb - med), 3, 'omitnan');
out = abs(F - med) ./ (1.4826 * madv + eps_v) > thresh;
out(isnan(F)) = false;
end
