function delta = parab(y1, y2, y3)
%PARAB Three-point parabolic sub-sample offset, clamped to half a sample.
den = double(y1) - 2 * double(y2) + double(y3);
delta = 0;
if abs(den) > eps
    delta = max(-0.5, min(0.5, 0.5 * (double(y1) - double(y3)) / den));
end
end
