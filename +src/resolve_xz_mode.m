function tf = resolve_xz_mode(value, ph)
%RESOLVE_XZ_MODE Resolve auto/on/off XZ compensation for a phantom.
if nargin < 1 || isempty(value)
    value = 'auto';
end
if islogical(value) && isscalar(value)
    tf = value;
    return;
end
if isnumeric(value) && isscalar(value) && isfinite(value) && ...
        (value == 0 || value == 1)
    tf = logical(value);
    return;
end
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('src:resolve_xz_mode:badValue', ...
        'compensate_xz must be auto, on, off, true, or false.');
end

switch lower(strtrim(char(value)))
    case 'on'
        tf = true;
    case 'off'
        tf = false;
    case 'auto'
        if nargin >= 2 && isstruct(ph) && isfield(ph, 'compensate_xz')
            tf = logical(ph.compensate_xz);
        else
            % Unknown/custom phantoms retain the original compensated path.
            tf = true;
        end
    otherwise
        error('src:resolve_xz_mode:badValue', ...
            'compensate_xz must be auto, on, off, true, or false.');
end
end
