function [headers, values] = parseH2StorageIPhreeqcCOMSelectedOutput(raw, context)
% Normalize a GetSelectedOutputArray result from IPhreeqcCOM.
%
% Headings are normalized by lower-casing and removing punctuation. This
% makes PHREEQC variants such as "mass_H2O", "mass H2O", and "massh2o"
% addressable through the same alias without relying on selected-output
% column positions.

if nargin < 2 || isempty(context)
    context = 'IPhreeqcCOM selected output';
end
assert(iscell(raw) && size(raw, 1) >= 2 && size(raw, 2) >= 1, ...
    '%s does not contain a header row and at least one result row.', context);

headers = cellfun(@(x) lower(regexprep(char(x), '[^a-zA-Z0-9]', '')), ...
    raw(1, :), 'UniformOutput', false);
assert(numel(unique(headers)) == numel(headers), ...
    '%s contains ambiguous selected-output headings.', context);

data = raw(2:end, :);
values = zeros(size(data));
for i = 1:numel(data)
    entry = data{i};
    if ischar(entry) || (isstring(entry) && isscalar(entry))
        entry = str2double(entry);
    end
    assert(isnumeric(entry) && isscalar(entry) && isfinite(entry), ...
        '%s contains a non-finite numeric entry.', context);
    values(i) = entry;
end
end
