function selected = getH2StorageIPhreeqcCOMSelectedOutputValue(headers, values, aliases, context)
% Return the last selected-output value matching one normalized alias.

if nargin < 4 || isempty(context)
    context = 'IPhreeqcCOM selected output';
end
assert(iscell(aliases) && ~isempty(aliases), ...
    'aliases must be a nonempty cell array.');
aliases = cellfun(@(x) lower(regexprep(char(x), '[^a-zA-Z0-9]', '')), ...
    aliases, 'UniformOutput', false);
column = find(ismember(headers, aliases), 1);
assert(~isempty(column), ...
    '%s is missing "%s". Available columns: %s', context, aliases{1}, ...
    strjoin(headers, ', '));
selected = values(end, column);
end
