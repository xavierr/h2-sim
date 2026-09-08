function state = checkPhreeqcElementBalance(state, inputInventory, outputInventory, ...
    elements, options, backend)
% Validate and store per-cell PHREEQC elemental conservation diagnostics.

assert(isequal(size(inputInventory), size(outputInventory)) && ...
    size(inputInventory, 2) == numel(elements), ...
    'PHREEQC input/output element inventories must be nc-by-nelement arrays.');
assert(all(isfinite(inputInventory(:))) && all(isfinite(outputInventory(:))), ...
    'PHREEQC element inventories must be finite.');

ne = numel(elements);
absTol = tolerance(options, 'phreeqcElementBalanceAbsoluteTolerance', 1e-7, ne);
relTol = tolerance(options, 'phreeqcElementBalanceRelativeTolerance', 1e-8, ne);
absoluteResidual = abs(outputInventory - inputInventory);
scale = max(abs(inputInventory), abs(outputInventory));
normalizedResidual = zeros(size(scale));
nonzero = scale > 0;
normalizedResidual(nonzero) = absoluteResidual(nonzero)./scale(nonzero);
normalizedResidual(~nonzero & absoluteResidual > 0) = inf;
limit = bsxfun(@plus, absTol, bsxfun(@times, scale, relTol));
pass = absoluteResidual <= limit;

state.phreeqcElementBalanceElements = elements;
state.phreeqcElementBalanceInput = inputInventory;
state.phreeqcElementBalanceOutput = outputInventory;
state.phreeqcElementBalanceAbsoluteResidual = absoluteResidual;
state.phreeqcElementBalanceNormalizedResidual = normalizedResidual;
state.phreeqcElementBalancePass = pass;

if ~all(pass(:))
    [cellNo, elementNo] = find(~pass, 1);
    warning('H2Biochem:PhreeqcElementBalance', ...
        ['PHREEQC elemental conservation failed for backend "%s", cell %d, ', ...
         'element %s: input=%.16g mol, output=%.16g mol, absolute ', ...
         'residual=%.16g mol, normalized residual=%.16g (limit %.16g mol).'], ...
        backend, cellNo, elements{elementNo}, inputInventory(cellNo, elementNo), ...
        outputInventory(cellNo, elementNo), absoluteResidual(cellNo, elementNo), ...
        normalizedResidual(cellNo, elementNo), limit(cellNo, elementNo));
end
end

function values = tolerance(options, field, default, ne)
if isfield(options, field)
    values = options.(field);
else
    values = default;
end
validateattributes(values, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative'}, ...
    mfilename, field);
if isscalar(values)
    values = repmat(values, 1, ne);
else
    values = reshape(values, 1, []);
    assert(numel(values) == ne, ...
        '%s must be scalar or contain one value per audited element.', field);
end
end
