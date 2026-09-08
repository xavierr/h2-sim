function feedback = buildSequentialH2BiochemPhreeqcChemistryFeedback(model, state)
% Build a fixed numerical PHREEQC snapshot for hybrid kinetic evaluation.

nc = model.G.cells.num;
rhoWater = model.EOSModel.rho_water;
feedback = struct( ...
    'pH', stateVector(state, 'phreeqcPH', nc, model.carbonateBufferPH), ...
    'carbonatePka1', stateVector(state, 'phreeqcCarbonatePka1', nc, ...
        model.carbonateBufferPka1), ...
    'hco3Molality', [], ...
    'co2Molality', [], ...
    'totalCarbonMolality', [], ...
    'tds', stateVector(state, 'phreeqcTDS', nc, 0), ...
    'sulfateMolality', stateVector(state, 'tracerSO4', nc, 0)./rhoWater);

if isfield(state, 'phreeqcHCO3Molality')
    feedback.hco3Molality = stateVector( ...
        state, 'phreeqcHCO3Molality', nc, 0);
else
    feedback.hco3Molality = stateVector( ...
        state, 'tracerHCO3', nc, 0)./rhoWater;
end
if isfield(state, 'phreeqcCO2Molality')
    feedback.co2Molality = stateVector( ...
        state, 'phreeqcCO2Molality', nc, 0);
else
    feedback.co2Molality = dissolvedCO2Molality(model, state, rhoWater);
end
if isfield(state, 'phreeqcTotalCarbon')
    feedback.totalCarbonMolality = stateVector( ...
        state, 'phreeqcTotalCarbon', nc, 0);
else
    feedback.totalCarbonMolality = feedback.co2Molality + ...
        feedback.hco3Molality;
end
end

function molality = dissolvedCO2Molality(model, state, rhoWater)
names = model.EOSModel.getComponentNames();
index = find(strcmpi(names, 'CO2'), 1);
assert(~isempty(index), ...
    'sequential-h2biochem-phreeqc requires an EOS CO2 component.');
x = value(state.x);
if iscell(x)
    xCO2 = x{index};
else
    xCO2 = x(:, index);
end
rhoMolar = value(state.pressure)./( ...
    value(state.Z_L).*8.314.*value(state.T));
molality = max(rhoMolar.*value(xCO2)./rhoWater, 0);
end

function values = stateVector(state, field, nc, default)
if isfield(state, field)
    values = value(state.(field));
else
    values = default;
end
if isscalar(values)
    values = repmat(values, nc, 1);
else
    values = values(:);
end
assert(numel(values) == nc && isreal(values) && all(isfinite(values)), ...
    ['sequential-h2biochem-phreeqc feedback field %s must be a ', ...
     'finite cell vector.'], field);
end
