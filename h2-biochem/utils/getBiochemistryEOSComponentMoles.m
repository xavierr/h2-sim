function componentMoles = getBiochemistryEOSComponentMoles(model, state)
%GETBIOCHEMISTRYEOSCOMPONENTMOLES Reconstruct EOS inventories per cell.

nc = model.G.cells.num;
pv = value(model.PVTPropertyFunctions.get(model, state, 'PoreVolume'));
pv = expandCellVector(pv, nc, 'pore volume');

s = value(state.s);
liquid = model.getLiquidIndex();
vapor = model.getVaporIndex();
if iscell(s)
    sL = s{liquid};
    sV = s{vapor};
else
    sL = s(:, liquid);
    sV = s(:, vapor);
end
sL = expandCellVector(sL, nc, 'liquid saturation');
sV = expandCellVector(sV, nc, 'vapor saturation');

propertyModel = model.EOSModel.PropertyModel;
rhoL = propertyModel.computeMolarDensity(model.EOSModel, ...
    value(state.pressure), value(state.x), value(state.Z_L), ...
    value(state.T), true);
rhoV = propertyModel.computeMolarDensity(model.EOSModel, ...
    value(state.pressure), value(state.y), value(state.Z_V), ...
    value(state.T), false);
rhoL = expandCellVector(value(rhoL), nc, 'liquid molar density');
rhoV = expandCellVector(value(rhoV), nc, 'vapor molar density');

totalMoles = pv.*(sL.*rhoL + sV.*rhoV);
components = value(state.components);
assert(isequal(size(components), ...
    [nc, model.EOSModel.getNumberOfComponents()]), ...
    'EOS components must contain one row per cell and column per component.');
componentMoles = bsxfun(@times, totalMoles, components);
assert(all(isfinite(componentMoles(:)) & componentMoles(:) >= 0), ...
    'EOS component inventories must be finite and nonnegative.');
end

function values = expandCellVector(values, nc, name)
if isscalar(values)
    values = repmat(values, nc, 1);
else
    values = values(:);
end
assert(numel(values) == nc && all(isfinite(values)), ...
    'Expected %d finite %s values.', nc, name);
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.
%}
