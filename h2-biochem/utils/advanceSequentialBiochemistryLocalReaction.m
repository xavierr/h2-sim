function [state, report] = advanceSequentialBiochemistryLocalReaction( ...
        model, state, dt)
%ADVANCESEQUENTIALBIOCHEMISTRYLOCALREACTION Bounded local kinetic update.
%
% Growth and decay are integrated analytically with rates frozen at the
% start of this local substep. Reaction extents are reduced before the
% update so H2, inorganic carbon, and sulfate cannot become negative.

validateattributes(dt, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'dt');
nc = model.G.cells.num;
bcrm = model.biochemFluid;
nreact = bcrm.nbioreact;
names = model.EOSModel.getComponentNames();
molarMass = model.EOSModel.CompositionalMixture.molarMass;
gamma = model.gammak;

state = clearStateFunctionCaches(model, state);
growth = value(model.FacilityModel.getProps(state, 'PsiGrowthRate'));
decay = value(model.FacilityModel.getProps(state, 'PsiDecayRate'));
bmass = value(model.PVTPropertyFunctions.get(model, state, 'BacterialMass'));
growth = asReactionMatrix(growth, nc, nreact, 'growth rate');
decay = asReactionMatrix(decay, nc, nreact, 'decay rate');
bmass = asReactionMatrix(bmass, nc, nreact, 'bacterial mass');

[pv, sL, rhoL] = liquidStorage(model, state);
if isfield(state, 'sequentialLocalReactionTargetComponentMoles')
    componentMoles = ...
        value(state.sequentialLocalReactionTargetComponentMoles);
else
    componentMoles = getBiochemistryEOSComponentMoles(model, state);
end
idxH2 = findComponent(names, {'H2', 'Hydrogen'}, true);
idxCO2 = findComponent(names, {'CO2', 'CarbonDioxide'}, false);
hco3Moles = tracerMoles(state, 'tracerHCO3', pv, sL);
so4Moles = tracerMoles(state, 'tracerSO4', pv, sL);
hsMoles = tracerMoles(state, 'tracerHS', pv, sL);
so4Moles = so4Moles + sulfateDissolutionExtent( ...
    model, state, sL, dt);
% Methanogens/acetogens may draw inorganic carbon from the dolomite
% inventory (2 mol C per mol CaMg(CO3)2), releasing the stoichiometric
% Ca/Mg -- element totals equivalent to PHREEQC dissolving the mineral
% itself, so the per-substep equilibrium re-precipitates any overdraw.
% This in-substep mass transfer is what the compositional backend gets
% from solving EQUILIBRIUM_PHASES with KINETICS: the aqueous pool alone
% cannot carry the carbon flux. It is safe ONLY because (a) the
% equilibrium/feedback refresh runs every substep and (b) the
% equilibrium call is aqueous-only and idempotent (no GAS_PHASE carbon
% pump), so the feedback DIC that throttles the carbon Monod factor
% stays at the honest dolomite-buffered level.
if isfield(state, 'phreeqcMineralDolomite')
    dolomiteMoles = max(value(state.phreeqcMineralDolomite), 0);
else
    dolomiteMoles = zeros(nc, 1);
end
dolomiteCarbon = 2.*dolomiteMoles;
state.carbonSubstrateDt = dt;
state.carbonSubstrateMoles = hco3Moles + dolomiteCarbon;
if ~isempty(idxCO2)
    state.carbonSubstrateMoles = ...
        state.carbonSubstrateMoles + componentMoles(:, idxCO2);
end

qbaseCandidate = integratedReactionBase( ...
    growth, decay, bmass, dt, bcrm.Y_H2);
demandH2 = zeros(nc, nreact);
demandCarbon = zeros(nc, nreact);
demandSO4 = zeros(nc, nreact);
for reactionNo = 1:nreact
    demandH2(:, reactionNo) = bcrm.nbactMax(reactionNo).* ...
        qbaseCandidate(:, reactionNo)./rhoL;
    if ~isempty(idxCO2) && gamma(reactionNo, idxCO2) < 0
        demandCarbon(:, reactionNo) = ...
            -gamma(reactionNo, idxCO2)./abs(gamma(reactionNo, idxH2)).* ...
            demandH2(:, reactionNo);
    end
    if strcmpi(bcrm.metabolicReaction{reactionNo}, ...
            'SulfateReducingBacteria')
        demandSO4(:, reactionNo) = ...
            -bcrm.gamrsub(reactionNo)./abs(bcrm.gamrH2(reactionNo)).* ...
            demandH2(:, reactionNo);
    end
end

extentScale = ones(nc, nreact);
extentScale = applySharedLimit(extentScale, demandH2, ...
    componentMoles(:, idxH2), true(1, nreact));
carbonReactions = demandCarbon(1, :) > 0 | ...
    any(demandCarbon > 0, 1);
if any(carbonReactions)
    carbonAvailable = hco3Moles + dolomiteCarbon;
    if ~isempty(idxCO2)
        carbonAvailable = carbonAvailable + componentMoles(:, idxCO2);
    end
    extentScale = applySharedLimit(extentScale, demandCarbon, ...
        carbonAvailable, carbonReactions);
end
srbReactions = demandSO4(1, :) > 0 | any(demandSO4 > 0, 1);
if any(srbReactions)
    extentScale = applySharedLimit(extentScale, demandSO4, ...
        so4Moles, srbReactions);
end

% Limit the integrated reaction extent, not the instantaneous growth rate.
% Applying extentScale to mu and exponentiating a second time suppresses
% high-rate reactions by orders of magnitude when a resource is limiting.
qbaseExtent = qbaseCandidate.*extentScale;
decayedBmass = bmass.*exp((-decay - 1e-10).*dt);
unlimitedBmass = bmass.*exp((growth - decay - 1e-10).*dt);
targetBmass = decayedBmass + ...
    extentScale.*(unlimitedBmass - decayedBmass);

deltaMoles = zeros(size(componentMoles));
h2Consumed = zeros(nc, nreact);
fH2S = value(model.EOSModel.fractionH2SVolatile( ...
    value(state.T), getStatePH(state, nc)));
for reactionNo = 1:nreact
    gammaNorm = bcrm.nbactMax(reactionNo).*gamma(reactionNo, :).* ...
        molarMass./abs(gamma(reactionNo, idxH2));
    for componentNo = 1:numel(names)
        contribution = gammaNorm(componentNo).* ...
            qbaseExtent(:, reactionNo)./rhoL./molarMass(componentNo);
        if strcmpi(bcrm.metabolicReaction{reactionNo}, ...
                'SulfateReducingBacteria') && ...
                strcmpi(names{componentNo}, 'H2S')
            contribution = contribution.*fH2S;
        end
        deltaMoles(:, componentNo) = ...
            deltaMoles(:, componentNo) + contribution;
    end
    h2Consumed(:, reactionNo) = max( ...
        -gammaNorm(idxH2).*qbaseExtent(:, reactionNo)./ ...
        rhoL./molarMass(idxH2), 0);
end

if ~isempty(idxCO2)
    carbonDemand = max(-deltaMoles(:, idxCO2), 0);
    hco3Consumed = min(carbonDemand, hco3Moles);
    hco3Moles = hco3Moles - hco3Consumed;
    % Dolomite supplies only the overdraft AFTER both aqueous pools
    % (HCO3- tracer and EOS CO2) are drained, releasing the
    % stoichiometric Ca and Mg -- the same element totals as if PHREEQC
    % had dissolved the mineral itself. Draining the aqueous CO2 first
    % matters: CO2(aq) is the carbonic acid, and its actual removal is
    % what lets the per-substep equilibrium raise pH toward the
    % dolomite-buffered state (pH ~7.7, DIC ~3e-5) whose low carbon
    % Monod factor rate-limits methanogenesis exactly as the
    % compositional backend's in-PHREEQC kinetics are rate-limited.
    % Refunding dolomite carbon into the CO2 pool instead keeps the acid
    % in solution, pins pH near its initial value, and overpredicts
    % methanogenic consumption several-fold.
    dolomiteCarbonUsed = min(max(carbonDemand - hco3Consumed - ...
        componentMoles(:, idxCO2), 0), dolomiteCarbon);
    dolomiteDissolved = dolomiteCarbonUsed./2;
    dolomiteMoles = max(dolomiteMoles - dolomiteDissolved, 0);
    deltaMoles(:, idxCO2) = deltaMoles(:, idxCO2) + hco3Consumed + ...
        dolomiteCarbonUsed;
    if isfield(state, 'phreeqcMineralDolomite') && any(dolomiteCarbonUsed > 0)
        liquidVolume = max(pv.*sL, 1e-30);
        state.phreeqcMineralDolomite = dolomiteMoles;
        state.tracerCa = (tracerMoles(state, 'tracerCa', pv, sL) + ...
            dolomiteDissolved)./liquidVolume;
        state.tracerMg = (tracerMoles(state, 'tracerMg', pv, sL) + ...
            dolomiteDissolved)./liquidVolume;
    end
end

idxSRB = find(strcmpi(bcrm.metabolicReaction, ...
    'SulfateReducingBacteria'), 1);
if ~isempty(idxSRB)
    h2ExtentSRB = h2Consumed(:, idxSRB);
    so4Consumed = -bcrm.gamrsub(idxSRB)./ ...
        abs(bcrm.gamrH2(idxSRB)).*h2ExtentSRB;
    sulfideProduced = bcrm.gamp2(idxSRB)./ ...
        abs(bcrm.gamrH2(idxSRB)).*h2ExtentSRB;
    so4Moles = max(so4Moles - so4Consumed, 0);
    hsMoles = hsMoles + (1 - fH2S).*sulfideProduced;
end

componentMoles = max(componentMoles + deltaMoles, 0);
state.tracerHCO3 = hco3Moles./max(pv.*sL, 1e-30);
state.tracerSO4 = so4Moles./max(pv.*sL, 1e-30);
state.tracerHS = hsMoles./max(pv.*sL, 1e-30);
state.sequentialLocalReactionTargetComponentMoles = componentMoles;
state.components = bsxfun(@rdivide, componentMoles, ...
    sum(componentMoles, 2));
state = clearStateFunctionCaches(model, state);
state = model.computeFlash(state, inf);

[pv, sL, rhoL] = liquidStorage(model, state);
nbact = targetBmass./max(pv.*sL.*rhoL, 1e-30);
nbact = min(max(nbact, model.bact_capProp), model.bact_maxProp);
state.nbact = nbact;
state.h2ConsumptionRate = h2Consumed./dt;
if isfield(state, ...
        'sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles')
    cumulative = value( ...
        state.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles);
else
    cumulative = zeros(nc, nreact);
end
state.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles = ...
    cumulative + h2Consumed;
state = clearStateFunctionCaches(model, state);

report = struct( ...
    'Converged', true, ...
    'Failure', false, ...
    'Method', 'bounded-local-ode', ...
    'AcceptedSteps', 1, ...
    'RejectedSteps', 0, ...
    'ExtentScaleMinimum', min(extentScale(:)), ...
    'H2ConsumedMoles', sum(h2Consumed(:)));
end

function qbaseExtent = integratedReactionBase(mu, decay, bmass, dt, yield)
netRateDt = (mu - decay - 1e-10).*dt;
averageFactor = exponentialAverage(netRateDt);
qbaseExtent = bsxfun(@rdivide, ...
    mu.*bmass.*averageFactor.*dt, reshape(yield, 1, []));
end

function average = exponentialAverage(x)
average = zeros(size(x));
small = abs(x) < 1e-7;
average(~small) = (exp(x(~small)) - 1)./x(~small);
xs = x(small);
average(small) = 1 + xs./2 + xs.^2./6;
end

function scales = applySharedLimit(scales, demand, available, reactions)
totalDemand = sum(demand(:, reactions), 2);
resourceScale = min(1, available./max(totalDemand, 1e-30));
scales(:, reactions) = min(scales(:, reactions), ...
    repmat(resourceScale, 1, nnz(reactions)));
end

function [pv, sL, rhoL] = liquidStorage(model, state)
pv = value(model.PVTPropertyFunctions.get(model, state, 'PoreVolume'));
s = value(state.s);
rho = value(model.PVTPropertyFunctions.get(model, state, 'Density'));
liquid = model.getLiquidIndex();
if iscell(s)
    sL = s{liquid};
else
    sL = s(:, liquid);
end
if iscell(rho)
    rhoL = rho{liquid};
else
    rhoL = rho(:, liquid);
end
pv = pv(:);
sL = max(sL(:), 1e-12);
rhoL = max(rhoL(:), 1e-12);
end

function moles = tracerMoles(state, field, pv, sL)
if isfield(state, field)
    moles = max(value(state.(field)), 0).*pv.*sL;
else
    moles = zeros(size(pv));
end
end

function extent = sulfateDissolutionExtent(model, state, sL, dt)
extent = zeros(model.G.cells.num, 1);
if ~model.sulfateReduction || ~model.enableSulfateSource
    return;
end
rate = model.FlowDiscretization.getStateFunction('SRBTracerConvRate');
so4 = max(value(state.tracerSO4), 0);
drive = 1 - so4./rate.C_eq;
smoothDrive = 0.5.*(drive + sqrt(drive.^2 + 1e-6));
rateDensity = rate.k_dissolve.*rate.specific_surface_area.* ...
    sL.^rate.sw_exponent.*smoothDrive;
extent = rateDensity.*model.G.cells.volumes.*dt;
end

function values = asReactionMatrix(values, nc, nreact, name)
if iscell(values)
    values = cell2mat(cellfun(@(x) value(x(:)), values, ...
        'UniformOutput', false));
end
assert(isequal(size(values), [nc, nreact]) && ...
    all(isfinite(values(:))), ...
    'Expected a finite %d-by-%d %s matrix.', nc, nreact, name);
end

function index = findComponent(names, aliases, required)
index = [];
for i = 1:numel(aliases)
    index = find(strcmpi(names, aliases{i}), 1);
    if ~isempty(index)
        break
    end
end
assert(~required || ~isempty(index), ...
    'Required EOS component %s was not found.', aliases{1});
end

function pH = getStatePH(state, nc)
if isfield(state, 'phreeqcPH')
    pH = value(state.phreeqcPH);
else
    pH = [];
end
if ~isempty(pH)
    pH = reshape(pH, nc, 1);
end
end

function state = clearStateFunctionCaches(model, state)
groups = model.getStateFunctionGroupings();
for i = 1:numel(groups)
    name = groups{i}.getStateFunctionContainerName();
    if isfield(state, name)
        state = rmfield(state, name);
    end
end
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.
%}
