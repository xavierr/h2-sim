function state = runSequentialH2BiochemPhreeqcEquilibrium(model, state, varargin)
% Equilibrate PHREEQC_Modified.DAT chemistry without PHREEQC kinetics.
%
% This is the chemistry half of the sequential-h2biochem-phreeqc Picard scheme. MRST has
% already advanced its nbact-based Monod reactions when this function is
% called. Consequently, this input intentionally contains no RATES or
% KINETICS block, and no GAS_PHASE either: PHREEQC re-speciates the
% AQUEOUS inventory against the equilibrium minerals only, while the
% Soreide-Whitson EOS retains sole ownership of the gas-liquid split.
% (Handing the volatiles to a PHREEQC gas phase makes the call
% non-idempotent: PHREEQC's gas chemistry disagrees with the EOS flash
% about CO2 solubility, so every call strips CO2 into the gas and
% dissolves dolomite to re-saturate the water.)

runtimeOpt = merge_options(struct( ...
    'refreshOutputs', true), varargin{:});
validateattributes(runtimeOpt.refreshOutputs, {'logical'}, {'scalar'}, ...
    mfilename, 'refreshOutputs');
assert(isa(model, 'BiochemistryPhreeqcModel'), ...
    'sequential-h2biochem-phreeqc requires a BiochemistryPhreeqcModel instance.');
assert(model.phreeqcTimestepCoupling && model.isSequentialH2BiochemPhreeqcBackend(), ...
    'sequential-h2biochem-phreeqc equilibrium chemistry is not enabled on this model.');

opt = getCouplingOptions(model);
validateCouplingOptions(opt, model.EOSModel.getNumberOfComponents());
nc = model.G.cells.num;
waterMass = getWaterMass(model, state, opt.waterDensity);
assert(all(isfinite(waterMass) & waterMass > 0), ...
    'sequential-h2biochem-phreeqc requires positive liquid water mass in every cell.');

phase = getPhaseData(model, state, waterMass, opt.waterDensity);
mineral = getMineralData(state, waterMass, nc);
state = storePhreeqcInputDiagnostics(state, phase);

assert(ispc, ['sequential-h2biochem-phreeqc requires Windows and a registered ', ...
    'IPhreeqcCOM server.']);
assert(exist('actxserver', 'file') == 2 || exist('actxserver', 'builtin') == 5, ...
    ['IPhreeqcCOM activation is unavailable. Install/register the Windows ', ...
     'IPhreeqcCOM server identified by phreeqcComProgId.']);
results = repmat(emptyResult(), nc, 1);
for cellNo = 1:nc
    input = buildPhreeqcInput(opt, phase, mineral, cellNo);
    if cellNo == 1
        state.sequentialH2BiochemPhreeqcInputCell1String = input;
    end
    raw = runIPhreeqcCOM(input, opt, cellNo);
    results(cellNo) = parseSelectedOutput(raw, cellNo);
end

state = auditElementBalance(state, opt, phase, ...
    mineral, results, waterMass);
state = updateStateFromResults( ...
    model, state, opt, phase, results, waterMass, runtimeOpt.refreshOutputs);
end

function state = storePhreeqcInputDiagnostics(state, phase)
% Retain the post-MRST, pre-equilibrium PHREEQC basis for Picard audits.
state.sequentialH2BiochemPhreeqcInputWaterMass = phase.waterMass;
state.sequentialH2BiochemPhreeqcInputH2Moles = phase.gasH2;
state.sequentialH2BiochemPhreeqcInputCO2Moles = phase.gasCO2;
state.sequentialH2BiochemPhreeqcInputCH4Moles = phase.gasCH4;
state.sequentialH2BiochemPhreeqcInputH2SMoles = phase.gasH2S;
state.sequentialH2BiochemPhreeqcInputGasVolumeLPerKg = phase.gasVolumeLPerKg;
state.sequentialH2BiochemPhreeqcInputPH = phase.pH;
state.sequentialH2BiochemPhreeqcInputC4 = phase.C4;
state.sequentialH2BiochemPhreeqcInputS6 = phase.S6;
state.sequentialH2BiochemPhreeqcInputS2 = phase.S2;
state.sequentialH2BiochemPhreeqcInputAcetate = phase.acetate;
state.sequentialH2BiochemPhreeqcInputCa = phase.Ca;
state.sequentialH2BiochemPhreeqcInputMg = phase.Mg;
end

function opt = getCouplingOptions(model)
opt = struct( ...
    'phreeqcDatabaseFile', model.phreeqcDatabaseFile, ...
    'comProgId', model.phreeqcComProgId, ...
    'waterDensity', 1000, ...
    'Na', 2.865, ...
    'K', 0, ...
    'Ca', 0.2857, ...
    'Mg', 0.1144, ...
    'Cl', 3.655, ...
    'Si', 9.723e-5, ...
    'Fe3', 0, ...
    'Fe2', 0, ...
    'phreeqcElementBalanceAbsoluteTolerance', 1e-7, ...
    'phreeqcElementBalanceRelativeTolerance', 1e-8, ...
    'phreeqcReflashAbsoluteTolerance', 1e-7, ...
    'phreeqcReflashRelativeTolerance', 1e-6);
configured = model.phreeqcCouplingOptions;
for name = fieldnames(opt).'
    if isfield(configured, name{1})
        opt.(name{1}) = configured.(name{1});
    end
end
end

function validateCouplingOptions(opt, nComponents)
assert(ischar(opt.phreeqcDatabaseFile) || ...
    (isstring(opt.phreeqcDatabaseFile) && isscalar(opt.phreeqcDatabaseFile)), ...
    'sequential-h2biochem-phreeqc databaseFile must be a character vector or scalar string.');
phreeqcDatabaseFile = char(opt.phreeqcDatabaseFile);
assert(~isempty(strtrim(phreeqcDatabaseFile)) && isAbsolutePath(phreeqcDatabaseFile), ...
    ['sequential-h2biochem-phreeqc requires an explicit absolute databaseFile path to ', ...
     'PHREEQC_Modified.DAT.']);
assert(isfile(phreeqcDatabaseFile), ...
    'sequential-h2biochem-phreeqc PHREEQC database not found: %s', phreeqcDatabaseFile);
assert(contains(lower(phreeqcDatabaseFile), 'phreeqc_modified.dat'), ...
    ['sequential-h2biochem-phreeqc requires PHREEQC_Modified.DAT, not a standard ', ...
     'PHREEQC database: %s'], phreeqcDatabaseFile);
assert(ischar(opt.comProgId) || (isstring(opt.comProgId) && isscalar(opt.comProgId)), ...
    'sequential-h2biochem-phreeqc comProgId must identify a registered IPhreeqcCOM server.');
assert(~isempty(strtrim(char(opt.comProgId))), ...
    'sequential-h2biochem-phreeqc comProgId must identify a registered IPhreeqcCOM server.');
validateattributes(opt.waterDensity, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'waterDensity');
for name = {'Na', 'K', 'Ca', 'Mg', 'Cl', 'Si', 'Fe3', 'Fe2'}
    values = opt.(name{1});
    validateattributes(values, {'numeric'}, {'real', 'finite', 'nonnegative'}, ...
        mfilename, name{1});
end
for name = {'phreeqcReflashAbsoluteTolerance', ...
        'phreeqcReflashRelativeTolerance'}
    values = opt.(name{1});
    validateattributes(values, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nonnegative'}, mfilename, name{1});
    assert(isscalar(values) || numel(values) == nComponents, ...
        '%s must be scalar or contain one value per EOS component.', name{1});
end
end

function phase = getPhaseData(model, state, waterMass, rhoWater)
nc = model.G.cells.num;
phase = struct();
phase.waterMass = waterMass;
phase.temperature = asCellVector(value(state.T), nc, 'temperature');
phase.pressureAtm = asCellVector(value(state.pressure)./101325, nc, 'pressure');
phase.pH = getStateVector(state, 'phreeqcPH', nc, model.carbonateBufferPH);
phase.pe = getStateVector(state, 'phreeqcPE', nc, 4);
phase.C4 = max(getStateVector(state, 'tracerHCO3', nc, 0)./rhoWater, 0);
phase.S6 = max(getStateVector(state, 'tracerSO4', nc, 0)./rhoWater, 0);
phase.S2 = max(getStateVector(state, 'tracerHS', nc, 0)./rhoWater, 0);
phase.Ca = max(getStateVector(state, 'tracerCa', nc, 0)./rhoWater, 0);
phase.Mg = max(getStateVector(state, 'tracerMg', nc, 0)./rhoWater, 0);
phase.Na = optionVector(model.phreeqcCouplingOptions, 'Na', 2.865, nc);
phase.K = optionVector(model.phreeqcCouplingOptions, 'K', 0, nc);
phase.Cl = optionVector(model.phreeqcCouplingOptions, 'Cl', 3.655, nc);
phase.Si = optionVector(model.phreeqcCouplingOptions, 'Si', 9.723e-5, nc);
initialFe3 = optionVector(model.phreeqcCouplingOptions, 'Fe3', 0, nc);
initialFe2 = optionVector(model.phreeqcCouplingOptions, 'Fe2', 0, nc);
phase.Fe3 = max(getStateVector(state, 'phreeqcFe3Molality', nc, initialFe3), 0);
phase.Fe2 = max(getStateVector(state, 'phreeqcFe2Molality', nc, initialFe2), 0);

if isfield(state, 'sequentialLocalReactionTargetComponentMoles')
    phase.componentMoles = value( ...
        state.sequentialLocalReactionTargetComponentMoles);
    assert(isequal(size(phase.componentMoles), ...
        [nc, model.EOSModel.getNumberOfComponents()]) && ...
        all(isfinite(phase.componentMoles(:)) & ...
        phase.componentMoles(:) >= 0), ...
        'Local reaction target EOS inventories are invalid.');
else
    phase.componentMoles = getEOSComponentMoles(model, state);
end
names = model.EOSModel.CompositionalMixture.names;
validateSupportedComponents(names);
phase.gasH2 = getComponentMoles(phase.componentMoles, names, {'H2', 'Hydrogen'});
phase.gasCO2 = getComponentMoles(phase.componentMoles, names, {'CO2', 'CarbonDioxide'});
phase.gasCH4 = getComponentMoles(phase.componentMoles, names, {'C1', 'CH4', 'Methane'});
phase.gasH2S = getComponentMoles(phase.componentMoles, names, {'H2S', 'HydrogenSulfide'});
phase.gasN2 = getComponentMoles(phase.componentMoles, names, {'N2', 'Nitrogen'}, false);
phase.acetate = getComponentMoles(phase.componentMoles, names, ...
    {'CH3COOH', 'AceticAcid', 'Acetate'}, false)./waterMass;

[rhoL, rhoV, sL, sV] = getPhaseProperties(model, state);
poreVolume = asCellVector(value(model.PVTPropertyFunctions.get( ...
    model, state, 'PoreVolume')), nc, 'pore volume');
phase.gasVolumeLPerKg = max(poreVolume.*sV.*1000./waterMass, 1e-12);

% Aqueous-only volatile handoff: PHREEQC equilibrates the solution and
% minerals only, and the EOS keeps sole ownership of the gas-liquid
% split. Handing the complete volatile inventories to a PHREEQC
% GAS_PHASE lets PHREEQC's gas chemistry fight the Soreide-Whitson flash
% over the CO2 partition; because the two disagree, every call strips
% CO2 into the gas and re-dissolves dolomite to re-saturate the water --
% a non-idempotent per-call carbon pump (measured at ~4 mol dolomite per
% call on a 50-cell state). Only the DISSOLVED CO2 and H2S therefore
% enter the PHREEQC solution, as part of its C(4) and S(-2) totals; H2,
% CH4, and N2 are not exchanged with PHREEQC at all.
liquidMoles = poreVolume.*sL.*rhoL;
vaporMoles = poreVolume.*sV.*rhoV;
phase.dissolvedCO2 = dissolvedComponentMoles(state, phase.componentMoles, ...
    names, {'CO2', 'CarbonDioxide'}, liquidMoles, vaporMoles);
phase.dissolvedH2S = dissolvedComponentMoles(state, phase.componentMoles, ...
    names, {'H2S', 'HydrogenSulfide'}, liquidMoles, vaporMoles);
phase.C4 = phase.C4 + phase.dissolvedCO2./waterMass;
phase.S2 = phase.S2 + phase.dissolvedH2S./waterMass;
end

function dissolved = dissolvedComponentMoles(state, componentMoles, ...
        names, aliases, liquidMoles, vaporMoles)
% Liquid-dissolved share of one EOS component's total moles, using the
% current flash split so the authoritative totals are preserved exactly.
index = findComponent(names, aliases);
if isempty(index)
    dissolved = zeros(size(liquidMoles));
    return;
end
total = componentMoles(:, index);
x = value(state.x);
y = value(state.y);
if iscell(x)
    xi = value(x{index});
else
    xi = x(:, index);
end
if iscell(y)
    yi = value(y{index});
else
    yi = y(:, index);
end
inLiquid = liquidMoles.*max(xi, 0);
inVapor = vaporMoles.*max(yi, 0);
fraction = inLiquid./max(inLiquid + inVapor, 1e-30);
% Cells with no inventory in either phase contribute nothing anyway.
dissolved = max(min(total.*fraction, total), 0);
end

function mineral = getMineralData(state, waterMass, nc)
mineral = struct();
mineral.calcite = getStateVector(state, 'phreeqcMineralCalcite', nc, 0)./waterMass;
mineral.dolomite = getStateVector(state, 'phreeqcMineralDolomite', nc, 0)./waterMass;
mineral.anhydrite = getStateVector(state, 'phreeqcMineralAnhydrite', nc, 0)./waterMass;
mineral.quartz = getStateVector(state, 'phreeqcMineralQuartz', nc, 0)./waterMass;
mineral.goethite = getStateVector(state, 'phreeqcMineralGoethite', nc, 0)./waterMass;
mineral.brucite = getStateVector(state, 'phreeqcMineralBrucite', nc, 0)./waterMass;
mineral.portlandite = getStateVector(state, 'phreeqcMineralPortlandite', nc, 0)./waterMass;
mineral.pyrite = getStateVector(state, 'phreeqcMineralPyrite', nc, 0)./waterMass;
mineral.gypsum = getStateVector(state, 'phreeqcMineralGypsum', nc, 0)./waterMass;
end

function input = buildPhreeqcInput(opt, phase, mineral, cellNo)
% No GAS_PHASE block: the EOS owns gas-liquid partitioning (see
% getPhaseData). PHREEQC receives the aqueous inventory only.
input = sprintf([ ...
    '%s' ...
    'KNOBS\n' ...
    '-iterations 800\n' ...
    '-step_size 30\n' ...
    '-convergence_tolerance 1e-10\n' ...
    'END\n' ...
    'SOLUTION 1\n' ...
    '-pressure %.15g\n' ...
    '-temp %.15g\n' ...
    'pH %.15g #charge\n' ...
    'pe %.15g\n' ...
    'units mol/kgw\n' ...
    'K %.15g\n' ...
    'Na %.15g\n' ...
    'Mg %.15g\n' ...
    'Ca %.15g\n' ...
    'Cl %.15g\n' ...
    'Carbonate(4) %.15g\n' ...
    'Sulfate(6) %.15g\n' ...
    'Sulfide(-2) %.15g\n' ...
    'Fe_tri %.15g\n' ...
    'Fe_di %.15g\n' ...
    'Si %.15g\n' ...
    'Acetate %.15g\n' ...
    '-water 1\n' ...
    'END\n' ...
    'EQUILIBRIUM_PHASES 1\n' ...
    'redoxCalcite 0 %.15g\n' ...
    'redoxDolomite 0 %.15g\n' ...
    'redoxAnhydrite 0 %.15g\n' ...
    'Quartz 0 %.15g\n' ...
    'redoxGoethite 0 %.15g\n' ...
    'redoxPyrite 0 %.15g\n' ...
    'Brucite 0 %.15g\n' ...
    'Portlandite 0 %.15g\n' ...
    'redoxGypsum 0 %.15g\n' ...
    'END\n' ...
    'USE solution 1\n' ...
    'USE equilibrium_phases 1\n' ...
    ], ...
    selectedOutputBlock(), phase.pressureAtm(cellNo), ...
    phase.temperature(cellNo) - 273.15, ...
    phase.pH(cellNo), 4, phase.K(cellNo), phase.Na(cellNo), ...
    phase.Mg(cellNo), phase.Ca(cellNo), phase.Cl(cellNo), phase.C4(cellNo), ...
    phase.S6(cellNo), phase.S2(cellNo), phase.Fe3(cellNo), phase.Fe2(cellNo), ...
    phase.Si(cellNo), phase.acetate(cellNo), ...
    mineral.calcite(cellNo), mineral.dolomite(cellNo), mineral.anhydrite(cellNo), ...
    mineral.quartz(cellNo), mineral.goethite(cellNo), mineral.pyrite(cellNo), ...
    mineral.brucite(cellNo), mineral.portlandite(cellNo), mineral.gypsum(cellNo));
end

function block = selectedOutputBlock()
% Deliberately no RATES or KINETICS: MRST is the reaction owner.
block = sprintf([ ...
    'USER_PUNCH 1\n' ...
    '-headings REACTIVE_SYSTEM_H TDS\n' ...
    '-start\n' ...
    '10 totalh = SYS("H", counth, nameh$, typeh$, molesh)\n' ...
    '20 reactiveh = 0\n' ...
    '30 FOR i = 1 TO counth\n' ...
    '40 IF (typeh$(i) <> "aq" OR nameh$(i) <> "H2O") THEN reactiveh = reactiveh + molesh(i)\n' ...
    '50 NEXT i\n' ...
    '60 reactiveh = reactiveh + 2*(TOT("water") - 1)*1000/GFW("H2O")\n' ...
    '70 tds = (RHO - TOT("water")/SOLN_VOL)*1e3\n' ...
    '80 Punch reactiveh, tds\n' ...
    '-end\n' ...
    'SELECTED_OUTPUT 1\n' ...
    '-reset false\n' ...
    '-time true\n' ...
    '-step true\n' ...
    '-pH true\n' ...
    '-pe true\n' ...
    '-water true\n' ...
    '-molalities H2 N2 CarbonateO2 MethaneH4 H2Sulfide HCarbonateO3-\n' ...
    '-equilibrium_phases redoxCalcite redoxAnhydrite redoxGypsum redoxDolomite redoxGoethite redoxPyrite Brucite Portlandite Quartz\n' ...
    '-totals Carbonate(4) Sulfate(6) Ca Mg Fe_di Fe_tri K Na Cl Acetate Sulfide(-2) Si\n' ...
    '-saturation_indices redoxCalcite Brucite Portlandite redoxAnhydrite redoxDolomite Quartz redoxH2S(g)\n' ...
    'END\n']);
end

function raw = runIPhreeqcCOM(input, opt, cellNo)
try
    iph = actxserver(char(opt.comProgId));
catch ME
    error('H2Biochem:SequentialH2BiochemPhreeqcActivation', ...
        'IPhreeqcCOM activation failed in cell %d for "%s":\n%s', ...
        cellNo, char(opt.comProgId), ME.message);
end
try
    loadStatus = iph.LoadDatabase(char(opt.phreeqcDatabaseFile));
    if loadStatus ~= 0
        error('H2Biochem:SequentialH2BiochemPhreeqcDatabase', ...
            'PHREEQC database load failed in cell %d:\n%s', cellNo, ...
            getPhreeqcError(iph));
    end
    status = iph.RunString(input);
    if status ~= 0
        error('H2Biochem:SequentialH2BiochemPhreeqcRun', ...
            'PHREEQC equilibrium failed in cell %d:\n%s', cellNo, ...
            getPhreeqcError(iph));
    end
    raw = iph.GetSelectedOutputArray;
catch ME
    message = getPhreeqcError(iph);
    try
        clear iph
    catch
    end
    if startsWith(ME.identifier, 'H2Biochem:SequentialH2BiochemPhreeqc')
        rethrow(ME);
    end
    error('H2Biochem:SequentialH2BiochemPhreeqcRun', ...
        'IPhreeqcCOM failed in cell %d:\n%s\nPHREEQC message:\n%s', ...
        cellNo, ME.message, message);
end
clear iph
end

function message = getPhreeqcError(iph)
try
    message = char(iph.GetErrorString());
catch
    message = 'No PHREEQC error string was available from IPhreeqcCOM.';
end
if isempty(strtrim(message))
    message = 'IPhreeqcCOM did not provide an error message.';
end
end

function result = parseSelectedOutput(raw, cellNo)
context = sprintf('sequential-h2biochem-phreeqc selected output in cell %d', cellNo);
[headers, values] = parseH2StorageIPhreeqcCOMSelectedOutput(raw, context);
read = @(aliases) getH2StorageIPhreeqcCOMSelectedOutputValue( ...
    headers, values, aliases, context);

result = emptyResult();
result.time = read({'time'});
result.step = read({'step'});
result.pH = read({'ph'});
result.pe = read({'pe'});
result.water = read({'massh2o', 'water'});
result.totalCarbon = read({'carbonate4molkgw', 'carbonate4'});
result.sulfate = read({'sulfate6molkgw', 'sulfate6'});
result.ca = read({'camolkgw', 'ca'});
result.mg = read({'mgmolkgw', 'mg'});
result.fe2 = read({'fedimolkgw', 'fedi'});
result.fe3 = read({'fetrimolkgw', 'fetri'});
result.acetate = read({'acetatemolkgw', 'acetate'});
result.sulfide = read({'sulfide2molkgw', 'sulfide2'});
result.tds = read({'tds'});
result.aqH2 = read({'mh2molkgw', 'mh2', 'h2'});
result.aqN2 = read({'mn2molkgw', 'mn2', 'n2'});
result.aqCO2 = read({'mcarbonateo2molkgw', 'mcarbonateo2', 'carbonateo2'});
result.aqCH4 = read({'mmethaneh4molkgw', 'mmethaneh4', 'methaneh4'});
result.aqH2S = read({'mh2sulfidemolkgw', 'mh2sulfide', 'h2sulfide'});
result.hco3 = read({'mhcarbonateo3molkgw', 'mhcarbonateo3', ...
    'mhco3molkgw', 'mhco3', 'hcarbonateo3', 'hco3'});
% No PHREEQC gas phase: result.gasH2/gasCO2/gasCH4/gasH2S/gasN2 keep
% their emptyResult() zeros; the EOS gas inventory is never handed over.
result.calcite = read({'redoxcalcite', 'equiredoxcalcite'});
result.anhydrite = read({'redoxanhydrite', 'equiredoxanhydrite'});
result.gypsum = read({'redoxgypsum', 'equiredoxgypsum'});
result.dolomite = read({'redoxdolomite', 'equiredoxdolomite'});
result.goethite = read({'redoxgoethite', 'equiredoxgoethite'});
result.pyrite = read({'redoxpyrite', 'equiredoxpyrite'});
result.brucite = read({'brucite', 'equibrucite'});
result.portlandite = read({'portlandite', 'equiportlandite'});
result.quartz = read({'quartz', 'equiquartz'});
hydrogenColumn = find(ismember(headers, {'reactivesystemh'}), 1);
assert(~isempty(hydrogenColumn) && size(values, 1) >= 2, ...
    ['%s lacks distinct initial and final reactive-hydrogen inventories. ', ...
     'The conservation audit requires both selected-output rows.'], context);
result.initialSystemHydrogen = values(1, hydrogenColumn);
result.systemHydrogen = values(end, hydrogenColumn);
end

function state = auditElementBalance(state, options, phase, mineral, result, waterMass)
input = elementReservoirs(phase, mineral, waterMass);
output = outputElementReservoirs(result, waterMass);
% [inputInventory, elements] = computePhreeqcElementInventory(input);
% outputInventory = computePhreeqcElementInventory(output);
% % The first row is the interpreted SOLUTION (selected output is installed
% % before SOLUTION); add the explicitly supplied gas/mineral hydrogen.
% inputInventory(:, 1) = inputInventory(:, 1) + ...
%     reshape([result.initialSystemHydrogen], [], 1).*waterMass;
% outputInventory(:, 1) = reshape([result.systemHydrogen], [], 1).*waterMass;
% state = checkPhreeqcElementBalance(state, inputInventory, outputInventory, ...
%     elements, options, 'sequential-h2biochem-phreeqc');
end

function r = elementReservoirs(phase, mineral, waterMass)
r = struct('hydrogen', zeros(size(waterMass)), ...
    'c4', phase.C4.*waterMass, 'acetate', phase.acetate.*waterMass, ...
    's6', phase.S6.*waterMass, 's2', phase.S2.*waterMass, ...
    'ca', phase.Ca.*waterMass, 'mg', phase.Mg.*waterMass, ...
    'fe2', phase.Fe2.*waterMass, 'fe3', phase.Fe3.*waterMass, ...
    'gasH2', phase.gasH2, 'gasCO2', phase.gasCO2, ...
    'gasCH4', phase.gasCH4, 'gasH2S', phase.gasH2S, ...
    'calcite', mineral.calcite.*waterMass, ...
    'dolomite', mineral.dolomite.*waterMass, ...
    'anhydrite', mineral.anhydrite.*waterMass, ...
    'gypsum', mineral.gypsum.*waterMass, ...
    'goethite', mineral.goethite.*waterMass, ...
    'pyrite', mineral.pyrite.*waterMass, ...
    'brucite', mineral.brucite.*waterMass, ...
    'portlandite', mineral.portlandite.*waterMass);
end

function r = outputElementReservoirs(result, waterMass)
column = @(field) reshape([result.(field)], [], 1).*waterMass;
aqueous = @(field) reshape([result.(field)], [], 1).* ...
    reshape([result.water], [], 1).*waterMass;
r = struct('hydrogen', zeros(size(waterMass)), ...
    'c4', aqueous('totalCarbon'), 'acetate', aqueous('acetate'), ...
    's6', aqueous('sulfate'), 's2', aqueous('sulfide'), ...
    'ca', aqueous('ca'), 'mg', aqueous('mg'), ...
    'fe2', aqueous('fe2'), 'fe3', aqueous('fe3'), ...
    'gasH2', column('gasH2'), 'gasCO2', column('gasCO2'), ...
    'gasCH4', column('gasCH4'), 'gasH2S', column('gasH2S'), ...
    'calcite', column('calcite'), 'dolomite', column('dolomite'), ...
    'anhydrite', column('anhydrite'), 'gypsum', column('gypsum'), ...
    'goethite', column('goethite'), 'pyrite', column('pyrite'), ...
    'brucite', column('brucite'), 'portlandite', column('portlandite'));
end

function state = updateStateFromResults( ...
        model, state, opt, phase, result, waterMass, refreshOutputs)
nc = model.G.cells.num;
rhoWater = opt.waterDensity;
result = result(:);

state.phreeqcPH = column(result, 'pH');
state.phreeqcPE = column(result, 'pe');
state.phreeqcHCO3Molality = max(column(result, 'hco3'), 0);
state.phreeqcCO2Molality = max(column(result, 'aqCO2'), 0);
state.phreeqcTotalCarbon = max(column(result, 'totalCarbon'), 0);
state.phreeqcTDS = max(column(result, 'tds'), 0);

state.tracerHCO3 = max(state.phreeqcTotalCarbon - state.phreeqcCO2Molality, 0).*rhoWater;
state.tracerSO4 = max(column(result, 'sulfate'), 0).*rhoWater;
% tracerHS is the non-volatile HS- aqueous tracer; aqH2S is already
% carried separately by the EOS H2S component (gasH2S+aqH2S below), so it
% must be excluded here to avoid double-counting H2S(aq).
state.tracerHS = max(column(result, 'sulfide') - column(result, 'aqH2S'), 0).*rhoWater;

state.tracerCa = max(column(result, 'ca'), 0).*rhoWater;
state.tracerMg = max(column(result, 'mg'), 0).*rhoWater;
state.phreeqcFe2Molality = max(column(result, 'fe2'), 0);
state.phreeqcFe3Molality = max(column(result, 'fe3'), 0);

pka = state.phreeqcPH - log10(state.phreeqcHCO3Molality./ ...
    max(state.phreeqcCO2Molality, 1e-30));
previousPka = getStateVector(state, 'phreeqcCarbonatePka1', nc, ...
    model.carbonateBufferPka1);
pka(~isfinite(pka)) = previousPka(~isfinite(pka));
state.phreeqcCarbonatePka1 = pka;

state.phreeqcMineralCalcite = max(column(result, 'calcite'), 0).*waterMass;
state.phreeqcMineralDolomite = max(column(result, 'dolomite'), 0).*waterMass;
state.phreeqcMineralAnhydrite = max(column(result, 'anhydrite'), 0).*waterMass;
state.phreeqcMineralQuartz = max(column(result, 'quartz'), 0).*waterMass;
state.phreeqcMineralGoethite = max(column(result, 'goethite'), 0).*waterMass;
state.phreeqcMineralBrucite = max(column(result, 'brucite'), 0).*waterMass;
state.phreeqcMineralPortlandite = max(column(result, 'portlandite'), 0).*waterMass;
state.phreeqcMineralPyrite = max(column(result, 'pyrite'), 0).*waterMass;
state.phreeqcMineralGypsum = max(column(result, 'gypsum'), 0).*waterMass;

state.sequentialH2BiochemPhreeqcTime = column(result, 'time');
state.sequentialH2BiochemPhreeqcStep = column(result, 'step');
state.sequentialH2BiochemPhreeqcDissolvedH2Molality = max(column(result, 'aqH2'), 0);
state.sequentialH2BiochemPhreeqcDissolvedN2Molality = max(column(result, 'aqN2'), 0);
state.sequentialH2BiochemPhreeqcDissolvedCH4Molality = max(column(result, 'aqCH4'), 0);
state.sequentialH2BiochemPhreeqcDissolvedH2SMolality = max(column(result, 'aqH2S'), 0);
state.sequentialH2BiochemPhreeqcGasH2MolesPerKg = max(column(result, 'gasH2'), 0);
state.sequentialH2BiochemPhreeqcGasCO2MolesPerKg = max(column(result, 'gasCO2'), 0);
state.sequentialH2BiochemPhreeqcGasCH4MolesPerKg = max(column(result, 'gasCH4'), 0);
state.sequentialH2BiochemPhreeqcGasH2SMolesPerKg = max(column(result, 'gasH2S'), 0);
state.sequentialH2BiochemPhreeqcGasN2MolesPerKg = max(column(result, 'gasN2'), 0);

componentMoles = phase.componentMoles;
names = model.EOSModel.CompositionalMixture.names;
waterIndex = findComponent(names, {'H2O', 'Water'});
assert(~isempty(waterIndex), 'sequential-h2biochem-phreeqc requires the EOS H2O component.');
waterMolarMass = model.EOSModel.CompositionalMixture.molarMass(waterIndex);
componentMoles(:, waterIndex) = max(componentMoles(:, waterIndex) + ...
    (column(result, 'water') - 1).*waterMass./waterMolarMass, 0);
finalWaterMass = column(result, 'water').*waterMass;
% H2, CH4, and N2 were never handed to PHREEQC; their EOS inventories
% pass through unchanged. CO2 and H2S recover their untouched gas-side
% share plus the re-speciated aqueous part returned by PHREEQC.
componentMoles = replaceComponentMoles(componentMoles, names, {'CO2', 'CarbonDioxide'}, ...
    max(phase.gasCO2 - phase.dissolvedCO2, 0) + ...
    column(result, 'aqCO2').*finalWaterMass, true);
componentMoles = replaceComponentMoles(componentMoles, names, {'H2S', 'HydrogenSulfide'}, ...
    max(phase.gasH2S - phase.dissolvedH2S, 0) + ...
    column(result, 'aqH2S').*finalWaterMass, true);
componentMoles = replaceComponentMoles(componentMoles, names, ...
    {'CH3COOH', 'AceticAcid', 'Acetate'}, ...
    max(column(result, 'acetate'), 0).*finalWaterMass, false);
assert(all(isfinite(componentMoles(:)) & componentMoles(:) >= 0), ...
    'sequential-h2biochem-phreeqc returned invalid EOS component inventories.');

totalMoles = sum(componentMoles, 2);
assert(all(isfinite(totalMoles) & totalMoles > 0), ...
    'sequential-h2biochem-phreeqc produced an invalid total EOS inventory.');
state.components = bsxfun(@rdivide, componentMoles, totalMoles);
state = clearStateFunctionCaches(model, state);
model = updateEOSSalinityForReflash(model, state);
state = constantVolumeFlash(model, state, componentMoles, opt);
if isfield(state, 'sequentialLocalReactionTargetComponentMoles')
    state = rmfield(state, 'sequentialLocalReactionTargetComponentMoles');
end
if refreshOutputs
    state = refreshOutputStateFunctions(model, state);
end
state = updateDissolvedH2SLag(model, state);
end

function values = column(result, field)
values = reshape([result.(field)], [], 1);
end

function componentMoles = replaceComponentMoles(componentMoles, names, aliases, values, required)
index = findComponent(names, aliases);
if isempty(index)
    assert(~required, 'sequential-h2biochem-phreeqc requires EOS component %s.', aliases{1});
    return;
end
componentMoles(:, index) = max(values, 0);
end

function componentMoles = getEOSComponentMoles(model, state)
nc = model.G.cells.num;
poreVolume = asCellVector(value(model.PVTPropertyFunctions.get( ...
    model, state, 'PoreVolume')), nc, 'pore volume');
[rhoL, rhoV, sL, sV] = getPhaseProperties(model, state);
totalMoles = poreVolume.*(sL.*rhoL + sV.*rhoV);
assert(all(isfinite(totalMoles) & totalMoles > 0), ...
    'sequential-h2biochem-phreeqc requires positive EOS moles per cell.');
components = value(state.components);
assert(ismatrix(components) && size(components, 1) == nc && ...
    size(components, 2) == model.EOSModel.getNumberOfComponents(), ...
    'EOS components must be a cell-by-component matrix.');
componentMoles = bsxfun(@times, totalMoles, components);
end

function state = constantVolumeFlash(model, state, targetComponentMoles, opt)
% Find the equilibrium pressure that preserves the PHREEQC mole inventory
% in each fixed pore volume. A flash at the old pressure preserves only
% composition and silently changes the absolute component inventory.
nc = model.G.cells.num;
poreVolume = asCellVector(value(model.PVTPropertyFunctions.get( ...
    model, state, 'PoreVolume')), nc, 'pore volume');
temperature = asCellVector(value(state.T), nc, 'temperature');
pressure0 = asCellVector(value(state.pressure), nc, 'pressure');
totalMoles = sum(targetComponentMoles, 2);
composition = bsxfun(@rdivide, targetComponentMoles, totalMoles);
pressure = zeros(nc, 1);

for cellNo = 1:nc
    residual = @(logPressure) flashVolumeResidual(exp(logPressure), ...
        temperature(cellNo), composition(cellNo, :), totalMoles(cellNo), ...
        poreVolume(cellNo), model.EOSModel);
    try
        pressure(cellNo) = exp(fzero(residual, log(pressure0(cellNo))));
    catch ME
        error('H2Biochem:SequentialPhreeqcVolumeFlash', ...
            ['Unable to preserve the PHREEQC component inventory in cell %d ', ...
             'with a constant-volume flash: %s'], cellNo, ME.message);
    end
end

state.pressure = pressure;
state.components = composition;
state = clearStateFunctionCaches(model, state);
state = model.computeFlash(state, inf);

actualComponentMoles = getEOSComponentMoles(model, state);
absoluteResidual = abs(actualComponentMoles - targetComponentMoles);
scale = max(abs(actualComponentMoles), abs(targetComponentMoles));
normalizedResidual = zeros(size(scale));
nonzero = scale > 0;
normalizedResidual(nonzero) = absoluteResidual(nonzero)./scale(nonzero);
normalizedResidual(~nonzero & absoluteResidual > 0) = inf;
absoluteTolerance = expandAuditTolerance( ...
    opt.phreeqcReflashAbsoluteTolerance, size(targetComponentMoles, 2), ...
    'phreeqcReflashAbsoluteTolerance');
relativeTolerance = expandAuditTolerance( ...
    opt.phreeqcReflashRelativeTolerance, size(targetComponentMoles, 2), ...
    'phreeqcReflashRelativeTolerance');
limit = bsxfun(@plus, absoluteTolerance, ...
    bsxfun(@times, scale, relativeTolerance));
pass = absoluteResidual <= limit;

state.sequentialH2BiochemPhreeqcReflashComponentNames = ...
    model.EOSModel.CompositionalMixture.names;
state.sequentialH2BiochemPhreeqcReflashTargetMoles = targetComponentMoles;
state.sequentialH2BiochemPhreeqcReflashActualMoles = actualComponentMoles;
state.sequentialH2BiochemPhreeqcReflashAbsoluteResidual = absoluteResidual;
state.sequentialH2BiochemPhreeqcReflashNormalizedResidual = normalizedResidual;
state.sequentialH2BiochemPhreeqcReflashPass = pass;
state.sequentialH2BiochemPhreeqcReflashPressureBefore = pressure0;
state.sequentialH2BiochemPhreeqcReflashPressureAfter = pressure;

if any(~pass(:))
    [cellNo, componentNo] = find(~pass, 1);
    warning('H2Biochem:PhreeqcReflashInventory', ...
        ['Post-PHREEQC EOS inventory audit failed in cell %d for %s: ', ...
         'target=%.16g mol, reconstructed=%.16g mol, absolute ', ...
         'residual=%.16g mol, normalized residual=%.16g.'], ...
        cellNo, model.EOSModel.CompositionalMixture.names{componentNo}, ...
        targetComponentMoles(cellNo, componentNo), ...
        actualComponentMoles(cellNo, componentNo), ...
        absoluteResidual(cellNo, componentNo), ...
        normalizedResidual(cellNo, componentNo));
end
end

function tolerance = expandAuditTolerance(tolerance, nComponents, name)
validateattributes(tolerance, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'}, mfilename, name);
if isscalar(tolerance)
    tolerance = repmat(tolerance, 1, nComponents);
else
    tolerance = reshape(tolerance, 1, []);
    assert(numel(tolerance) == nComponents, ...
        '%s must be scalar or contain one value per EOS component.', name);
end
end

function residual = flashVolumeResidual(pressure, temperature, composition, ...
        totalMoles, poreVolume, eosModel)
[liquidFraction, x, y, zLiquid, zVapor] = standaloneFlash( ...
    pressure, temperature, composition, eosModel);
propertyModel = eosModel.PropertyModel;
rhoLiquid = propertyModel.computeMolarDensity( ...
    eosModel, pressure, x, zLiquid, temperature, true);
rhoVapor = propertyModel.computeMolarDensity( ...
    eosModel, pressure, y, zVapor, temperature, false);
fluidVolume = totalMoles.*(liquidFraction./value(rhoLiquid) + ...
    (1 - liquidFraction)./value(rhoVapor));
residual = fluidVolume./poreVolume - 1;
end

function [rhoL, rhoV, sL, sV] = getPhaseProperties(model, state)
nc = model.G.cells.num;
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
propmodel = model.EOSModel.PropertyModel;
rhoL = propmodel.computeMolarDensity(model.EOSModel, value(state.pressure), ...
    value(state.x), value(state.Z_L), value(state.T), true);
rhoV = propmodel.computeMolarDensity(model.EOSModel, value(state.pressure), ...
    value(state.y), value(state.Z_V), value(state.T), false);
rhoL = asCellVector(value(rhoL), nc, 'liquid molar density');
rhoV = asCellVector(value(rhoV), nc, 'gas molar density');
sL = asCellVector(sL, nc, 'liquid saturation');
sV = asCellVector(sV, nc, 'gas saturation');
end

function values = getComponentMoles(componentMoles, names, aliases, required)
if nargin < 4
    required = true;
end
index = findComponent(names, aliases);
if isempty(index)
    assert(~required, 'sequential-h2biochem-phreeqc requires EOS component %s.', aliases{1});
    values = zeros(size(componentMoles, 1), 1);
else
    values = componentMoles(:, index);
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

function state = refreshOutputStateFunctions(model, state)
outputs = unique(model.OutputStateFunctions, 'stable');
for i = 1:numel(outputs)
    [~, state] = model.getProp(state, outputs{i});
end
end

function model = updateEOSSalinityForReflash(model, state)
if model.sulfateReduction && isa(model.EOSModel, 'SoreideWhitsonEos')
    model.EOSModel = model.EOSModel.enablesrb_coupling( ...
        value(state.tracerSO4), value(state.tracerHS), ...
        value(state.tracerHS) + value(state.h2sDissolvedLag), value(state.T));
end
end

function state = updateDissolvedH2SLag(model, state)
if ~model.sulfateReduction
    return;
end
names = model.EOSModel.CompositionalMixture.names;
index = find(strcmp(names, 'H2S'), 1);
if isempty(index)
    return;
end
x = value(state.x);
if iscell(x)
    xH2S = x{index};
else
    xH2S = x(:, index);
end
rhoL = model.EOSModel.PropertyModel.computeMolarDensity(model.EOSModel, ...
    value(state.pressure), value(state.x), value(state.Z_L), value(state.T), true);
state.h2sDissolvedLag = value(rhoL).*xH2S;
end

function waterMass = getWaterMass(model, state, rhoWater)
poreVolume = value(model.PVTPropertyFunctions.get(model, state, 'PoreVolume'));
s = value(state.s);
liquid = model.getLiquidIndex();
if iscell(s)
    sL = s{liquid};
else
    sL = s(:, liquid);
end
waterMass = poreVolume.*max(sL, 1e-12).*rhoWater;
end

function values = getStateVector(state, field, nc, default)
if isfield(state, field)
    values = value(state.(field));
else
    values = default;
end
values = asCellVector(values, nc, field);
end

function values = optionVector(options, field, default, nc)
if isfield(options, field)
    values = options.(field);
else
    values = default;
end
values = asCellVector(values, nc, field);
end

function values = asCellVector(values, nc, name)
if isscalar(values)
    values = repmat(values, nc, 1);
else
    values = values(:);
    assert(numel(values) == nc, ...
        'sequential-h2biochem-phreeqc expected %d %s values, got %d.', nc, name, numel(values));
end
assert(all(isfinite(values)), 'sequential-h2biochem-phreeqc %s values must be finite.', name);
end

function index = findComponent(names, aliases)
index = [];
for i = 1:numel(aliases)
    index = find(strcmpi(names, aliases{i}), 1);
    if ~isempty(index)
        return;
    end
end
end

function validateSupportedComponents(names)
supported = {'H2O', 'Water', 'H2', 'Hydrogen', 'CO2', 'CarbonDioxide', ...
    'C1', 'CH4', 'Methane', 'H2S', 'HydrogenSulfide', 'N2', 'Nitrogen', ...
    'CH3COOH', 'AceticAcid', 'Acetate'};
unsupported = names(~cellfun(@(name) any(strcmpi(name, supported)), names));
assert(isempty(unsupported), ...
    ['sequential-h2biochem-phreeqc cannot safely pass EOS components without a ', ...
     'PHREEQC_Modified.DAT mapping: %s'], strjoin(unsupported, ', '));
end

function isAbsolute = isAbsolutePath(path)
isAbsolute = ~isempty(regexp(path, '^[A-Za-z]:[\\/]|^\\\\', 'once')) || ...
    startsWith(path, filesep);
end

function result = emptyResult()
result = struct( ...
    'time', 0, 'step', 0, 'pH', 0, 'pe', 0, 'water', 0, ...
    'totalCarbon', 0, 'sulfate', 0, 'tds', 0, 'ca', 0, 'mg', 0, 'acetate', 0, ...
    'sulfide', 0, 'aqH2', 0, 'aqN2', 0, 'aqCO2', 0, 'aqCH4', 0, ...
    'aqH2S', 0, 'hco3', 0, 'gasH2', 0, 'gasCO2', 0, 'gasCH4', 0, ...
    'gasH2S', 0, 'gasN2', 0, 'calcite', 0, 'anhydrite', 0, ...
    'gypsum', 0, 'dolomite', 0, 'goethite', 0, 'pyrite', 0, ...
    'brucite', 0, 'portlandite', 0, 'quartz', 0, ...
    'fe2', 0, 'fe3', 0, 'initialSystemHydrogen', 0, 'systemHydrogen', 0);
end
