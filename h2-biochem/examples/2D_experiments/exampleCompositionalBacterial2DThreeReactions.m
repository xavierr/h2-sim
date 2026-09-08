%% 2D Compositional Hydrogen Storage with Three Reactions (MET, ACE, SRB)
% =========================================================================
% This example simulates hydrogen storage in a 2D dome-shaped saline aquifer.
% It includes three microbial reactions:
%   - Methanogenesis (MET):  4 H2 + CO2  -> CH4 + 2 H2O
%   - Acetogenesis (ACE):    4 H2 + 2 CO2 -> CH3COOH + 2 H2O
%   - Sulfate Reduction (SRB): 4 H2 + SO4^2- -> H2S + 2 H2O
%
% All kinetic parameters are taken from the default database
% (bioChemFluidsStructs.m) – no overrides.
%
% Sulfate (SO4) and bisulfide (HS) are treated as aqueous tracers.
% H2S is a volatile EOS component.
%
% Scenarios:
%   1) With bacteria and bio-clogging
%   2) With bacteria but without clogging
%   3) Abiotic (no bacteria)
%
% References:
%   - Original 2D case: Ahmed et al., 2024
%   - Three‑reaction model: Shojaee et al., 2025 (compositional PHREEQC)
% =========================================================================

clearvars;
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr test-suite spe10
mrstModule add compositional h2-biochem h2store

%% Define case identifiers and Eclipse deck
baseName = 'H2_STORAGE_DOME_TRAP_3RXN';
dataPath = getDatasetPath('h2storage');
dataFile = fullfile(dataPath, 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

%% Warn about computational cost
warning('ComputationalCost:High', ...
    'This is a multiple-cycle example; consider reducing cycles for faster runs.');

%% Set up black-oil model and schedule
[~, ~, state0Bo, modelBo, scheduleBo, ~] = modelForSimple2DAquifer(deck, 'numcycles', 10);

%% Convert black-oil to compositional model
model = convertBlackOilModelToCompositionalModel(modelBo);
state0 = convertBlackOilStateToCompositional(modelBo, state0Bo);

%% Define compositional fluid with 6 components (including H2S and AceticAcid)
compFluid = TableCompositionalMixture(...
    {'Water', 'Hydrogen', 'CarbonDioxide', 'Methane', 'HydrogenSulfide', 'AceticAcid'}, ...
    {'H2O',   'H2',       'CO2',           'C1',      'H2S',             'CH3COOH'});

%% Define biochemical reactions – three reactions, default parameters
reactNames = {'MethanogenicArchae', 'AcetogenicBacteria', 'SulfateReducingBacteria'};
biomassNames = {'bactM', 'bactA', 'bactS'};
biochemFluid = TableBioChemMixture(reactNames, biomassNames);
% No overrides – uses default values from bioChemFluidsStructs.

%% EOS – Soreide‑Whitson
initialSO4 = 0.1;
EOS = SoreideWhitsonEos(model.G, compFluid, ...
    'msalt', 0, ...
    'pH', 7.2, ...
    'initial_NaCl', 0, ...
    'initial_SO4', initialSO4, ...
    'rho_water', 1000);
model.EOSModel = EOS;

nc = model.G.cells.num;
T0 = 273.15 + 44.35;  % Initial temperature (from original deck)

%% Initial global composition (from original deck: 84.8% H2O, 15.3% CH4, trace CO2/H2)
% We need to adjust for the additional components (H2S and AceticAcid are zero initially)
comp0 = repmat([0.8480, 1.0e-5, 1.0e-5, 0.1530, 0.0, 0.0], nc, 1);

%% Bio‑clogging parameters
bacteriamodel = true;

%% Setup BiochemistryModel with three reactions
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
arg = {model.G, model.rock, model.fluid, compFluid, biochemFluid, ...
    false, diagonal_backend, 'oil', true, 'gas', true, ...
    'bacteriamodel', bacteriamodel, ...
    'bactDiffusion', false, ...
    'chemotaxisEffect', false, ...
    'molecularDiffusion', true, ...
    'molecularDispersion', true, ...,
    'liquidPhase', 'O', ...
    'vaporPhase', 'G'};
model = BiochemistryModel(arg{:});
model.OutputStateFunctions{end+1} = 'ComponentPhaseDensity';
model.gravity = modelBo.gravity;

clogModel = true;  % will be toggled per scenario
nbact0 = [15 15 15]; % Initial bacteria (normalized)
nc_bact = [120, 120, 120];
cp = [1.0, 1.0, 1.0];
modelWithClog = setupBioCloggingModel(model, nbact0, nc_bact, cp, clogModel);

%% Initialize compositional state (with tracers for sulfate)
state0 = initCompositionalStateBacteria( ...
    modelWithClog, state0.pressure, T0, state0.s, comp0, nbact0, EOS);

% Add sulfate and bisulfide tracers (if SRB is active)
% Note: modelWithClog.sulfateReduction will be set true because SRB reaction exists.
if isa(modelWithClog.EOSModel, 'SoreideWhitsonEos') && modelWithClog.sulfateReduction
    rho_water = 1000;    % kg/m3
    state0.tracerSO4 = repmat(initialSO4 * rho_water, nc, 1);
    state0.tracerHS  = zeros(nc, 1);
    state0.h2sDissolvedLag = zeros(nc, 1);
end

%% Update schedule controls for compositional injection
schedule = scheduleBo;
for i = 1:numel(schedule.control)
    schedule.control(i).W.compi = [0, 1]; % Well components
    if strcmp(schedule.control(i).W.name, 'cushion') && i<11
        % Cushion gas injection: 10% H2, 90% CO2 (as in original)
        schedule.control(i).W.components = [0.0, 0.1, 0.9, 0.0, 0.0, 0.0];
    else
        % Main H2 injection: 95% H2, 5% CO2
        schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0, 0.0, 0.0];
    end
    schedule.control(i).W.T = T0; % Temperature
    schedule.control(i).bc = [];   % Remove BCs from controls
end
modelWithClog.outputFluxes = false;

%% Plot initial grid, porosity, and permeability
fig = paperFigure([24, 10], '2D dome aquifer - rock properties');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout); axes(ax); %#ok<LAXES>
plotCellData(modelWithClog.G, model.rock.poro); hold(ax, 'on');
plotGrid(modelWithClog.G, schedule.control(1).W.cells, ...
    'FaceColor', 'red', 'LineStyle', 'none');
title(ax, 'Porosity'); axis(ax, 'off', 'tight'); colorbar(ax);
ax = nexttile(layout); axes(ax); %#ok<LAXES>
plotCellData(modelWithClog.G, log10(model.rock.perm(:,1))); hold(ax, 'on');
plotGrid(modelWithClog.G, schedule.control(1).W.cells, ...
    'FaceColor', 'red', 'LineStyle', 'none');
title(ax, 'Permeability (log_{10} mD)'); axis(ax, 'off', 'tight'); colorbar(ax);
title(layout, '2D dome-shaped aquifer with three reactions');
paperExport(fig, 'threeReactions_2D_rock_properties');

%% Initialize solvers
[nls,~] = setupOptimizedLinearSolver(modelWithClog, 'complexityLevel', 'high', ...
    'solverTolerance', 1e-3, ...
    'maxNonlinIter', 10, ...
    'cprDamp', []);

%% Pack and run simulations
% --- Scenario 1: With bacteria and bio-clogging
caseNameWithClogging = [baseName '_WITH_CLOGGING'];
% Enable clogging
problemWithClogging = packSimulationProblem(state0, modelWithClog, schedule, caseNameWithClogging, 'NonLinearSolver', nls);
simulatePackedProblem(problemWithClogging);

% --- Scenario 2: With bacteria but without clogging
caseNameNoClogging = [baseName '_NO_CLOGGING_DIFF__DISP'];
modelNoClog  = setupBioCloggingModel(model, nbact0, nc_bact, cp, false);
[nlsNoClog, ~] = setupOptimizedLinearSolver(modelNoClog, ...
    'complexityLevel', 'high', ...
    'solverTolerance', 1e-3, ...
    'maxNonlinIter', 10, ...
    'cprDamp', []);
problemNoClogging = packSimulationProblem( ...
    state0, modelNoClog, schedule, caseNameNoClogging, ...
    'NonLinearSolver', nlsNoClog);
simulatePackedProblem(problemNoClogging);

% --- Scenario 3: Abiotic (no bacteria)
caseNameNoBact = [baseName '_NO_BACT_DIFF__DISP'];
modelNoBact = modelNoClog;
modelNoBact.bacteriamodel = false;
state0NoBact = initCompositionalState(modelNoBact, state0.pressure, T0, state0.s, comp0, EOS);
% The abiotic model has a different primary-variable layout.
[nlsNoBact, ~] = setupOptimizedLinearSolver(modelNoBact, 'complexityLevel', 'high', ...
    'solverTolerance', 1e-3, ...
    'maxNonlinIter', 10, ...
    'cprDamp', []);
problemNoBact = packSimulationProblem(state0NoBact, modelNoBact, schedule, caseNameNoBact, 'NonLinearSolver', nlsNoBact);
simulatePackedProblem(problemNoBact);

%% Get and compare results
[wsWithClog, statesWithClog] = getPackedSimulatorOutput(problemWithClogging);
[wsNoClog, statesNoClog] = getPackedSimulatorOutput(problemNoClogging);
[wsNoBact, statesNoBact] = getPackedSimulatorOutput(problemNoBact);

% Interactive state explorers (plotToolbar / plotWellSols are GUIs and are
% not exported by paperExport).
figure('Name', 'States - with clogging');    plotToolbar(modelWithClog.G, statesWithClog);
figure('Name', 'States - no clogging');      plotToolbar(modelWithClog.G, statesNoClog);
figure('Name', 'States - abiotic');          plotToolbar(modelWithClog.G, statesNoBact);
figure('Name', 'Well solutions');
plotWellSols({wsWithClog, wsNoClog, wsNoBact}, ...
    'datasetnames', {'with clogging', 'no clogging', 'abiotic'});

%% Copyright notice
% <html>
% <p><font size="-1">
% Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% ...
% </html>