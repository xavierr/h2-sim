%% 3D anticlinal-dome aquifer (coarsened), compositional H2 storage, 3 reactions
% =========================================================================
% Three-reaction (MET, ACE, SRB) compositional bio-reactive simulation on
% the 3D anticlinal-dome aquifer of the h2store module, using the same
% black-oil -> compositional conversion as
% exampleCompositionalBacterial2DThreeReactions.m.
%
% The full dome grid (51x51x32 ~ 83k cells) is far too large for a
% compositional bio-reactive run, so this script (only -- the h2store setup
% is untouched) builds a fresh coarse Cartesian grid over the same bounding
% box, reproduces the anticline with a Gaussian bump, and samples the
% layered rock from the fine model by nearest neighbour. Adjust
% `coarsenFactor` for resolution.
%
% Black-oil setup:   modelForDome3DAquifer / simulateStorageAnticlinalDomeAquifer
%   Ahmed, E., et al. (2024), Adv. Water Resour. 191, 104772.
% Three-reaction model: Shojaee et al. (2025), compositional PHREEQC benchmark.
%
% Single simulation: compositional flow + microbial growth/decay only.
% No molecular diffusion, dispersion, microbial diffusion, chemotaxis or
% bio-clogging.
%
%   MET : 4 H2 +  CO2      -> CH4     + 2 H2O
%   ACE : 4 H2 + 2 CO2     -> CH3COOH + 2 H2O
%   SRB : 4 H2 +  SO4^2-   -> H2S     + 2 H2O    (SO4/HS are aqueous tracers)
% =========================================================================

clearvars;
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui vemmech ...
    test-suite compositional h2store h2-biochem

coarsenFactor = [4, 4, 4];   % 51x51x32 -> 13x13x8 (~1350 cells)
caseName = 'H2_3D_DOME_3RXN_COARSE';

%% -------------------------------------------------- Black-oil dome case --
% RS-only deck: the compositional conversion cannot handle the PVTG/RV of
% DOME_RSRV.DATA.
dataFile = fullfile(getDatasetPath('anticlinal_dome'), 'DOME_RS.DATA');
deck = readEclipseDeck(dataFile);
[~, options, state0Bo, modelBo, scheduleBo] = modelForDome3DAquifer(deck, ...
    'numCycles', 1);
compModelFine = convertBlackOilModelToCompositionalModel(modelBo);

%% --------------------------------------------------- Coarsen the grid ---
% A fresh coarse Cartesian grid over the fine bounding box, with the
% anticline reproduced by a Gaussian node-z bump, and the layered rock
% sampled from the fine model by nearest neighbour. This keeps a normal
% grid (with nodes) so addWell / plotGrid work.
Gf = modelBo.G;
if isfield(Gf, 'cartDims') && numel(Gf.cartDims) == 3 && all(Gf.cartDims > 1)
    cdim = max(2, round(Gf.cartDims ./ coarsenFactor));
else
    cdim = [13, 13, 8];   % fallback if the fine grid carries no cartDims
end
lo = min(Gf.nodes.coords);
hi = max(Gf.nodes.coords);
G = cartGrid(cdim, hi - lo);
G.nodes.coords = G.nodes.coords + lo;
xc = G.nodes.coords(:, 1);  yc = G.nodes.coords(:, 2);
xm = mean([lo(1), hi(1)]);   ym = mean([lo(2), hi(2)]);
zf = (G.nodes.coords(:, 3) - lo(3)) ./ max(hi(3) - lo(3), eps);   % 0 top .. 1 base
bump = 100 .* exp(-((xc - xm).^2 + (yc - ym).^2) ./ (2*500^2));
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) - bump .* (1 - zf);   % anticline
G = computeGeometry(G);
nc = G.cells.num;

% nearest fine cell for each coarse cell (looped to avoid a huge Delaunay)
idx = zeros(nc, 1);
for c = 1:nc
    d2 = sum((Gf.cells.centroids - G.cells.centroids(c, :)).^2, 2);
    [~, idx(c)] = min(d2);
end
rock = struct('poro', modelBo.rock.poro(idx), ...
    'perm', max(modelBo.rock.perm(idx, :), 1e-6*milli*darcy));
fprintf('Coarsened %d -> %d cells (%dx%dx%d).\n', Gf.cells.num, nc, ...
    cdim(1), cdim(2), cdim(3));

%% --------------------------------------------------------------- Fluid ---
compFluid = TableCompositionalMixture( ...
    {'Water', 'Hydrogen', 'CarbonDioxide', 'Methane', ...
     'HydrogenSulfide', 'AceticAcid'}, ...
    {'H2O',   'H2',       'CO2',           'C1',      'H2S', 'CH3COOH'});
reactNames   = {'MethanogenicArchae', 'AcetogenicBacteria', 'SulfateReducingBacteria'};
biomassNames = {'bactM', 'bactA', 'bactS'};
biochemFluid = TableBioChemMixture(reactNames, biomassNames);   % default kinetics

T0 = options.tempCharge;              % 323.15 K
initialSO4 = 0.1;
EOS = SoreideWhitsonEos(G, compFluid, ...
    'msalt', 0, 'pH', 7.2, 'initial_NaCl', 0, ...
    'initial_SO4', initialSO4, 'rho_water', 1000);

%% ------------------------------------------------------------- Model ----
backend = DiagonalAutoDiffBackend('modifyOperators', true);
model = BiochemistryModel(G, rock, compModelFine.fluid, ...
    compFluid, biochemFluid, false, backend, ...
    'oil', true, 'gas', true, ...
    'bacteriamodel', true, ...
    'bactDiffusion', false, ...
    'chemotaxisEffect', false, ...
    'molecularDiffusion', false, ...
    'molecularDispersion', false, ...
    'liquidPhase', 'O', 'vaporPhase', 'G');
model.EOSModel = EOS;
model.gravity = modelBo.gravity;
model.outputFluxes = false;
model.OutputStateFunctions{end+1} = 'ComponentPhaseDensity';

%% --------------------------------------------------------- Initial state --
nbact0 = [15, 15, 15];
comp0  = [0.8480, 1.0e-5, 1.0e-5, 0.1530, 0.0, 0.0];   % H2O H2 CO2 C1 H2S Ac
% Sample the black-oil hydrostatic initial state onto the coarse grid.
p0 = state0Bo.pressure(idx);
s0 = state0Bo.s(idx, :);
state0 = initCompositionalStateBacteria(model, p0, T0, s0, comp0, nbact0, EOS);
if isa(model.EOSModel, 'SoreideWhitsonEos') && model.sulfateReduction
    state0.tracerSO4       = repmat(initialSO4*1000, nc, 1);   % kg/m3
    state0.tracerHS        = zeros(nc, 1);
    state0.h2sDissolvedLag = zeros(nc, 1);
end

%% ------------------------- Remap wells + schedule to the coarse grid ----
% Single vertical well at the crest of the dome, fully penetrating.
wi = max(1, ceil(cdim(1)/2));
wj = max(1, ceil(cdim(2)/2));
schedule = scheduleBo;
for i = 1:numel(schedule.control)
    Wf = schedule.control(i).W;
    W  = verticalWell([], G, rock, wi, wj, 1:cdim(3), ...
        'Name', Wf.name, 'Type', Wf.type, 'Val', Wf.val, 'sign', Wf.sign, ...
        'comp_i', [0, 1], 'Radius', 0.5);
    if strcmpi(Wf.name, 'cushion')
        W.components = [0.0, 0.10, 0.90, 0.0, 0.0, 0.0];   % CO2-rich cushion
    else
        W.components = [0.0, 0.95, 0.05, 0.0, 0.0, 0.0];   % H2-rich charge
    end
    W.T = T0;
    schedule.control(i).W  = W;
    schedule.control(i).bc = [];
end

% Subdivide every step to <= 3 days and ramp the first control block up.
maxDt = 3*day;
dt = schedule.step.val(:);  ctl = schedule.step.control(:);
nSub = max(1, ceil(dt./maxDt));
newDt  = cell2mat(arrayfun(@(k) repmat(dt(k)/nSub(k), nSub(k), 1), ...
    (1:numel(dt)).', 'UniformOutput', false));
newCtl = cell2mat(arrayfun(@(k) repmat(ctl(k), nSub(k), 1), ...
    (1:numel(dt)).', 'UniformOutput', false));
first = newCtl(1);  i1 = find(newCtl == first);
ramp  = rampupTimesteps(sum(newDt(i1)), maxDt, 8);
schedule.step.val     = [ramp; newDt(i1(end)+1:end)];
schedule.step.control = [repmat(first, numel(ramp), 1); newCtl(i1(end)+1:end)];

%% --------------------------------------------- Rock-property figure -----
yCut = G.cells.centroids(:, 2) > mean(G.cells.centroids(:, 2));
fig = paperFigure([24, 10], '3D dome aquifer (coarse) - rock properties');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout); axes(ax); %#ok<LAXES>
plotCellData(G, rock.poro, yCut, 'EdgeAlpha', 0.2); view(3); axis tight
title(ax, 'Porosity'); colorbar(ax);
ax = nexttile(layout); axes(ax); %#ok<LAXES>
plotCellData(G, log10(rock.perm(:,1)/(milli*darcy)), yCut, 'EdgeAlpha', 0.2);
view(3); axis tight
title(ax, 'Permeability (log_{10} mD)'); colorbar(ax);
title(layout, '3D anticlinal-dome aquifer (coarsened), three reactions');
paperExport(fig, 'threeReactions_3D_rock_properties');

%% ---------------------------------------------------------- Simulation --
[nls, ~] = setupOptimizedLinearSolver(model, ...
    'complexityLevel', 'high', 'solverTolerance', 1e-3, ...
    'maxNonlinIter', 12, 'cprDamp', []);
nls.maxTimestepCuts         = 15;
nls.useRelaxation           = true;
nls.enforceResidualDecrease = false;
problem = packSimulationProblem(state0, model, schedule, caseName, ...
    'NonLinearSolver', nls);
simulatePackedProblem(problem);
[ws, states] = getPackedSimulatorOutput(problem);
assert(all(~cellfun(@isempty, states)), ...
    'The simulation did not complete all schedule steps.');

%% ----------------------------------------------------------- Results ----
timeDays = cumsum(schedule.step.val(:))./day;
nR = model.biochemFluid.nbioreact;
cumMol = zeros(numel(states), nR);
for r = 1:nR
    [~, perCell] = computeH2Consumption(states, schedule, model, r);
    cumMol(:, r) = sum(perCell, 1).';
end
fprintf('\nCumulative microbial H2 consumption (mol):\n');
fprintf('  MET %.3e   ACE %.3e   SRB %.3e   total %.3e\n', ...
    cumMol(end, 1), cumMol(end, 2), cumMol(end, 3), sum(cumMol(end, :)));

fig = paperFigure([20, 11], '3D three-reaction H2 consumption over time');
ax = axes(fig); hold(ax, 'on');
lab = {'MET', 'ACE', 'SRB'};
col = paperColors(nR + 1);
for r = 1:nR
    plot(ax, timeDays, cumMol(:, r), 'LineWidth', 1.8, ...
        'Color', col(r, :), 'DisplayName', lab{r});
end
plot(ax, timeDays, sum(cumMol, 2), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'total');
xlabel(ax, 'Time (days)');
ylabel(ax, 'Cumulative H_2 consumed (mol)');
title(ax, '3D anticlinal-dome benchmark (coarse), three reactions');
legend(ax, 'Location', 'northwest');
styleAxes(ax);
paperExport(fig, 'threeReactions_3D_H2_consumption_over_time');

%% ---------------------------------------------- Interactive explorers ---
figure('Name', 'States - three-reaction 3D dome (coarse)');
plotToolbar(model.G, states); view(3); axis tight
figure('Name', 'Well solutions');
plotWellSols(ws);

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%}
