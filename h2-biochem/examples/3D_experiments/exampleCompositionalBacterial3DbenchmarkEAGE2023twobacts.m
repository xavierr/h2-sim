%% MRST Simulation for Hydrogen Storage with Two-Species Bacterial Growth
%
% Description:
% Two-reaction (methanogenic archaea + acetogenic bacteria) extension of
% exampleCompositionalBacterial3DbenchmarkEAGE2023.m: same 3D dome-shaped
% domain, with and without microbial activity, based on:
%
% Khoshnevis, N., Hogeweg, S., Goncalves Machado, C., and Hagemann, B.
% "Numerical Modeling of Bio-Reactive Transport During Underground
% Hydrogen Storage – a Benchmark Study", EAGE, vol. 1, pp. 1-5, 2023.
% https://doi.org/10.3997/2214-4609.202321087
%
% Features:
% - Two-phase system (liquid 'O' and gas 'G')
% - Five components (H2O, H2, CO2, CH4, CH3COOH)
% - With/without microbial activity of MethanogenicArchae + AcetogenicBacteria
% - Cyclic injection/production schedule
% - Base compositional flow + bacterial growth/decay only: no molecular
%   diffusion, mechanical dispersion, or microbial diffusion are enabled
%   at this 3D scale (1800 cells) -- those transport extras are exercised
%   separately on the 1D benchmark
%   (examples/Sensitivity_analysis/runTransportEffectsSensitivity.m).
%
% Author: Stéphanie Delage Santacreu
% Date: 23/09/2025
% Organization: Université de Pau et des Pays de l'Adour, E2S UPPA,
%               CNRS, LFCR, UMR5150, Pau, France
% -------------------------------------------------------------------------

%% Initialization
clear; clc;
mrstModule add h2-biochem compositional ad-blackoil ad-core ad-props mrst-gui

gravity reset on

%% ============ Grid and Rock Properties =====================
[nx, ny, nz] = deal(15, 15, 8);
[Lx, Ly, Lz] = deal(1525, 1525, 50);

G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 1000; % Reservoir depth
G = computeGeometry(G);

perm = [100, 100, 10] .* milli*darcy;
rock = makeRock(G, perm, 0.2);

%% Fluid Properties
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane', 'AceticAcid'}, ...
    {'H2O', 'H2', 'CO2', 'C1', 'CH3COOH'});
biochemFluid = TableBioChemMixture({'MethanogenicArchae','AcetogenicBacteria'},{'bactM','bactA'});


[rhow, rhog]   = deal(999.7 * kilogram/meter^3, 1.2243 * kilogram/meter^3);
[viscow, viscog] = deal(1.3059 * centi*poise, 0.01763 * centi*poise);
[cfw, cfg]     = deal(5.0e-5/barsa, 1.0/barsa);

[srw, src] = deal(0.2, 0.05);
P0 = 100 * barsa;

fluid = initSimpleADIFluid('phases', 'OG', ...
    'mu', [viscow, viscog], ...
    'rho', [rhow, rhog], ...
    'pRef', P0, ...
    'c', [cfw, cfg], ...
    'n', [2, 2], ...
    'smin', [srw, src]);

Pe = 0.1 * barsa;
fluid.pcOG = @(sg) Pe * max((1 - sg - srw) ./ (1 - srw), 1e-5).^(-1/2);

%% Wells
W = [];
n1 = floor(0.5*nx) + 1;
n2 = floor(0.5*ny) + 1;

% Injector (CO2-rich)
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'comp_i', [0, 1], 'Radius', 0.5, 'name', 'Injector', ...
    'type', 'rate', 'Val', 1e6*meter^3/day, 'sign', 1);
W(end).components = [0.0, 0.6, 0.4, 0.0, 0.0];

% Rest (shut-in)
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'compi', [0, 1], 'Radius', 0.5, 'name', 'Rest', ...
    'type', 'rate', 'Val', 0.0, 'sign', 1);
W(end).components = [0.0, 0.95, 0.05, 0.0, 0.0];

% Injector (H2-rich)
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'comp_i', [0, 1], 'Radius', 0.5, 'name', 'InjectorH2', ...
    'type', 'rate', 'Val', 1e6*meter^3/day, 'sign', 1);
W(end).components = [0.0, 0.95, 0.05, 0.0, 0.0];

% Producer
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'compi', [0, 1], 'Radius', 0.5, 'name', 'Producer', ...
    'type', 'rate', 'Val', -1e6*meter^3/day, 'sign', -1);
W(end).components = [0.0, 0.95, 0.05, 0.0, 0.0];

%% Time Schedule
ncycles    = 1;
deltaT     = 5*day;
nbj_buildUp = 60*day; nbj_rest = 20*day;
nbj_inject  = 30*day; nbj_idle = 20*day;
nbj_prod    = 30*day; nbj_idle1 = 20*day;

[schedule, TotalTime, nbuildUp, nrest, ninject, nidle, nprod, nidle1] = ...
    createCyclicScenario(deltaT, ncycles, nbj_buildUp, nbj_rest, ...
    nbj_inject, nbj_idle, nbj_prod, nbj_idle1, W);
% smaller step for initialization
schedule.step.val(1) = 1*day;

%% EOS and Model
compEOS = SoreideWhitsonEos(G, compFluid);

backend = DiagonalAutoDiffBackend('modifyOperators', true);

%% Initial State
T0 = 40 + 273.15;
s0 = [0.2, 0.8];
z0 = [0.7, 0.0, 0.02, 0.28, 0.0];

%% --- Simulation 1: Without bacteria ---
arg = {G, rock, fluid, compFluid, biochemFluid, true, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', false, ...
    'liquidPhase', 'O', 'vaporPhase', 'G'};

model_nobact = BiochemistryModel(arg{:});
model_nobact.outputFluxes = false;
model_nobact.EOSModel = compEOS;

state0_nobact = initCompositionalState(model_nobact, P0, T0, s0, z0);

lsolve = selectLinearSolverAD(model_nobact);
nls = NonLinearSolver(); nls.LinearSolver = lsolve;

problem_nobact = packSimulationProblem(state0_nobact, model_nobact, schedule, ...
    'Benchmark_NoBacteria3', 'NonLinearSolver', nls);

%% --- Simulation 2: With bacteria ---
model_bact = BiochemistryModel(G, rock, fluid, compFluid, biochemFluid, true, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', true, ...
    'liquidPhase', 'O', 'vaporPhase', 'G');
model_bact.outputFluxes = false;
model_bact.EOSModel = compEOS;
model_bact.bact_capProp=1.e-3;

nbioreact=model_bact.biochemFluid.nbioreact;
nbact0 = ones(1, nbioreact);
nbactMax = ones(1, nbioreact) * 1e8;
model_bact.biochemFluid.nbactMax = nbactMax;
state0_bact = initCompositionalStateBacteria(model_bact, P0, T0, s0, z0, nbact0, compEOS);

lsolve = selectLinearSolverAD(model_bact, 'useAMGCLCPR', true);
nls.LinearSolver = lsolve;

problem_bact = packSimulationProblem(state0_bact, model_bact, schedule, ...
    'Benchmark_Bacteria2', 'NonLinearSolver', nls);
simulatePackedProblem(problem_bact,'restartStep', 1);
[ws_bact, states_bact] = getPackedSimulatorOutput(problem_bact);
results_bact = postProcessResults(states_bact, ws_bact, model_bact, 'bact');

%% --- Efficiency and Comparison ---
eff_nobact = calculateH2Efficiency(results_nobact.H2_well, ...
    nbuildUp, nrest, ninject, nidle, nprod, ncycles,nidle1);
fprintf('H2 Production Efficiency (No bacteria): %.2f%%\n', eff_nobact);

eff_bact = calculateH2Efficiency(results_bact.H2_well, ...
    nbuildUp, nrest, ninject, nidle, nprod, ncycles,nidle1);
fprintf('H2 Production Efficiency (With bacteria): %.2f%%\n', eff_bact);

% Component changes
H2_loss = (abs(results_nobact.totMassH2 - results_bact.totMassH2) ./ results_nobact.totMassH2) * 100;
CO2_loss = (abs(results_nobact.totMassCO2 - results_bact.totMassCO2) ./ results_nobact.totMassCO2) * 100;
C1_gain  = (abs(results_nobact.totMassC1 - results_bact.totMassC1) ./ results_nobact.totMassC1) * 100;

fprintf('Total H2 loss due to bacteria: %.2f%%\n', H2_loss(end));
fprintf('Total CO2 loss due to bacteria: %.2f%%\n', CO2_loss(end));
fprintf('Total C1 gain due to bacteria:  %.2f%%\n', C1_gain(end));

%% --- Plot Benchmark Results ---
% Choose plotting function based on whether we have only MethanogenicArchae
hasOnlyMethanogenic = all(strcmp(model_bact.biochemFluid.metabolicReaction, 'MethanogenicArchae'));
if nbioreact==1 || hasOnlyMethanogenic
    plotBenchmarckAEGE2023(numel(states_nobact), ...
        results_nobact.pressure, results_bact.pressure, ...
        H2_loss, results_nobact.totMassH2, results_bact.totMassH2, ...
        G, states_bact, nbact0);
else
    plotBenchmarckAEGE2023MetAcet(numel(states_nobact), ...
        results_nobact.pressure, results_bact.pressure, ...
        H2_loss, results_nobact.totMassH2, results_bact.totMassH2, ...
        results_bact.totMassCH3COOH,G, states_bact, nbact0);
end
%% Post-processing Function
function results = postProcessResults(states, ws, model, caseType)
% Post-process simulation results
%
% Parameters:
%   states - Cell array of state variables
%   ws - Well solution data
%   model - Simulation model
%   caseType - 'bact' or 'nobact'

% Get component indices
namecp = model.EOSModel.getComponentNames();
ncomp = model.EOSModel.getNumberOfComponents();
nbioreact=model.biochemFluid.nbioreact;
nT = numel(states);

indH2 = find(strcmp(namecp, 'H2'));
indCO2 = find(strcmp(namecp, 'CO2'));
indC1 = find(strcmp(namecp, 'C1'));
indCH3COOH = find(strcmp(namecp, 'CH3COOH'));

% Initialize result structure with fields valid for all configurations
results = struct();
baseFields = {'xH2', 'yH2', 'xCO2', 'yCO2', 'pressure', 'H2_well', 'CO2_well', ...
    'totMassH2', 'totMassCO2', 'FractionMassH2', 'FractionMassCO2', 'totMassComp'};

% Add reactor-specific fields
resultFields = baseFields;
for i = 1:nbioreact
    bactName = model.biochemFluid.bactnames{i};
    % Get product component for this reactor
    prodComp = model.biochemFluid.p2{i};

    if ~any(strcmp(resultFields, ['y', prodComp]))
        resultFields = [resultFields, {['y', prodComp]}];
    end
    if ~any(strcmp(resultFields, ['totMass', prodComp]))
        resultFields = [resultFields, {['totMass', prodComp], ['FractionMass', prodComp]}];
    end
    if ~any(strcmp(resultFields, [prodComp, '_well']))
        resultFields = [resultFields, {[prodComp, '_well']}];
    end
end

for i = 1:numel(resultFields)
    results.(resultFields{i}) = zeros(nT, 1);
end

% Process results for each time step
for i = 1:nT
    results.xH2(i) = max(states{i}.x(:, indH2));
    results.yH2(i) = max(states{i}.y(:, indH2));
    results.xCO2(i) = max(states{i}.x(:, indCO2));
    results.yCO2(i) = max(states{i}.y(:, indCO2));

    results.pressure(i) = mean(states{i}.pressure(:));
    results.H2_well(i) = ws{i}.H2;
    results.CO2_well(i) = ws{i}.CO2;

    % Calculate mass balances
    for j = 1:ncomp
        results.totMassComp(i) = results.totMassComp(i) + ...
            sum(states{i}.FlowProps.ComponentTotalMass{j});
    end

    results.totMassH2(i) = sum(states{i}.FlowProps.ComponentTotalMass{indH2});
    results.totMassCO2(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCO2});
    results.FractionMassH2(i) = results.totMassH2(i) / results.totMassComp(i);
    results.FractionMassCO2(i) = results.totMassCO2(i) / results.totMassComp(i);

end

% Process reactor-specific results
for i = 1:nT
    for j = 1:nbioreact
        % Get product component for this reactor
        prodComp = model.biochemFluid.p2{j};
        prodIdx = find(strcmp(namecp, prodComp));

        % Store phase mole fractions - try to store both phases if available
        if ~isempty(prodIdx)
            % Try liquid phase
            resultsFieldLiq = ['x', prodComp];
            if isfield(results, resultsFieldLiq) && size(states{i}.x, 2) >= prodIdx
                results.(resultsFieldLiq)(i) = max(states{i}.x(:, prodIdx));
            end
            % Try gas phase
            resultsFieldGas = ['y', prodComp];
            if isfield(results, resultsFieldGas) && size(states{i}.y, 2) >= prodIdx
                results.(resultsFieldGas)(i) = max(states{i}.y(:, prodIdx));
            end
        end

        % Store well flux
        wellField = [prodComp, '_well'];
        if isfield(ws{i}, prodComp) && isfield(results, wellField)
            results.(wellField)(i) = ws{i}.(prodComp);
        end

        % Store mass balance
        massField = ['totMass', prodComp];
        fracField = ['FractionMass', prodComp];
        if isfield(results, massField)
            results.(massField)(i) = sum(states{i}.FlowProps.ComponentTotalMass{prodIdx});
            if isfield(results, fracField)
                results.(fracField)(i) = results.(massField)(i) / results.totMassComp(i);
            end
        end
    end
end

results.caseType = caseType;
end

%% Efficiency Calculation Function
function efficiency = calculateH2Efficiency(H2_well, nbuildUp, nrest, ninject, nidle, nprod, ncycles,nidle1)
% Calculate H2 production efficiency
%
% Parameters:
%   H2_well - H2 well data over time
%   nbuildUp, nrest, ninject, nidle, nprod - Time step counts
%   ncycles - Number of cycles

ndeb0 = nbuildUp + nrest;
njcycle = ninject + nidle + nprod+nidle1;

mH2_injected = sum(H2_well(1:nbuildUp));
mH2_produced = 0.0;

for cycle = 1:ncycles
    ndebi = ndeb0 + (cycle-1)*njcycle;
    nj1 = ndebi + ninject;

    mH2_injected = mH2_injected + sum(H2_well(ndebi+1:nj1));
    mH2_produced = mH2_produced + sum(H2_well(nj1+nidle+1:nj1+nidle+nprod));
end

efficiency = abs(mH2_produced / mH2_injected) * 100;
end

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>