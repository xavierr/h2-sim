%% Simple comparison of bacterial diffusion in a bio-clogged H2 storage
% ===========================================================================
% Two cases:
%   1) No bacterial diffusion (bactDiffusion = false)
%   2) With bacterial diffusion (bactDiffusion = true)
%
% Both include bio-clogging and the methanogenic reaction.
% The script shows the effect on biomass distribution and porosity.
% ===========================================================================

mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr test-suite spe10
mrstModule add compositional h2-biochem h2store

%% 1. Load and setup the base model (common to both scenarios)
baseName = 'H2_STORAGE_DOME_TRAP';
dataPath = getDatasetPath('h2storage');
dataFile = fullfile(dataPath, 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

% Black‑oil model and schedule (2 injection cycles to speed up)
[~, ~, state0Bo, modelBo, scheduleBo, ~] = modelForSimple2DAquifer(deck, 'numcycles', 10);

% Convert to compositional model
model0 = convertBlackOilModelToCompositionalModel(modelBo);
state0 = convertBlackOilStateToCompositional(modelBo, state0Bo);

% Fluid and EOS
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
    {'H2O', 'H2', 'CO2', 'C1'});
biochemFluid = TableBioChemMixture({'MethanogenicArchae'});
EOS = SoreideWhitsonEos([], compFluid);
model0.EOSModel = EOS;

nc = model0.G.cells.num;
T0 = 273.15 + 44.35;  % Initial temperature [K]
comp0 = repmat([0.8480, 1.0e-5, 1.0e-5, 0.1530], nc, 1);

%% 2. Bio‑clogging parameters (identical for both cases)
bacteriamodel = true;
clogModel     = false;
nbact0 = 5;                 % initial bacteria concentration
nc_bact = 120;              % critical concentration
cp = 1.0;                   % clogging coefficient

% Setup base bio‑clogging model and store original porosity/permeability
[model0, poro0, perm0] = setupBioCloggingModel(model0, nbact0, nc_bact, cp, clogModel);

%% 3. Update schedule controls for compositional injection
schedule = scheduleBo;
bc  = addBC([], schedule.control(1).bc.face, 'flux', 0, 'sat',[0.1514 0.8486] );    % no-flow for the flow eqs
nc = 4; % number of components
nf = numel(bc.face); % number of boundary faces
bc.components = repmat([0.0; 0.95; 0.05; 0.0], 1, nf)';
bc  = addBacterialBC(bc, bc.face, nbact0);          % Dirichlet nbact at those faces
for i = 1:numel(schedule.control)
    schedule.control(i).W.compi = [0, 1];  % well components (water, oil)
    if strcmp(schedule.control(i).W.name, 'cushion') && i < 3
        schedule.control(i).W.components = [0.0, 0.1, 0.9, 0.0];
    else
        schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0];
    end
    schedule.control(i).W.T = T0;
    schedule.control(i).bc = [];

end
%% 4. Common initial state (with bacteria)
backend = DiagonalAutoDiffBackend('modifyOperators', true);

%% 5. Helper to create a BiochemistryModel with optional bacterial diffusion

createModel = @(doBactDiff) createBactDiffModel(model0, compFluid, biochemFluid, backend, doBactDiff);

%% 6. Create and solve the two scenarios
% 6.1 No bacterial diffusion
model_noBactDiff = createModel(false);
state0 = initCompositionalStateBacteria(model_noBactDiff, state0.pressure, T0, ...
    state0.s, comp0, nbact0, EOS);
nls = NonLinearSolver();
nls.LinearSolver = selectLinearSolverAD(model_noBactDiff);
linsolve = nls.LinearSolver;
disp(linsolve)
nls = NonLinearSolver('LinearSolver',linsolve);
prob_noBactDiff = packSimulationProblem(state0, model_noBactDiff, schedule, ...
    [baseName '_bio_noDiff_noDisp'], 'NonLinearSolver', nls);
[ws_noDiff, states_noDiff] = getPackedSimulatorOutput(prob_noBactDiff);
mrstVerbose true;
fprintf('\n===== Case 1 Complete (No Bacterial Diffusion) =====\n\n');
% 6.2 With bacterial diffusion
model_withBactDiff = createModel(true);
state0_with = initCompositionalStateBacteria(model_withBactDiff, state0.pressure, T0, ...
    state0.s, comp0, nbact0, EOS);
nls.LinearSolver = selectLinearSolverAD(model_withBactDiff);
linsolve = nls.LinearSolver;
disp(linsolve)
nls = NonLinearSolver('LinearSolver',linsolve);
prob_withBactDiff = packSimulationProblem(state0_with, model_withBactDiff, schedule, ...
    [baseName '_withBactDiff'], 'NonLinearSolver', nls);
simulatePackedProblem(prob_withBactDiff,'RestartStep',1);
[ws_withDiff, states_withDiff] = getPackedSimulatorOutput(prob_withBactDiff);
fprintf('\n===== Case 2 Complete (With Bacterial Diffusion + AMGCL Tuning) =====\n\n');

%% 7. Visualise biomass at final time
figure;
subplot(1,2,1);
plotCellData(model0.G, states_noDiff{end}.nbact); title('No Bacterial Diffusion');
colorbar; view(0,90);
subplot(1,2,2);
plotCellData(model0.G, states_withDiff{end}.nbact); title('With Bacterial Diffusion');
colorbar; view(0,90);
sgtitle('Biomass concentration (cells/m³)');

%% 8. Porosity at final time
figure;
subplot(1,2,1);
plotCellData(model0.G, states_noDiff{end}.rock.poro); title('No Bacterial Diffusion');
colorbar; view(0,90);
subplot(1,2,2);
plotCellData(model0.G, states_withDiff{end}.rock.poro); title('With Bacterial Diffusion');
colorbar; view(0,90);
sgtitle('Final porosity');

%% 9. Well outputs comparison
ws_combined = {ws_noDiff, ws_withDiff};
figure; plotWellSols(ws_combined, cumsum(schedule.step.val)/day, ...
    'datasetnames', {'No Bact Diffusion', 'With Bact Diffusion'});

% ===========================================================================
% Helper function to add bacterial diffusion to a BiochemistryModel
% ===========================================================================
function model = createBactDiffModel(model0, compFluid, biochemFluid, backend, bactDiffusion)
    % Create the base biochemistry model (without internal flags for diffusion)
    model = BiochemistryModel(...
        model0.G, model0.rock, model0.fluid, compFluid, biochemFluid, ...
        false, backend, ...
        'oil', true, 'gas', true, ...
        'bacteriamodel', true, ...
        'molecularDiffusion', false, ...
        'bactDiffusion', bactDiffusion,...
        'molecularDispersion', false, ...
        'liquidPhase', 'O', 'vaporPhase', 'G');

    % Ensure the bacteria mass balance equation exists (it should be in model.Equations)
    % Add bacterial diffusion flux and extra term to the equation if requested
   
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}