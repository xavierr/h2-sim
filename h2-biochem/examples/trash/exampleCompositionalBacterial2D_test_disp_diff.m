%% 2D Compositional Hydrogen Storage with Bio‑Clogging: Diffusion & Dispersion Analysis
% ===========================================================================
% This script compares four scenarios of a bio‑clogged hydrogen storage
% simulation with bacteria effects enabled:
%   1) No molecular diffusion, no mechanical dispersion (pure advection)
%   2) Molecular diffusion only
%   3) Mechanical dispersion only
%   4) Both molecular diffusion and mechanical dispersion
%
% The aim is to isolate the impact of these two transport mechanisms on
% the evolution of hydrogen, CO2, CH4 distributions, and on biomass‑induced
% porosity/permeability changes.
% ===========================================================================

mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr test-suite spe10
mrstModule add compositional h2-biochem h2store

%% 1. Load and setup the base model (common to all scenarios)
baseName = 'H2_STORAGE_DOME_TRAP';
dataPath = getDatasetPath('h2storage');
dataFile = fullfile(dataPath, 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

% Black‑oil model and schedule (10 injection cycles)
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

%% 2. Bio‑clogging parameters (identical for all four cases)
bacteriamodel = true;       % bacteria always present
clogModel     = false;       % clogging always active
nbact0 = 5;                 % initial bacteria concentration
nc_bact = 120;              % critical concentration
cp = 1.0;                   % clogging coefficient

% Setup base bio‑clogging model and store original porosity/permeability
[model0, poro0, perm0] = setupBioCloggingModel(model0, nbact0, nc_bact, cp, clogModel);

%% 3. Update schedule controls for compositional injection
schedule = scheduleBo;
for i = 1:numel(schedule.control)
    schedule.control(i).W.compi = [0, 1];  % well components (water, oil)
    if strcmp(schedule.control(i).W.name, 'cushion') && i < 11
        schedule.control(i).W.components = [0.0, 0.1, 0.9, 0.0];
    else
        schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0];
    end
    schedule.control(i).W.T = T0;          % injection temperature
    schedule.control(i).bc = [];           % remove boundary conditions
end

%% 4. Common initial state (with bacteria)
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);

%% 5. Helper to create BiochemistryModel with given diffusion/dispersion flags
createModel = @(doDiff, doDisp) BiochemistryModel( ...
    model0.G, model0.rock, model0.fluid, compFluid, biochemFluid, ...
    false, diagonal_backend, ...               % false = don't re‑initialise from deck
    'oil', true, 'gas', true, ...
    'bacteriamodel', true, ...
    'bactDiffusion', true, ...
    'molecularDiffusion', doDiff, ...
    'molecularDispersion', doDisp, ...
    'liquidPhase', 'O', 'vaporPhase', 'G');

%% 6. Setup optimized solver for compositional biochemical models with diffusion
model_NoNo = createModel(false, false);
[nls, ~] = setupOptimizedLinearSolver(model_NoNo, ...
    'complexityLevel', 'medium', ...
    'solverTolerance', 1e-4, ...
    'maxNonlinIter', 15, ...
    'cprDamp', 0.7);
state0 = initCompositionalStateBacteria(model_NoNo, state0.pressure, T0, state0.s, comp0, nbact0, EOS);

% Initialize timing storage
cpu_times = struct();

fprintf('\n========================================================\n');
fprintf('Starting simulations with RestartStep = 1\n');
fprintf('========================================================\n\n');

% 6.1 No diffusion, no dispersion (simplest case: lowest AMG complexity)
fprintf('>>> Scenario 1: No Diffusion, No Dispersion\n');
[nls_00, ~] = setupOptimizedLinearSolver(model_NoNo, ...
    'complexityLevel', 'low', ...
    'solverTolerance', 1e-4, ...
    'maxNonlinIter', 15);
prob_00 = packSimulationProblem(state0, model_NoNo, schedule, ...
     [baseName '_noDiff_noDisp'], 'NonLinearSolver', nls_00);
tic;
simulatePackedProblem(prob_00,'RestartStep',1);
cpu_times.noDiff_noDisp = toc;
[ws00, states00] = getPackedSimulatorOutput(prob_00);
fprintf('   CPU Time: %.2f seconds\n\n', cpu_times.noDiff_noDisp);

% 6.2 Diffusion only (medium complexity)
fprintf('>>> Scenario 2: Diffusion Only\n');
model_D = createModel(true, false);
state0_D = initCompositionalStateBacteria(model_D, state0.pressure, T0, state0.s, comp0, nbact0, EOS);
[nls_D, ~] = setupOptimizedLinearSolver(model_D, ...
    'complexityLevel', 'medium', ...
    'solverTolerance', 1e-4, ...
    'maxNonlinIter', 15);
prob_D = packSimulationProblem(state0_D, model_D, schedule, ...
    [baseName '_diffOnly'], 'NonLinearSolver', nls_D);
tic;
simulatePackedProblem(prob_D,'RestartStep',1);
cpu_times.diffOnly = toc;
[wsD, statesD] = getPackedSimulatorOutput(prob_D);
fprintf('   CPU Time: %.2f seconds\n\n', cpu_times.diffOnly);

% 6.3 Dispersion only (medium complexity)
fprintf('>>> Scenario 3: Dispersion Only\n');
model_Dp = createModel(false, true);
state0_Dp = initCompositionalStateBacteria(model_Dp, state0.pressure, T0, state0.s, comp0, nbact0, EOS);
[nls_Dp, ~] = setupOptimizedLinearSolver(model_Dp, ...
    'complexityLevel', 'medium', ...
    'solverTolerance', 1e-4, ...
    'maxNonlinIter', 15);
prob_Dp = packSimulationProblem(state0_Dp, model_Dp, schedule, ...
    [baseName '_dispOnly'], 'NonLinearSolver', nls_Dp);
tic;
simulatePackedProblem(prob_Dp, 'RestartStep',1);
cpu_times.dispOnly = toc;
%[wsDp, statesDp] = getPackedSimulatorOutput(prob_Dp);
fprintf('   CPU Time: %.2f seconds\n\n', cpu_times.dispOnly);

% 6.4 Both diffusion and dispersion (highest complexity: conservative AMG settings)
fprintf('>>> Scenario 4: Both Diffusion & Dispersion\n');
model_DD = createModel(true, true);
state0_DD = initCompositionalStateBacteria(model_DD, state0.pressure, T0, state0.s, comp0, nbact0, EOS);
[nls_DD, ~] = setupOptimizedLinearSolver(model_DD, ...
    'complexityLevel', 'high', ...
    'solverTolerance', 1e-4, ...
    'maxNonlinIter', 15);
prob_DD = packSimulationProblem(state0_DD, model_DD, schedule, ...
    [baseName '_bothDiffDisp'], 'NonLinearSolver', nls_DD);
tic;
simulatePackedProblem(prob_DD,'RestartStep',1);
cpu_times.bothDiffDisp = toc;
[wsDD, statesDD] = getPackedSimulatorOutput(prob_DD);
fprintf('   CPU Time: %.2f seconds\n\n', cpu_times.bothDiffDisp);
% 
fprintf('========================================================\n');
fprintf('All simulations completed!\n');
fprintf('========================================================\n\n');

%% 7. Display CPU Times
fprintf('\n==================== CPU TIMES ====================\n');
fprintf('Scenario                            Time (seconds)\n');
fprintf('--------------------------------------------------\n');
fprintf('1. No Diffusion, No Dispersion:     %8.2f s\n', cpu_times.noDiff_noDisp);
fprintf('2. Diffusion Only:                  %8.2f s\n', cpu_times.diffOnly);
fprintf('3. Dispersion Only:                 %8.2f s\n', cpu_times.dispOnly);
fprintf('4. Both Diffusion & Dispersion:     %8.2f s\n', cpu_times.bothDiffDisp);
fprintf('==================================================\n');

% Calculate speedup
base_time = cpu_times.noDiff_noDisp;
fprintf('\nSpeedup relative to No Diffusion/No Dispersion:\n');
fprintf('  Diffusion Only:                  %.2fx\n', base_time/cpu_times.diffOnly);
fprintf('  Dispersion Only:                 %.2fx\n', base_time/cpu_times.dispOnly);
fprintf('  Both Diffusion & Dispersion:     %.2fx\n', base_time/cpu_times.bothDiffDisp);

% Save CPU times
save('cpu_times.mat', 'cpu_times');
fprintf('\nCPU times saved to cpu_times.mat\n');

%% 8. Visualise results – states (toolbar plots)
figure; plotToolbar(model0.G, states00); title('No Diffusion, No Dispersion');
figure; plotToolbar(model0.G, statesD);  title('Diffusion Only');
figure; plotToolbar(model0.G, statesDp); title('Dispersion Only');
figure; plotToolbar(model0.G, statesDD); title('Both Diffusion & Dispersion');

%% 9. Visualise well outputs – comparison
% Collect well solutions for all scenarios
ws_combined = {ws00, wsD, wsDp, wsDD};
figure;
plotWellSols(ws_combined, cumsum(schedule.step.val)/day, ...
    'datasetnames', {'No Diff/No Disp', 'Diffusion Only', 'Dispersion Only', 'Both Diff & Disp'});
title('Well Solutions Comparison');

%% 10. Porosity comparison at final state
figure;
subplot(2,2,1); plotCellData(model0.G, states00{end}.rock.poro);
title('No Diff/No Disp'); colorbar;
subplot(2,2,2); plotCellData(model0.G, statesD{end}.rock.poro);
title('Diffusion Only'); colorbar;
subplot(2,2,3); plotCellData(model0.G, statesDp{end}.rock.poro);
title('Dispersion Only'); colorbar;
subplot(2,2,4); plotCellData(model0.G, statesDD{end}.rock.poro);
title('Both Diff & Disp'); colorbar;
sgtitle('Final Porosity Comparison');

%% 11. H2 Loss Analysis (if function exists)
try
    figure;
    plotH2Loss(model0, schedule, states00, ws00);
    title('H2 Loss - No Diff/No Disp');

    figure;
    plotH2Loss(model_D, schedule, statesD, wsD);
    title('H2 Loss - Diffusion Only');

    figure;
    plotH2Loss(model_Dp, schedule, statesDp, wsDp);
    title('H2 Loss - Dispersion Only');

    figure;
    plotH2Loss(model_DD, schedule, statesDD, wsDD);
    title('H2 Loss - Both Diff & Disp');
catch
    fprintf('Note: plotH2Loss function not available. Skipping.\n');
end

%% 12. Plot CPU Time Comparison
figure;
scenario_names = {'No Diff/No Disp', 'Diffusion Only', 'Dispersion Only', 'Both Diff & Disp'};
cpu_values = [cpu_times.noDiff_noDisp, cpu_times.diffOnly, cpu_times.dispOnly, cpu_times.bothDiffDisp];
bar(cpu_values);
set(gca, 'XTickLabel', scenario_names);
ylabel('CPU Time (seconds)');
title('CPU Time Comparison for Different Transport Mechanisms');
grid on;

% Add values on top of bars
for i = 1:length(cpu_values)
    text(i, cpu_values(i) + 0.5, sprintf('%.1f s', cpu_values(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

%% 13. Summary Statistics
fprintf('\n==================== SUMMARY STATISTICS ====================\n');
fprintf('Total simulation time (all scenarios): %.2f seconds (%.2f minutes)\n', ...
    sum(cpu_values), sum(cpu_values)/60);
fprintf('Average time per scenario: %.2f seconds\n', mean(cpu_values));
fprintf('============================================================\n');
% 
% %% Copyright notice
% % <html>
% % <p><font size="-1">
% % Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.
% % </font></p>
% % </html>