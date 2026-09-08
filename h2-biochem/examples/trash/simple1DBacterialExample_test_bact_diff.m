%% 1D Bio‑clogging: Effect of bacterial diffusion on biomass and composition
% ===========================================================================
% Two cases:
%   1) No bacterial diffusion (bactDiffusion = false)
%   2) With bacterial diffusion (bactDiffusion = true)
%
% Both include bio‑clogging and the methanogenic reaction.
% The script shows the effect on biomass distribution, porosity, and
% component profiles.
% ===========================================================================
mrstVerbose true;
mrstModule add ad-core ad-props deckformat mrst-gui
mrstModule add compositional h2-biochem
mrstVerbose off;

%% 1. Setup base compositional model (common geometry, fluids, wells)
[~, model0, schedule, ~] = setupSimpleCompositionalExample(false);

% Increase porosity in boundary cells to mimic well effects
model0.rock.poro([1, end]) = 0.90;
model0.rock.perm = model0.rock.perm .* (model0.rock.poro ./ model0.rock.poro).^3;

model0.fluid.rhoGS = 0.08988;   % kg/m³
model0.fluid.rhoOS = 999.0140;

% Shorten schedule
schedule.step.val     = schedule.step.val / 3;
schedule.step.val     = schedule.step.val(1:400);
schedule.step.control = schedule.step.control(1:400);

%% 2. Common physics: EOS, relperm, bio‑clogging, shut wells
model0.EOSModel = SoreideWhitsonEos( ...
    model0.G, model0.EOSModel.CompositionalMixture, 'msalt', 5);

% Residual saturations
swc = 0.20; sor = 0.20; sgr = 0.10;
swmax = 0.90; sgmax = 0.80;
model0.fluid = modifyRelPermForResidualSaturations( ...
    model0.fluid, swc, swmax, sgr, sor, sgmax);

% Initial composition, temperature, pressure
initComp  = [0.90, 0.045, 0.005, 0.05];
initTemp  = 273.15 + 40;   % K
initPress = 82 * barsa;

% Bio‑clogging (same for both cases)
initBact = 9;   % normalized bacteria concentration
clogModel = false;
nc       = 180; % critical concentration
cp       = 0.0; % scaling coefficient
[model0, poro0, perm0] = setupBioCloggingModel(model0, initBact, nc, cp, clogModel);

% Shut wells – transport driven by initial gradients only
schedule.control.W(1).components = [0.001, 0.958, 0.001, 0.05];
schedule.control.W(2).components = [0.001, 0.998, 0.001, 0.05];
schedule.control.W(1).type = 'rate'; schedule.control.W(1).val = 0;
schedule.control.W(2).type = 'rate'; schedule.control.W(2).val = 0;

%% 3. Create models with/without bacterial diffusion
compFluid = TableCompositionalMixture( ...
    {'Water','Hydrogen','Methane','CarbonDioxide'}, ...
    {'H2O','H2','C1','CO2'});
biochemFluid = TableBioChemMixture({'MethanogenicArchae'});
% Ensure diffusion coefficient is set (required by DiffusiveBactFlux)
if ~isfield(biochemFluid, 'bactdiff') || isempty(biochemFluid.bactdiff)
    biochemFluid.bactdiff = 1e-8;   % m²/s
end

backend = DiagonalAutoDiffBackend('modifyOperators', true);

% Helper to build a BiochemistryModel with optional bacterial diffusion
createModel = @(bactDiff) createBactDiffModel(model0, compFluid, biochemFluid, backend, bactDiff);

model_noDiff = createModel(false);
model_withDiff = createModel(true);

%% 4. Initialise state (identical for both)
state0 = initCompositionalStateBacteria(model_noDiff, initPress, initTemp, [0,1], ...
    initComp, initBact, model_noDiff.EOSModel);

%% 5. Solve
nls = NonLinearSolver();
nls.LinearSolver = selectLinearSolverAD(model_noDiff);
nls.LinearSolver.verbose = true;
% No bacterial diffusion
prob_noDiff = packSimulationProblem(state0, model_noDiff, schedule, 'bio_bactDiff_off',  'NonLinearSolver', nls);
prob_noDiff.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(prob_noDiff);
[st_noDiff, ws_noDiff] = getPackedSimulatorOutput(prob_noDiff);

% With bacterial diffusion
prob_withDiff = packSimulationProblem(state0, model_withDiff, schedule, 'bio_bactDiff_on',  'NonLinearSolver', nls);
prob_withDiff.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(prob_withDiff,'RestartStep',1);
[st_withDiff, ws_withDiff] = getPackedSimulatorOutput(prob_withDiff);

%% 6. Compare results
timeDays = cumsum(schedule.step.val) / day;
cnames   = model_noDiff.EOSModel.getComponentNames();
idxH2    = find(strcmp(cnames, 'H2'));
idxCO2   = find(strcmp(cnames, 'CO2'));
idxCH4   = find(strcmp(cnames, 'C1'));
G = model0.G;
xc = G.cells.centroids(:,1);

% Final H2 profile
figure('Name','Final H2 profile');
plot(xc, ws_noDiff{end}.x(:, idxH2), 'b-', 'LineWidth', 1.5); hold on;
plot(xc, ws_withDiff{end}.x(:, idxH2), 'r--', 'LineWidth', 1.5);
hold off; legend('No bact diffusion','With bact diffusion','Location','best');
xlabel('Distance [m]'); ylabel('H_2 mole fraction');

% Final CO2 profile
figure('Name','Final CO2 profile');
plot(xc, ws_noDiff{end}.x(:, idxCO2), 'b-', 'LineWidth', 1.5); hold on;
plot(xc, ws_withDiff{end}.x(:, idxCO2), 'r--', 'LineWidth', 1.5);
hold off; legend('No bact diffusion','With bact diffusion','Location','best');
xlabel('Distance [m]'); ylabel('CO_2 mole fraction');

% Biomass (nbact) profile
figure('Name','Biomass (nbact)');
plot(xc, ws_noDiff{end}.nbact, 'b-', 'LineWidth', 1.5); hold on;
plot(xc, ws_withDiff{end}.nbact, 'r--', 'LineWidth', 1.5);
hold off; legend('No bact diffusion','With bact diffusion','Location','best');
xlabel('Distance [m]'); ylabel('Normalized bacteria concentration');

% Porosity profile (computed from dynamic model)
if clogModel
    phi_noDiff = model0.rock.poro(ws_noDiff{end}.pressure, ws_noDiff{end}.nbact);
    phi_withDiff = model0.rock.poro(ws_withDiff{end}.pressure, ws_withDiff{end}.nbact);
    figure('Name','Porosity');
    plot(xc, phi_noDiff, 'b-', 'LineWidth', 1.5); hold on;
    plot(xc, phi_withDiff, 'r--', 'LineWidth', 1.5);
    hold off; legend('No bact diffusion','With bact diffusion','Location','best');
    xlabel('Distance [m]'); ylabel('Porosity');
end
%% 7. Helper function: build model with optional bacterial diffusion
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