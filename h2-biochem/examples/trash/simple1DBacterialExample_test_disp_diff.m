%% 1D Compositional Simulation with Bio-Clogging: Analysis of Diffusion & Dispersion
% ===========================================================================
% This script compares four scenarios of a bio-clogged compositional
% simulation with bacteria effects enabled:
%   1) No molecular diffusion, no mechanical dispersion (pure advection)
%   2) Molecular diffusion only
%   3) Mechanical dispersion only
%   4) Both molecular diffusion and mechanical dispersion
%
% The goal is to isolate the impact of these two physical mechanisms on
% component transport (e.g. hydrogen, methane, CO2) and biomass evolution.
% The comparison uses a 1D two-phase (liquid/gas) composition model based
% on the Soreide-Whitson EOS and includes bio-clogging feedback.
%
% Reference: https://www.sciencedirect.com/science/article/pii/S0360319925039473
% ===========================================================================

mrstModule add ad-core ad-props deckformat mrst-gui
mrstModule add compositional h2-biochem
mrstVerbose off;

%% 1. Setup base compositional model (common geometry, fluids, wells)
[~, model0, schedule, ~] = setupSimpleCompositionalExample(false);

% Increase porosity in boundary cells to mimic well effects
model0.rock.poro([1, end]) = 0.90;
model0.rock.perm = model0.rock.perm .* (model0.rock.poro ./ model0.rock.poro).^3; % simple perm update

model0.fluid.rhoGS = 0.08988;   % gas density at STP [kg/m3]
model0.fluid.rhoOS = 999.0140;  % liquid density at STP [kg/m3]

% Shorten schedule for faster runs
schedule.step.val     = schedule.step.val / 3;
schedule.step.val     = schedule.step.val(1:400);
schedule.step.control = schedule.step.control(1:400);

%% 2. Common physics: EOS, relperm, initial state, bio-clogging
% EOS
model0.EOSModel = SoreideWhitsonEos( ...
    model0.G, model0.EOSModel.CompositionalMixture, 'msalt', 5);

% Residual saturations and relperm
swc   = 0.20; sor = 0.20; sgr = 0.10;
swmax = 0.90; sgmax = 0.80;
model0.fluid = modifyRelPermForResidualSaturations( ...
    model0.fluid, swc, swmax, sgr, sor, sgmax);

% Initial overall composition [H2O, H2, C1, CO2]
initComp  = [0.90, 0.045, 0.005, 0.05];
initTemp  = 273.15 + 40;   % K
initPress = 82 * barsa;

% Bio-clogging model (same for all cases; we keep bacteria true)
initBact = 9;   % normalized bacteria concentration
nc       = 180; % critical concentration
cp       = 0.0; % scaling coefficient
clogModel = false;
[model0, poro0, perm0] = setupBioCloggingModel(model0, initBact, nc, cp, clogModel);

% Shut wells (no injection/production, flow driven by initial gradients)
schedule.control.W(1).components = [0.001, 0.958, 0.001, 0.05];
schedule.control.W(2).components = [0.001, 0.998, 0.001, 0.05];
schedule.control.W(1).type = 'rate'; schedule.control.W(1).val = 0;
schedule.control.W(2).type = 'rate'; schedule.control.W(2).val = 0;

%% 3. Define the four scenarios via model flags (bacteria always true)
% We create a fresh BiochemistryModel for each scenario to ensure correct
% state function registration.

compFluid = TableCompositionalMixture( ...
    {'Water','Hydrogen','Methane','CarbonDioxide'}, ...
    {'H2O','H2','C1','CO2'});
backend = DiagonalAutoDiffBackend('modifyOperators', true);
biochemFluid = TableBioChemMixture({'MethanogenicArchae'});

scenarios = {};

% Helper function to create model with given flags
createModel = @(doDiff, doDisp) BiochemistryModel( ...
    model0.G, model0.rock, model0.fluid, compFluid, biochemFluid, ...
    true, backend, 'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', true, ...
    'molecularDiffusion', doDiff, ...
    'molecularDispersion', doDisp, ...
    'liquidPhase','O', 'vaporPhase','G');

% 1) No diffusion, no dispersion
model_00 = createModel(false, false);
% Common initial state (bacteria concentration same)
state0 = initCompositionalStateBacteria(model_00, initPress, initTemp, [0,1], ...
    initComp, initBact, model_00.EOSModel);
prob00 = packSimulationProblem(state0, model_00, schedule, 'bio_noDiff_noDisp');
prob00.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(prob00, 'RestartStep',1);
[~, st00] = getPackedSimulatorOutput(prob00);
scenarios{end+1} = struct('name','No Diff, No Disp', 'states',{st00}, 'color',[0.7 0 0], 'line','-');

% 2) Diffusion only
model_Diff = createModel(true, false);
probDiff = packSimulationProblem(state0, model_Diff, schedule, 'bio_diffOnly');
probDiff.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(probDiff,'RestartStep',1);
[~, stDiff] = getPackedSimulatorOutput(probDiff);
scenarios{end+1} = struct('name','Diffusion Only', 'states',{stDiff}, 'color',[0 0.7 0.7], 'line','--');

% 3) Dispersion only
model_Disp = createModel(false, true);
probDisp = packSimulationProblem(state0, model_Disp, schedule, 'bio_dispOnly');
probDisp.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(probDisp,'RestartStep',1);
[~, stDisp] = getPackedSimulatorOutput(probDisp);
scenarios{end+1} = struct('name','Dispersion Only', 'states',{stDisp}, 'color',[0 0.5 0], 'line','-.');

% 4) Both
model_Both = createModel(true, true);
probBoth = packSimulationProblem(state0, model_Both, schedule, 'bio_diffAndDisp');
probBoth.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(probBoth,'RestartStep',1);
[~, stBoth] = getPackedSimulatorOutput(probBoth);
scenarios{end+1} = struct('name','Both Diff & Disp', 'states',{stBoth}, 'color',[0 0 0], 'line','-');

%% 4. Compare scenarios
timeDays = cumsum(schedule.step.val) / day;
timeYrs  = timeDays / 365;
compNames = model_00.EOSModel.getComponentNames();
idxH2  = find(strcmp(compNames,'H2'));
idxCO2 = find(strcmp(compNames,'CO2'));
idxCH4 = find(strcmp(compNames,'C1'));
% Simple plot: final time step overall mole fraction profiles
figure('Name','Final H2 profile');
for i = 1:numel(scenarios)
    x = scenarios{i}.states{100}.x(:,idxH2);
    plot(model0.G.cells.centroids(:,1), x, 'Color', scenarios{i}.color, ...
        'LineStyle', scenarios{i}.line, 'LineWidth', 1.5);
    hold on;
end
hold off; legend({scenarios.name}, 'Location','best');
xlabel('Distance [m]'); ylabel('H2 mole fraction');

figure('Name','Final CO2 profile');
for i = 1:numel(scenarios)
    x = scenarios{i}.states{end}.x(:,idxCO2);
    plot(model0.G.cells.centroids(:,1), x, 'Color', scenarios{i}.color, ...
        'LineStyle', scenarios{i}.line, 'LineWidth', 1.5);
    hold on;
end
hold off; legend({scenarios.name}, 'Location','best');
xlabel('Distance [m]'); ylabel('CO2 mole fraction');

% Porosity profile (bio-clogging effect)
figure('Name','Porosity profile');
for i = 1:numel(scenarios)
    phi = model_00.rock.poro(scenarios{i}.states{end}.pressure,scenarios{i}.states{end}.nbact);
    plot(model0.G.cells.centroids(:,1), phi, 'Color', scenarios{i}.color, ...
        'LineStyle', scenarios{i}.line, 'LineWidth', 1.5);
    hold on;
end
hold off; legend({scenarios.name}, 'Location','best');
xlabel('Distance [m]'); ylabel('Porosity');

%% 5. Copyright notice
% <html>
% <p><font size="-1">
% Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% ...
% </html>