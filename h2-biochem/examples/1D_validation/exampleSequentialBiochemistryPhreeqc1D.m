function summary = exampleSequentialBiochemistryPhreeqc1D(varargin)
% Coarse-flow/local-reaction split applied to the 1D H2-storage benchmark.
%
% This experimental example demonstrates SequentialBiochemistryPhreeqcModel,
% a general coarse-flow/local-reaction sequential model: ONE global
% compositional flow-and-transport solve per schedule/control step (with
% biological reaction sources disabled), followed by local MRST
% biochemical-reaction substeps (no spatial flux/well/BC forcing), each
% followed by a PHREEQC equilibrium/chemistry-feedback update. Unlike
% simulateSequentialH2BiochemPhreeqc, there is no outer Picard iteration
% between the two stages.
%
% The setup matches the high-rate sequential h2-biochem-PHREEQC case used
% by runThreeBackendComparison (paperBiomassKinetics, N/N0 = 1, pH 6.24,
% HCO3- = 1.370e-3, zero initial overall CO2, high-rate kinetics, complete
% injection/storage/production schedule), run with a coarse schedule and
% two-day flow steps and 0.4 day maximum local-reaction substeps. It is run with the plain
% simulateScheduleAD driver -- no bespoke outer coupling loop is required.
%
% Note: runThreeBackendComparison is intentionally left unchanged; this
% script is a separate, standalone entry point for the new model.
%
% SYNOPSIS:
%   summary = exampleSequentialBiochemistryPhreeqc1D();
%   summary = exampleSequentialBiochemistryPhreeqc1D( ...
%       'phreeqcDatabaseFile', db, 'reactionSubstepMaxDt', 0.4*day);
%
% OPTIONAL PARAMETERS:
%   phreeqcDatabaseFile  - Absolute path to PHREEQC_Modified.DAT. Falls
%                          back to the PHREEQC_DATABASE_FILE environment
%                          variable, then to `which('PHREEQC_Modified.DAT')`.
%   reactionSubstepMaxDt - Maximum duration of one local reaction/PHREEQC
%                          substep (default 0.4 day, matching the
%                          hybrid Picard reference used elsewhere in this
%                          module).
%   flowTimestepMaxDt     - Maximum coarse compositional-flow timestep
%                          (default 2 days, giving 125 flow steps over the
%                          complete 250-day benchmark).
%   maximumFlowSteps      - Optional positive integer limit for short
%                          validation runs (default Inf: complete schedule).
%   reactionSolverVerbose - Print local reaction Newton residuals
%                          (default false).
%
% RETURNS:
%   summary - Struct with the model, schedule, states, well solutions,
%             report, and H2 consumption/runtime metrics.
%
% SEE ALSO:
%   SequentialBiochemistryPhreeqcModel,
%   convertToSequentialBiochemistryPhreeqcModel,
%   setupH2StorageExampleWithSRB_benchmark, runThreeBackendComparison

    mrstModule add ad-core ad-props compositional deckformat h2-biochem

    opt = struct( ...
        'phreeqcDatabaseFile', '', ...
        'flowTimestepMaxDt', 2*day, ...
        'reactionSubstepMaxDt', 0.4*day, ...
        'maximumFlowSteps', inf, ...
        'reactionSolverVerbose', false);
    opt = merge_options(opt, varargin{:});
    validateattributes(opt.reactionSubstepMaxDt, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'reactionSubstepMaxDt');
    validateattributes(opt.flowTimestepMaxDt, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'flowTimestepMaxDt');
    assert(isnumeric(opt.maximumFlowSteps) && isscalar(opt.maximumFlowSteps) && ...
        isreal(opt.maximumFlowSteps) && opt.maximumFlowSteps > 0 && ...
        (isinf(opt.maximumFlowSteps) || ...
         (isfinite(opt.maximumFlowSteps) && ...
          opt.maximumFlowSteps == floor(opt.maximumFlowSteps))), ...
        'maximumFlowSteps must be a positive integer or Inf.');
    validateattributes(opt.reactionSolverVerbose, {'logical', 'numeric'}, ...
        {'scalar', 'real', 'finite'}, mfilename, 'reactionSolverVerbose');
    databaseFile = resolvePhreeqcDatabaseFile(opt.phreeqcDatabaseFile);

    % Matched high-rate sequential h2-biochem-PHREEQC setup (see
    % runThreeBackendComparison's "Sequential h2-biochem PHREEQC" case).
    [~, baseModel, schedule, state0] = setupH2StorageExampleWithSRB_benchmark( ...
        'rate', 'highrate', ...
        'scheduleMode', 'complete', ...
        'injectionCO2', 0, ...
        'bacteriamodel', true, ...
        'bactDiffusion', false, ...
        'chemotaxisEffect', false, ...
        'molecularDiffusion', false, ...
        'molecularDispersion', false, ...
        'bioClogging', false, ...
        'carbonateBuffer', true, ...
        'carbonateBufferPH', 6.24, ...
        'initialHCO3', 1.370e-3, ...
        'initialOverallCO2', 0, ...
        'equilibrateInitialCO2', false, ...
        'paperBiomassKinetics', true, ...
        'nbact0', 1, ...
        'phreeqcTimestepCoupling', true, ...
        'phreeqcBackend', 'sequential-h2biochem-phreeqc', ...
        'phreeqcDatabaseFile', databaseFile);
    assert(baseModel.bacterialDecayOrder == 1, ...
        'The UGFACT-matched hybrid case requires first-order biomass decay.');
    assert(max(abs(value(state0.nbact) - 1), [], 'all') < 1e-12, ...
        'The UGFACT-matched hybrid case requires initial N/N0 = 1.');

    % Use the compositional reference's coarse flow resolution while the
    % sequential model independently subcycles local reactions.
    schedule = coarsenScheduleByControl(schedule, opt.flowTimestepMaxDt);
    if isfinite(opt.maximumFlowSteps)
        nSteps = min(opt.maximumFlowSteps, numel(schedule.step.val));
        schedule.step.val = schedule.step.val(1:nSteps);
        schedule.step.control = schedule.step.control(1:nSteps);
    end

    % Convert the fully configured benchmark model into the general
    % coarse-flow/local-reaction split model.
    model = convertToSequentialBiochemistryPhreeqcModel(baseModel, ...
        'reactionSubstepMaxDt', opt.reactionSubstepMaxDt);
    model.reactionNonLinearSolver.verbose = ...
        logical(opt.reactionSolverVerbose);

    solver = NonLinearSolver();
    solver.maxTimestepCuts = 12;
    timer = tic();
    [ws, states, report] = simulateScheduleAD(state0, model, schedule, ...
        'nonlinearSolver', solver);
    elapsedSeconds = toc(timer);
    if isstruct(report) && isfield(report, 'Failure')
        assert(~report.Failure, ...
            'SequentialBiochemistryPhreeqcModel simulation reported failure.');
    end

    nReactions = model.biochemFluid.nbioreact;
    cumulative = zeros(numel(states), nReactions);
    for reaction = 1:nReactions
        [~, reactionCumulative] = computeH2Consumption( ...
            states, schedule, model, reaction);
        cumulative(:, reaction) = sum(reactionCumulative, 1).';
    end
    consumedH2Moles = sum(cumulative(end, :), 2);
    injectedH2Moles = prescribedInjectedH2(schedule, model);

    fprintf('\nSequentialBiochemistryPhreeqcModel, 1D benchmark:\n');
    fprintf('  Control steps            : %d\n', numel(schedule.step.val));
    fprintf('  Reaction substep max dt  : %.3g day\n', ...
        opt.reactionSubstepMaxDt/day);
    fprintf('  Injected H2              : %.6g mol\n', injectedH2Moles);
    fprintf('  Consumed H2              : %.6g mol (%.3g %%)\n', ...
        consumedH2Moles, 100*consumedH2Moles/injectedH2Moles);
    fprintf('  Runtime                  : %.3g s\n', elapsedSeconds);

    summary = struct( ...
        'model', model, ...
        'schedule', schedule, ...
        'state0', state0, ...
        'states', {states}, ...
        'wellSols', {ws}, ...
        'report', report, ...
        'injectedH2Moles', injectedH2Moles, ...
        'consumedH2Moles', consumedH2Moles, ...
        'h2LossPercent', 100*consumedH2Moles/injectedH2Moles, ...
        'elapsedSeconds', elapsedSeconds);
end

function injected = prescribedInjectedH2(schedule, model)
    names = model.EOSModel.CompositionalMixture.names;
    idxH2 = findComponentIndex(names, {'H2', 'Hydrogen'});
    gasIndex = model.getVaporIndex();
    injected = 0;
    for step = 1:numel(schedule.step.val)
        wells = schedule.control(schedule.step.control(step)).W;
        for w = 1:numel(wells)
            if wells(w).sign <= 0 || ...
                    (isfield(wells, 'status') && ~wells(w).status)
                continue;
            end

            if strcmpi(wells(w).type, 'grat')
                gasRate = wells(w).val;
            elseif strcmpi(wells(w).type, 'rate')
                gasRate = wells(w).val*wells(w).compi(gasIndex);
            else
                continue;
            end
            molarRate = gasRate*model.FacilityModel.pressure/ ...
                (8.314462618*model.FacilityModel.T);
            injected = injected + molarRate*wells(w).components(idxH2)* ...
                schedule.step.val(step);
        end
    end
    assert(injected > 0, 'No prescribed H2 injection was found.');
end

function index = findComponentIndex(names, aliases)
    for i = 1:numel(aliases)
        candidate = find(strcmpi(names, aliases{i}), 1);
        if ~isempty(candidate)
            index = candidate;
            return;
        end
    end
    error('Component "%s" was not found.', strjoin(aliases, '" or "'));
end

function databaseFile = resolvePhreeqcDatabaseFile(databaseFile)
    assert(ischar(databaseFile) || ...
        (isstring(databaseFile) && isscalar(databaseFile)), ...
        'phreeqcDatabaseFile must be a character vector or scalar string.');
    databaseFile = char(databaseFile);
    if isempty(strtrim(databaseFile))
        databaseFile = getenv('PHREEQC_DATABASE_FILE');
    end
    if isempty(strtrim(databaseFile))
        databaseFile = which('PHREEQC_Modified.DAT');
    end
    if isempty(strtrim(databaseFile)) || ~isfile(databaseFile)
        error('exampleSequentialBiochemistryPhreeqc1D:MissingPhreeqcDatabase', ...
            ['PHREEQC_Modified.DAT was not found. Pass its absolute path ', ...
            'using ''phreeqcDatabaseFile'', databaseFile.']);
    end
    [~, name, extension] = fileparts(databaseFile);
    assert(strcmpi([name, extension], 'PHREEQC_Modified.DAT'), ...
        'The database must be PHREEQC_Modified.DAT.');
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
