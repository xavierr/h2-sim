classdef SequentialBiochemistryPhreeqcModel < BiochemistryPhreeqcModel
    % General coarse-flow/local-reaction sequential meta-model.
    %
    % SYNOPSIS:
    %   model = SequentialBiochemistryPhreeqcModel(G, rock, fluid, ...
    %       compFluid, biochemFluid, includeWater, backend, ...
    %       'reactionSubstepMaxDt', 0.4*day, 'pn1', vn1, ...)
    %
    % DESCRIPTION:
    %   For every control-step timestep dt, this model:
    %
    %     1. Solves ONE complete compositional flow-and-transport stage
    %        (wells, BC, source terms, compositional transport, aqueous
    %        tracer transport, and -- if enabled -- spatial bacterial
    %        transport) with all biological reaction source terms and the
    %        PHREEQC equilibrium call disabled.
    %     2. Partitions dt into nSub = ceil(dt/reactionSubstepMaxDt) equal
    %        substeps (summing exactly to dt) and, for each substep,
    %        solves a LOCAL biochemical reaction stage (MRST's bacterial
    %        Monod kinetics and aqueous tracer reaction sources) with all
    %        spatial flux divergence and external forcing (wells/BC/src)
    %        disabled. Each substep is followed by a PHREEQC equilibrium
    %        call that updates speciation, minerals, and the chemistry
    %        feedback (dolomite-buffered DIC/pH) that throttles the
    %        carbon-consuming Monod kinetics.
    %
    %   This is a flow-reaction operator split, not the pressure-transport
    %   split implemented by autodiff/sequential's
    %   SequentialPressureTransportModel, but it reuses that model's
    %   composition pattern: this meta-model overrides `stepFunction` and
    %   sets `stepFunctionIsLinear = true` so the outer NonLinearSolver
    %   calls it exactly once per control step, with NO outer Picard
    %   iterations. The expensive global flow solve therefore happens
    %   exactly once per control step; the reaction substeps are purely
    %   local (no spatial coupling between cells).
    %
    %   The `flowStageModel`/`reactionStageModel` sub-models used
    %   internally are plain `BiochemistryPhreeqcModel` instances (NOT
    %   `SequentialBiochemistryPhreeqcModel`), so their inner
    %   `NonLinearSolver.solveTimestep` calls use the ordinary inherited
    %   Newton `stepFunction` -- this avoids any risk of the nested solves
    %   recursively re-entering this class's overridden `stepFunction`
    %   (MRST models are value classes, so merely copying `model` would
    %   otherwise preserve its class and recurse).
    %
    % REQUIRED PARAMETERS:
    %   G, rock, fluid, compFluid, biochemFluid, includeWater, backend -
    %       Identical to BiochemistryPhreeqcModel.
    %
    % OPTIONAL PARAMETERS:
    %   reactionSubstepMaxDt    - Maximum duration of one local
    %                             reaction/PHREEQC substep. Required: must
    %                             be a positive, finite scalar. Any flow
    %                             timestep dt is split into
    %                             nSub = ceil(dt/reactionSubstepMaxDt)
    %                             equal substeps of size dt/nSub.
    %   flowNonLinearSolver     - NonLinearSolver used for the flow-stage
    %                             global solve (default: NonLinearSolver()).
    %   reactionNonLinearSolver - Retained for API compatibility. Local
    %                             reaction substeps use the bounded
    %                             positivity-preserving integrator.
    %   Any other property/value pair accepted by BiochemistryPhreeqcModel.
    %
    % NOTE:
    %   phreeqcBackend must be 'sequential-h2biochem-phreeqc' and
    %   phreeqcTimestepCoupling must be true: the split is only meaningful
    %   for the MRST-Monod-kinetics + sequential-PHREEQC-equilibrium
    %   hybrid backend that simulateSequentialH2BiochemPhreeqc's outer
    %   Picard driver also targets. Unlike that driver, this model uses no
    %   outer Picard iterations: it advances the flow and reaction stages
    %   exactly once each per control step (with nSub inner reaction
    %   substeps).
    %
    % SEE ALSO:
    %   BiochemistryPhreeqcModel, simulateSequentialH2BiochemPhreeqc,
    %   runSequentialH2BiochemPhreeqcEquilibrium,
    %   convertToSequentialBiochemistryPhreeqcModel,
    %   autodiff/sequential/models/SequentialPressureTransportModel

    properties
        % Maximum duration of one local reaction/PHREEQC substep. Any
        % control-step timestep dt is split into
        % nSub = ceil(dt/reactionSubstepMaxDt) equal substeps.
        reactionSubstepMaxDt

        % Plain BiochemistryPhreeqcModel sub-models, built internally by
        % buildStageModels (called from the constructor and from
        % validateModel, so a factory helper that further customizes the
        % meta-model after construction and then calls validateModel --
        % as simulateScheduleAD always does -- still gets sub-models that
        % match the final configuration).
        flowStageModel
        reactionStageModel

        % NonLinearSolver instances used for the flow-stage global solve
        % and for each local reaction substep, respectively.
        flowNonLinearSolver
        reactionNonLinearSolver
    end

    methods
        %-----------------------------------------------------------------%
        function model = SequentialBiochemistryPhreeqcModel(G, rock, fluid, compFluid, biochemFluid, includeWater, backend, varargin)
            baseOptions = getSequentialBiochemistryModelOptions(varargin);
            model = model@BiochemistryPhreeqcModel(G, rock, fluid, compFluid, ...
                biochemFluid, includeWater, backend, baseOptions{:});
            model = merge_options(model, varargin{:});

            % No outer Picard/nonlinear iterations: stepFunction is called
            % exactly once per control step (see NonLinearSolver.solveMinistep).
            model.stepFunctionIsLinear = true;

            validateattributes(model.reactionSubstepMaxDt, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'positive'}, mfilename, ...
                'reactionSubstepMaxDt');

            assert(model.phreeqcTimestepCoupling && model.isSequentialH2BiochemPhreeqcBackend(), ...
                ['SequentialBiochemistryPhreeqcModel requires phreeqcBackend=', ...
                 '''sequential-h2biochem-phreeqc'' and phreeqcTimestepCoupling=true: ', ...
                 'the coarse-flow/local-reaction split is only meaningful for the ', ...
                 'MRST-Monod-kinetics + sequential-PHREEQC-equilibrium hybrid backend.']);

            if isempty(model.flowNonLinearSolver)
                model.flowNonLinearSolver = NonLinearSolver();
            end
            assert(isa(model.flowNonLinearSolver, 'NonLinearSolver'), ...
                'flowNonLinearSolver must be a NonLinearSolver instance.');
            model.flowNonLinearSolver.identifier = 'SEQ-BIOCHEM-FLOW';

            if isempty(model.reactionNonLinearSolver)
                model.reactionNonLinearSolver = NonLinearSolver();
            end
            assert(isa(model.reactionNonLinearSolver, 'NonLinearSolver'), ...
                'reactionNonLinearSolver must be a NonLinearSolver instance.');
            model.reactionNonLinearSolver.identifier = 'SEQ-BIOCHEM-REACTION';

            % Build once at construction time for direct/simple usage;
            % validateModel rebuilds unconditionally so any further
            % post-construction customization (e.g. by
            % convertToSequentialBiochemistryPhreeqcModel) is picked up
            % before the first simulated step.
            model = model.buildStageModels();
        end

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            % Validate as a regular BiochemistryPhreeqcModel, then (re)build
            % the internal flow-stage/reaction-stage sub-models so they
            % stay in sync with any option changed after construction (see
            % buildStageModels).
            model = validateModel@BiochemistryPhreeqcModel(model, varargin{:});
            model = model.buildStageModels();
        end

        %-----------------------------------------------------------------%
        function model = buildStageModels(model)
            % (Re)build the internal flow-stage and reaction-stage
            % sub-models from the meta-model's CURRENT configuration.
            %
            % Plain BiochemistryPhreeqcModel instances are used (not
            % SequentialBiochemistryPhreeqcModel) so that their inner
            % NonLinearSolver.solveTimestep calls use the ordinary
            % inherited Newton stepFunction -- reusing this class would
            % recurse back into this very stepFunction, since MRST models
            % are value classes and copying "model" preserves its class.
            flow = BiochemistryPhreeqcModel(model.G, model.rock, model.fluid, ...
                model.compFluid, model.biochemFluid, model.water, model.AutoDiffBackend);
            flow = copyBiochemistryModelProperties(model, flow);
            % Flow stage: full spatial transport (compositional, aqueous
            % tracers, and -- when enabled -- bacterial diffusion), wells,
            % BC/src, no biological reaction sources, no PHREEQC call.
            flow.reactionsEnabled      = false;
            flow.localReactionMode     = false;
            flow.sequentialSplitActive = true;
            % CarbonLimitedGrowthRate requires the pre-step aqueous
            % carbon inventory (state.carbonSubstrateDt/Moles) that
            % getModelEquations only populates when reactionsEnabled is
            % true (see BiochemistryPhreeqcModel.getModelEquations).
            % Requesting it as post-convergence diagnostic output on the
            % (reaction-disabled) flow stage would therefore error; it is
            % meaningful only on the reaction stage.
            flow.OutputStateFunctions = setdiff(flow.OutputStateFunctions, ...
                {'CarbonLimitedGrowthRate'});
            % copyBiochemistryModelProperties copies the meta-model's
            % FacilityModel verbatim. If buildStageModels runs after the
            % meta-model's own FacilityModel has already been populated
            % with real wells (e.g. simulateScheduleAD calls
            % model.validateModel(fstruct, ...) once up front, which
            % cascades into FacilityModel.setupWells on the OUTER model
            % before this method runs), that stale well count would leak
            % into the stage sub-model. Reset to empty so each stage
            % (re)establishes its own well count from scratch, fresh, the
            % first time it is actually driven inside stepFunction --
            % avoiding FacilityModel's "number of wells has changed"
            % assertion when the inherited count does not match what
            % stepFunction subsequently passes in.
            flow.FacilityModel = [];
            flow = flow.validateModel();

            reaction = BiochemistryPhreeqcModel(model.G, model.rock, model.fluid, ...
                model.compFluid, model.biochemFluid, model.water, model.AutoDiffBackend);
            reaction = copyBiochemistryModelProperties(model, reaction);
            % Reaction stage: biological component/biomass/tracer
            % sources, no spatial flux, no wells/BC/src forcing. Reset
            % FacilityModel for the same reason as the flow stage above;
            % the reaction stage must always end up with zero wells, and
            % its own updateForChangedControls call inside stepFunction
            % (with empty forces.W) only succeeds if it starts from an
            % empty/zero-well FacilityModel rather than one inherited
            % from the meta-model.
            reaction.reactionsEnabled      = true;
            reaction.localReactionMode     = true;
            reaction.sequentialSplitActive = true;
            reaction.FacilityModel = [];
            reaction = reaction.validateModel();

            model.flowStageModel     = flow;
            model.reactionStageModel = reaction;
        end

        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok<INUSD>
            % The flow-stage and reaction-stage sub-models each already
            % ran their own updateAfterConvergence exactly once per
            % accepted solve inside stepFunction (once for the single
            % flow solve, once per accepted reaction substep) -- this is
            % where H2/reaction accounting, the PHREEQC equilibrium
            % coupling guard, and sulfate-lag caching for this timestep
            % are actually performed (see BiochemistryPhreeqcModel).
            %
            % Re-running the inherited BiochemistryPhreeqcModel
            % implementation here on the meta-model itself would be both
            % redundant and incorrect: it would evaluate FacilityModel
            % state functions (e.g. PhaseFlux) using THIS model's own
            % FacilityModel (populated with the real wells, for
            % reporting/output only -- see buildStageModels) against a
            % state produced by the well-less reaction stage, which
            % never populates the matching state.FacilityState. No
            % additional physics update is required at this level.
            % Mirrors SequentialPressureTransportModel's pattern in
            % autodiff/sequential.
            report = [];
        end

        %-----------------------------------------------------------------%
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nls, iteration, varargin) %#ok<INUSD>
            % Override of the standard Newton stepFunction, implementing
            % the coarse-flow/local-reaction operator split in two stages
            % run once per control step (stepFunctionIsLinear = true, so
            % the outer NonLinearSolver calls this exactly once, with no
            % outer Newton/Picard iteration at this level -- each stage
            % still iterates internally to its own convergence):
            %   1. One global flow+transport solve on flowStageModel
            %      (reactionsEnabled = false: no biological source terms).
            %   2. nSub local reaction substeps on reactionStageModel
            %      (localReactionMode = true: no spatial flux/diffusion),
            %      each closed by a PHREEQC equilibrium call that refreshes
            %      the chemistry feedback entering the next substep's
            %      kinetics.
            % See buildStageModels for why the two stage models are plain
            % BiochemistryPhreeqcModel instances rather than instances of
            % this class.
            timer = tic();
            flowModel     = model.flowStageModel;
            reactionModel = model.reactionStageModel;
            assert(~isempty(flowModel) && ~isempty(reactionModel), ...
                'SequentialBiochemistryPhreeqcModel:StageModelsNotBuilt', ...
                ['Flow-stage/reaction-stage sub-models have not been built. ', ...
                 'Call model.validateModel() before solving a timestep ', ...
                 '(simulateScheduleAD does this automatically).']);

            %% Stage 1: single global flow + compositional transport solve
            forceArg = flowModel.getDrivingForces(drivingForces);
            flowState0 = flowModel.validateState(state0);
            [flowModel, flowState0] = flowModel.updateForChangedControls(flowState0, drivingForces);
            [flowState, flowReport] = model.flowNonLinearSolver.solveTimestep( ...
                flowState0, dt, flowModel, 'initialGuess', state, forceArg{:});
            flowOk = flowReport.Converged || model.flowNonLinearSolver.continueOnFailure;

            stageReport = struct('FlowSolves', 1, 'FlowConverged', flowOk, ...
                'FlowStageReport', {flowReport}, ...
                'ReactionSubsteps', 0, 'ReactionSubstepSize', 0, ...
                'PhreeqcCalls', 0, ...
                'ReactionSubstepReports', {{}}, 'ReactionStageOk', false, ...
                'WallTime', toc(timer));

            if ~flowOk
                report = model.makeStepReport('Converged', false, 'Failure', true, ...
                    'FailureMsg', 'Coarse flow/transport stage failed to converge.');
                stageReport.WallTime = toc(timer);
                report.SequentialBiochemistryPhreeqcStages = stageReport;
                return
            end

            %% Stage 2: local biochemical substeps + one PHREEQC update
            % nSub equal substeps summing exactly to dt.
            nSub  = max(1, ceil(dt/model.reactionSubstepMaxDt));
            subDt = dt/nSub;

            % The reaction stage has no wells (no external forcing), so
            % its FacilityModel.WellModels is empty. flowState.wellSol
            % still holds the flow stage's (possibly non-empty) well
            % solution; clearing it here lets reactionModel.validateState
            % rebuild a wellSol sized to the well-less reaction stage
            % instead of mismatching model.WellModels in
            % updateWellSolAfterStep. The flow-stage wellSol is restored
            % onto the final state below.
            flowStateForReaction = flowState;
            flowStateForReaction.wellSol = [];
            reactionState0 = reactionModel.validateState(flowStateForReaction);
            [reactionModel, reactionState0] = reactionModel.updateForChangedControls( ...
                reactionState0, reactionModel.getValidDrivingForces());
            chemistryFeedback = ...
                buildSequentialH2BiochemPhreeqcChemistryFeedback( ...
                reactionModel, reactionState0);
            reactionModel = ...
                reactionModel.setSequentialH2BiochemPhreeqcChemistryFeedback( ...
                chemistryFeedback);
            reactionModel.FacilityModel = ...
                reactionModel.FacilityModel.setReservoirModel(reactionModel);

            reactionState  = reactionState0;
            substepReports = cell(nSub, 1);
            reactionOk     = true;
            phreeqcCalls   = 0;
            for i = 1:nSub
                [reactionState, subReport] = ...
                    advanceSequentialBiochemistryLocalReaction( ...
                    reactionModel, reactionState, subDt);
                substepReports{i} = subReport;
                reactionOk = subReport.Converged;
                if ~reactionOk
                    break
                end
                % One PHREEQC equilibrium per substep: dolomite
                % dissolution (the carbon resupply for methanogens and
                % acetogens) and the dolomite-buffered DIC/pH that
                % throttles the carbon Monod factor must both be
                % refreshed at the substep scale. A coarse-step-stale
                % feedback lets the early high-DIC transient compound
                % through biomass growth and overpredicts methanogenic
                % consumption, while a coarse-step carbon refill starves
                % it.
                reactionState = runSequentialH2BiochemPhreeqcEquilibrium( ...
                    reactionModel, reactionState);
                phreeqcCalls = phreeqcCalls + 1;
                chemistryFeedback = ...
                    buildSequentialH2BiochemPhreeqcChemistryFeedback( ...
                    reactionModel, reactionState);
                reactionModel = ...
                    reactionModel.setSequentialH2BiochemPhreeqcChemistryFeedback( ...
                    chemistryFeedback);
                reactionModel.FacilityModel = ...
                    reactionModel.FacilityModel.setReservoirModel(reactionModel);
            end

            stageReport.ReactionSubsteps       = nSub;
            stageReport.ReactionSubstepSize     = subDt;
            stageReport.ReactionSubstepReports = substepReports(1:i);
            stageReport.ReactionStageOk        = reactionOk;
            stageReport.ReactionIntegrationMethod = 'bounded-local-ode';
            stageReport.ReactionRejectedSteps = 0;
            stageReport.PhreeqcCalls = phreeqcCalls;

            if ~reactionOk
                report = model.makeStepReport('Converged', false, 'Failure', true, ...
                    'FailureMsg', sprintf( ...
                    'Local reaction substep %d of %d failed to converge.', i, nSub));
                report.SequentialBiochemistryPhreeqcStages = stageReport;
                return
            end

            stageReport.WallTime = toc(timer);

            % Wells are not part of the local reaction stage; restore the
            % flow-stage well solution for reporting/well-output purposes.
            reactionState.wellSol = flowState.wellSol;
            state = reactionState;

            report = model.makeStepReport('Converged', true, 'Failure', false);
            report.SequentialBiochemistryPhreeqcStages = stageReport;
        end
    end
end

function baseOptions = getSequentialBiochemistryModelOptions(options)
% Keep this class's own options away from the BiochemistryPhreeqcModel
% constructor.
ownOptions = {'reactionSubstepMaxDt', 'flowNonLinearSolver', 'reactionNonLinearSolver'};
assert(mod(numel(options), 2) == 0, ...
    'SequentialBiochemistryPhreeqcModel options must be property/value pairs.');
baseOptions = {};
for i = 1:2:numel(options)
    if ~any(strcmpi(options{i}, ownOptions))
        baseOptions(end + 1:end + 2) = options(i:i + 1); %#ok<AGROW>
    end
end
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
