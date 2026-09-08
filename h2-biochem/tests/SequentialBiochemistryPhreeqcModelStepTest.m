classdef SequentialBiochemistryPhreeqcModelStepTest < matlab.unittest.TestCase
% Timestep partitioning/report and physical-state checks for one step of
% SequentialBiochemistryPhreeqcModel.
%
% A single arbitrary (i.e. NOT a clean multiple of the maximum reaction
% substep duration) control-step timestep is solved once in
% TestClassSetup; every Test method below only inspects the resulting
% report/state, so no test method contains control-flow logic of its own.
%
% This requires a registered IPhreeqcCOM server and an on-disk
% PHREEQC_Modified.DAT (both are needed because every accepted coarse
% step runs one real PHREEQC equilibrium call); the whole class
% is skipped cleanly via assumeTrue in TestClassSetup when either is
% unavailable.
%
% SEE ALSO:
%   SequentialBiochemistryPhreeqcModel, convertToSequentialBiochemistryPhreeqcModel

    properties
        model
        state0
        dt
        maxDt
        resultState
        report
    end

    methods (TestClassSetup)
        function setupModelAndRunOneStep(testCase)
            mrstModule add ad-core ad-props compositional deckformat h2-biochem
            databaseFile = resolveStepTestDatabaseFile();
            testCase.assumeTrue(~isempty(databaseFile) && isIPhreeqcComAvailable(), ...
                ['A registered IPhreeqcCOM server and an on-disk ', ...
                 'PHREEQC_Modified.DAT are required for this test.']);

            % Small grid and a short single control period keep this test
            % fast; options otherwise match the high-rate sequential
            % h2-biochem-PHREEQC benchmark setup (see
            % runThreeBackendComparison / exampleSequentialBiochemistryPhreeqc1D).
            [~, baseModel, schedule, state0] = setupH2StorageExampleWithSRB_benchmark( ...
                'rate', 'highrate', ...
                'scheduleMode', 'injection', ...
                'gridCells', 4, ...
                'domainLength', 5, ...
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

            testCase.maxDt = 0.4*day;
            % Deliberately NOT a clean multiple of maxDt, to prove the
            % substep partitioning is general rather than hardcoded to
            % 2 day / 0.4 day steps.
            testCase.dt = 0.65*day;
            testCase.model = convertToSequentialBiochemistryPhreeqcModel(baseModel, ...
                'reactionSubstepMaxDt', testCase.maxDt);
            testCase.state0 = testCase.model.validateState(state0);
            forces = schedule.control(1);

            [testCase.resultState, testCase.report] = testCase.model.stepFunction( ...
                testCase.state0, testCase.state0, testCase.dt, forces, [], [], 1);
        end
    end

    methods (Test)
        %% Report structure: one flow solve, expected reaction substeps.
        function testStepConverged(testCase)
            testCase.verifyTrue(testCase.report.Converged);
        end

        function testExactlyOneFlowSolve(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyEqual(stages.FlowSolves, 1);
        end

        function testFlowStageConverged(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyTrue(stages.FlowConverged);
        end

        function testExpectedReactionSubstepCount(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyEqual(stages.ReactionSubsteps, ...
                ceil(testCase.dt/testCase.maxDt));
        end

        function testReactionSubstepsAreEqualAndSumToDt(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            total = stages.ReactionSubsteps*stages.ReactionSubstepSize;
            testCase.verifyEqual(total, testCase.dt, 'AbsTol', 1e-9*testCase.dt);
        end

        function testReactionSubstepReportCountMatchesSubstepCount(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyEqual(numel(stages.ReactionSubstepReports), ...
                stages.ReactionSubsteps);
        end

        function testReactionStageConverged(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyTrue(stages.ReactionStageOk);
        end

        function testBoundedIntegratorRejectsNoSubsteps(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyEqual(stages.ReactionIntegrationMethod, ...
                'bounded-local-ode');
            testCase.verifyEqual(stages.ReactionRejectedSteps, 0);
        end

        function testOnePhreeqcCallPerFlowStep(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyEqual(stages.PhreeqcCalls, 1);
        end

        function testStageTimingIsReported(testCase)
            stages = testCase.report.SequentialBiochemistryPhreeqcStages;
            testCase.verifyGreaterThanOrEqual(stages.WallTime, 0);
        end

        %% Flow stage: full transport, no biological sources, no PHREEQC.
        function testFlowStageHasReactionsDisabled(testCase)
            testCase.verifyFalse(testCase.model.flowStageModel.reactionsEnabled);
        end

        function testFlowStageHasSpatialTransportEnabled(testCase)
            testCase.verifyFalse(testCase.model.flowStageModel.localReactionMode);
        end

        %% Reaction stage: local biological sources, no spatial transport.
        function testReactionStageHasReactionsEnabled(testCase)
            testCase.verifyTrue(testCase.model.reactionStageModel.reactionsEnabled);
        end

        function testReactionStageHasSpatialTransportDisabled(testCase)
            testCase.verifyTrue(testCase.model.reactionStageModel.localReactionMode);
        end

        %% Physical state after the step.
        function testComponentsAreFiniteAndNonnegative(testCase)
            z = value(testCase.resultState.components);
            testCase.verifyTrue(all(isfinite(z(:))));
            testCase.verifyGreaterThanOrEqual(z(:), 0);
        end

        function testPressureIsFiniteAndPositive(testCase)
            p = value(testCase.resultState.pressure);
            testCase.verifyTrue(all(isfinite(p(:))));
            testCase.verifyGreaterThan(p(:), 0);
        end

        function testBiomassIsFiniteAndNonnegative(testCase)
            nbact = value(testCase.resultState.nbact);
            testCase.verifyTrue(all(isfinite(nbact(:))));
            testCase.verifyGreaterThanOrEqual(nbact(:), 0);
        end

        function testCumulativeH2ConsumptionIsFiniteAndNonnegative(testCase)
            cumulative = value( ...
                testCase.resultState.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles);
            testCase.verifyTrue(all(isfinite(cumulative(:))));
            testCase.verifyGreaterThanOrEqual(cumulative(:), 0);
        end
    end
end

function databaseFile = resolveStepTestDatabaseFile()
% Best-effort PHREEQC_Modified.DAT lookup. Returns '' (rather than
% throwing) when the database cannot be located, so calling tests can
% gate on it via assumeTrue and skip cleanly.
    databaseFile = getenv('PHREEQC_DATABASE_FILE');
    if isempty(strtrim(databaseFile))
        databaseFile = which('PHREEQC_Modified.DAT');
    end
    if isempty(strtrim(databaseFile)) || ~isfile(databaseFile)
        databaseFile = '';
    end
end

function tf = isIPhreeqcComAvailable()
% Best-effort probe for a registered IPhreeqcCOM server. Never throws:
% returns false on any failure (non-Windows, unregistered server, ...) so
% calling tests can gate on it via assumeTrue and skip cleanly.
    tf = false;
    if ~ispc
        return
    end
    try
        server = actxserver('IPhreeqcCOM.Object');
        delete(server);
        tf = true;
    catch
        tf = false;
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
