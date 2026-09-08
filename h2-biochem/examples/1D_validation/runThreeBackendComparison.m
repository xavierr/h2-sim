function scenarios = runThreeBackendComparison(varargin)
% Compare the original, compositional-PHREEQC, and coarse-flow/local-
% reaction sequential backends.
%
% Aqueous H2 mole fraction and total dissolved inorganic carbon (C(4)
% molality, mol/kgw) profiles are shown at the end of injection, storage,
% and production. DIC is compared instead of the aqueous CO2 mole fraction
% because CO2(aq) is slaved to pH and alkalinity, so the coupled backends
% and the UGFACT solution only agree on the same quantity when the full
% carbonate speciation is summed. Additional figures compare H2 loss over
% time and its final spatial distribution. Results remain in memory and no
% files are written.
%
% The three compared backends are:
%   1. Original h2-biochem (BiochemistryModel, no PHREEQC).
%   2. Sequential compositional PHREEQC (BiochemistryPhreeqcModel,
%      phreeqcBackend='sequential-compositional-phreeqc'): PHREEQC owns
%      MET/ACE/SRB kinetics entirely, applied once per transport timestep.
%   3. Sequential (SequentialBiochemistryPhreeqcModel): MRST retains
%      microbial kinetics/biomass, split into one coarse global flow step
%      followed by local reaction substeps each closed by a PHREEQC
%      equilibrium update -- see convertToSequentialBiochemistryPhreeqcModel.
%
% The older outer-Picard hybrid driver (utils/simulateSequentialH2BiochemPhreeqc.m,
% phreeqcBackend='sequential-h2biochem-phreeqc' run without the coarse-flow
% split) is intentionally not part of this comparison for now.
%
% EXAMPLE:
%   db = '\\wsl.localhost\Ubuntu\path\to\PHREEQC_Modified.DAT';
%   scenarios = runThreeBackendComparison('phreeqcDatabaseFile', db);
%   scenarios = runThreeBackendComparison('phreeqcDatabaseFile', db, ...
%       'referenceUseSoreideWhitsonEOS', true);

    mrstModule add ad-props compositional deckformat h2-biochem

    opt = struct( ...
        'phreeqcDatabaseFile', '', ...
        'flowTimestepMaxDt', 2*day, ...
        'reactionSubstepMaxDt', 0.4*day, ...
        'referenceUseSoreideWhitsonEOS', false);
    opt = merge_options(opt, varargin{:});
    validateattributes(opt.referenceUseSoreideWhitsonEOS, ...
        {'logical', 'numeric'}, {'scalar', 'real', 'finite'}, ...
        mfilename, 'referenceUseSoreideWhitsonEOS');
    opt.referenceUseSoreideWhitsonEOS = ...
        logical(opt.referenceUseSoreideWhitsonEOS);
    validateattributes(opt.flowTimestepMaxDt, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'flowTimestepMaxDt');
    validateattributes(opt.reactionSubstepMaxDt, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'reactionSubstepMaxDt');
    databaseFile = resolvePhreeqcDatabaseFile(opt.phreeqcDatabaseFile);

    commonOptions = { ...
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
        'initialOverallCO2', 0};

    cases = [ ...
        makeCase("MRST-h2biochem", { ...
            'equilibrateInitialCO2', false, ...
            'phreeqcTimestepCoupling', false}, 'plain'), ...
        makeCase("MRST-compositional PHREEQC", { ...
            'equilibrateInitialCO2', false, ...
            'phreeqcTimestepCoupling', true, ...
            'phreeqcBackend', 'sequential-compositional-phreeqc', ...
            'phreeqcDatabaseFile', databaseFile}, 'plain'), ...
        makeCase("MRST-h2biochem PHREEQC", { ...
            'equilibrateInitialCO2', false, ...
            'paperBiomassKinetics', true, ...
            'nbact0', 1, ...
            'phreeqcTimestepCoupling', true, ...
            'phreeqcBackend', 'sequential-h2biochem-phreeqc', ...
            'phreeqcDatabaseFile', databaseFile}, 'coarseSplit')];

    scenarios = initializeScenarios(numel(cases));
    for caseNo = 1:numel(cases)
        current = cases(caseNo);
        fprintf('\nBackend %d/%d: %s\n', ...
            caseNo, numel(cases), current.name);
        [~, model, schedule, state0] = ...
            setupH2StorageExampleWithSRB_benchmark( ...
            commonOptions{:}, current.options{:});

        if strcmp(current.driver, 'coarseSplit')
            assert(model.bacterialDecayOrder == 1, ...
               'The UGFACT-matched sequential case requires first-order biomass decay.');
            assert(max(abs(value(state0.nbact) - 1), [], 'all') < 1e-12, ...
               'The UGFACT-matched sequential case requires initial N/N0 = 1.');
            schedule = coarsenScheduleByControl(schedule, opt.flowTimestepMaxDt);
            model = convertToSequentialBiochemistryPhreeqcModel(model, ...
                'reactionSubstepMaxDt', opt.reactionSubstepMaxDt);
        end

        solver = NonLinearSolver();
        solver.maxTimestepCuts = 12;
        timer = tic();
        [ws, states, report] = simulateScheduleAD( ...
            state0, model, schedule, 'nonlinearSolver', solver);
        if isstruct(report) && isfield(report, 'Failure')
            assert(~report.Failure, ...
                'Simulation report indicates failure for "%s".', current.name);
        end

        metrics = collectMetrics(states, schedule, model);
        scenarios(caseNo).name = current.name;
        scenarios(caseNo).model = model;
        scenarios(caseNo).schedule = schedule;
        scenarios(caseNo).states = states;
        scenarios(caseNo).ws = ws;
        scenarios(caseNo).state0 = state0;
        scenarios(caseNo).timeDays = metrics.timeDays;
        scenarios(caseNo).snapshotIndices = metrics.snapshotIndices;
        scenarios(caseNo).snapshotLabels = metrics.snapshotLabels;
        scenarios(caseNo).aqueousH2 = metrics.aqueousH2;
        scenarios(caseNo).aqueousDIC = metrics.aqueousDIC;
        scenarios(caseNo).lossPercent = metrics.lossPercent;
        scenarios(caseNo).finalSpatialConsumptionMoles = ...
            metrics.finalSpatialConsumptionMoles;
        scenarios(caseNo).xDimensionless = metrics.xDimensionless;
        scenarios(caseNo).injectedH2Moles = metrics.injectedH2Moles;
        scenarios(caseNo).consumedH2Moles = metrics.consumedH2Moles;
        scenarios(caseNo).runtimeSeconds = toc(timer);
    end

    ugfactRoot = fileparts(fileparts(fileparts(databaseFile)));
    scenarios(end + 1) = runUGFACTReference( ...
        ugfactRoot, opt.referenceUseSoreideWhitsonEOS);
    printSummary(scenarios);
    createComparisonFigures(scenarios);
end

function item = makeCase(name, options, driver)
    % driver: 'plain' (simulateScheduleAD) or 'coarseSplit' (coarse-flow/
    % local-reaction split via convertToSequentialBiochemistryPhreeqcModel).
    item = struct('name', name, 'options', {options}, 'driver', driver);
end

function scenarios = initializeScenarios(nCases)
    template = struct( ...
        'name', "", 'model', [], 'schedule', [], 'states', [], ...
        'ws', [], 'state0', [], ...
        'timeDays', [], 'snapshotIndices', [], 'snapshotLabels', [], ...
        'aqueousH2', [], 'aqueousDIC', [], 'lossPercent', [], ...
        'finalSpatialConsumptionMoles', [], 'xDimensionless', [], ...
        'injectedH2Moles', nan, 'consumedH2Moles', nan, ...
        'runtimeSeconds', nan);
    scenarios = repmat(template, nCases, 1);
end

function metrics = collectMetrics(states, schedule, model)
    assert(numel(states) == numel(schedule.step.val), ...
        'Saved states and schedule must contain the same number of steps.');
    componentNames = model.EOSModel.getComponentNames();
    idxH2 = findComponentIndex(componentNames, {'H2', 'Hydrogen'});
    idxCO2 = findComponentIndex(componentNames, {'CO2', 'CarbonDioxide'});

    controls = schedule.step.control(:);
    snapshotIndices = [find(controls == 1, 1, 'last'), ...
        find(controls == 2, 1, 'last'), numel(states)];
    assert(all(snapshotIndices > 0), ...
        'The comparison requires injection, storage, and production periods.');
    snapshotLabels = ["End injection", "End storage", "End simulation"];

    aqueousH2 = zeros(model.G.cells.num, numel(snapshotIndices));
    aqueousDIC = zeros(model.G.cells.num, numel(snapshotIndices));
    for i = 1:numel(snapshotIndices)
        state = states{snapshotIndices(i)};
        aqueousH2(:, i) = componentColumn(state.x, idxH2);
        aqueousDIC(:, i) = aqueousDICColumn(state, model, componentNames, idxCO2);
    end

    nReactions = model.biochemFluid.nbioreact;
    cumulative = zeros(numel(states), nReactions);
    finalSpatial = zeros(model.G.cells.num, nReactions);
    for reaction = 1:nReactions
        [~, reactionCumulative] = computeH2Consumption( ...
            states, schedule, model, reaction);
        cumulative(:, reaction) = sum(reactionCumulative, 1).';
        finalSpatial(:, reaction) = reactionCumulative(:, end);
    end
    totalConsumed = sum(cumulative, 2);
    injected = prescribedInjectedH2(schedule, model);
    x = model.G.cells.centroids(:, 1);
    x = (x - min(x))./max(max(x) - min(x), eps);

    metrics = struct( ...
        'timeDays', cumsum(schedule.step.val(:))./day, ...
        'snapshotIndices', snapshotIndices, ...
        'snapshotLabels', snapshotLabels, ...
        'aqueousH2', aqueousH2, ...
        'aqueousDIC', aqueousDIC, ...
        'lossPercent', 100.*totalConsumed./injected, ...
        'finalSpatialConsumptionMoles', sum(finalSpatial, 2), ...
        'xDimensionless', x, ...
        'injectedH2Moles', injected, ...
        'consumedH2Moles', totalConsumed(end));
end

function scenario = runUGFACTReference(ugfactRoot, useSoreideWhitsonEOS)
    referenceFile = fullfile(ugfactRoot, 'examples', 'H2Storage1D.m');
    assert(isfile(referenceFile), ...
        'UGFACT reference driver was not found at %s.', referenceFile);
    mrstPath('register', 'UGFACT2', ugfactRoot);

    source = fileread(referenceFile);
    assert(contains(source, 'rateCase = ''highrate'''), ...
        'UGFACT H2Storage1D must explicitly select the high-rate case.');
    eosExpression = 'useSoreideWhitsonEOS\s*=\s*(true|false)\s*;';
    assert(~isempty(regexp(source, eosExpression, 'once')), ...
        ['The UGFACT driver EOS selector changed. Update the reference ', ...
         'override before running this comparison.']);
    source = strrep(source, 'clear;clc;close all', '');
    eosValue = char(string(useSoreideWhitsonEOS));
    source = regexprep(source, eosExpression, ...
        ['useSoreideWhitsonEOS = ', eosValue, ';'], 'once');

    if useSoreideWhitsonEOS
        eosName = "Soreide-Whitson EOS";
    else
        eosName = "original UGFACT EOS";
    end
    fprintf('\nBackend 4/4: UGFACT H2Storage1D reference (high rate, %s)\n', ...
        eosName);
    timer = tic();
    evalin('base', source);
    rateCase = evalin('base', 'rateCase');
    model = evalin('base', 'model');
    schedule = evalin('base', 'schedule');
    states = evalin('base', 'states');
    state0 = evalin('base', 'state0');
    wellSol = evalin('base', 'wellSol');
    assert(strcmp(rateCase, 'highrate'), ...
        'UGFACT run did not use the high-rate kinetic case.');
    assert(model.Kinetic.mu_MET == 4.1 && ...
        model.Kinetic.mu_ACE == 1.9 && model.Kinetic.mu_SRB == 5.5, ...
        'UGFACT high-rate kinetic constants are incorrect.');
    actualEOSChoice = evalin('base', 'useSoreideWhitsonEOS');
    assert(logical(actualEOSChoice) == useSoreideWhitsonEOS, ...
        'UGFACT run did not use the requested EOS.');

    metrics = collectUGFACTMetrics(states, schedule, model);
    scenario = initializeScenarios(1);
    scenario.name = "UGFACT";
    scenario.model = model;
    scenario.schedule = schedule;
    scenario.states = states;
    scenario.ws = wellSol;
    scenario.state0 = state0;
    scenario.timeDays = metrics.timeDays;
    scenario.snapshotIndices = metrics.snapshotIndices;
    scenario.snapshotLabels = metrics.snapshotLabels;
    scenario.aqueousH2 = metrics.aqueousH2;
    scenario.aqueousDIC = metrics.aqueousDIC;
    scenario.lossPercent = metrics.lossPercent;
    scenario.finalSpatialConsumptionMoles = ...
        metrics.finalSpatialConsumptionMoles;
    scenario.xDimensionless = metrics.xDimensionless;
    scenario.injectedH2Moles = metrics.injectedH2Moles;
    scenario.consumedH2Moles = metrics.consumedH2Moles;
    scenario.runtimeSeconds = toc(timer);
end

function metrics = collectUGFACTMetrics(states, schedule, model)
    assert(numel(states) == numel(schedule.step.val), ...
        'UGFACT states and schedule must contain the same number of steps.');
    componentNames = model.EOSModel.getComponentNames();
    idxH2 = findComponentIndex(componentNames, {'H2', 'Hydrogen'});
    idxCO2 = findComponentIndex(componentNames, {'CO2', 'CarbonDioxide'});

    controls = schedule.step.control(:);
    snapshotIndices = [find(controls == 1, 1, 'last'), ...
        find(controls == 2, 1, 'last'), numel(states)];
    snapshotLabels = ["End injection", "End storage", "End simulation"];
    aqueousH2 = zeros(model.G.cells.num, 3);
    aqueousDIC = zeros(model.G.cells.num, 3);
    for i = 1:3
        state = states{snapshotIndices(i)};
        aqueousH2(:, i) = componentColumn(state.x, idxH2);
        aqueousDIC(:, i) = aqueousDICColumn(state, model, componentNames, idxCO2);
    end

    cumulative = zeros(numel(states), 1);
    finalSpatial = zeros(model.G.cells.num, 1);
    for step = 1:numel(states)
        state = states{step};
        solution = state.Solution;
        assert(all(isfield(solution, ...
            {'MET_Rate', 'ACE_Rate', 'SRB_Rate', 'Water'})), ...
            'UGFACT state %d is missing reaction-rate diagnostics.', step);
        rateMolesPerDay = (solution.MET_Rate + solution.ACE_Rate + ...
            solution.SRB_Rate).*solution.Water./1000;
        increment = rateMolesPerDay.*schedule.step.val(step)./day;
        finalSpatial = finalSpatial + increment;
        cumulative(step) = sum(finalSpatial);
    end

    injected = prescribedInjectedH2(schedule, model);
    x = model.G.cells.centroids(:, 1);
    x = (x - min(x))./max(max(x) - min(x), eps);
    metrics = struct( ...
        'timeDays', cumsum(schedule.step.val(:))./day, ...
        'snapshotIndices', snapshotIndices, ...
        'snapshotLabels', snapshotLabels, ...
        'aqueousH2', aqueousH2, ...
        'aqueousDIC', aqueousDIC, ...
        'lossPercent', 100.*cumulative./injected, ...
        'finalSpatialConsumptionMoles', finalSpatial, ...
        'xDimensionless', x, ...
        'injectedH2Moles', injected, ...
        'consumedH2Moles', cumulative(end));
end

function dic = aqueousDICColumn(state, model, componentNames, idxCO2)
    % Total dissolved inorganic carbon (C(4)) molality, mol/kgw, so the
    % PHREEQC-coupled backends and the UGFACT solution are compared on the
    % same quantity instead of the pH-slaved aqueous CO2 mole fraction.
    if isfield(state, 'phreeqcTotalCarbon') && ~isempty(state.phreeqcTotalCarbon)
        dic = reshape(value(state.phreeqcTotalCarbon), [], 1);
        return;
    end
    if isfield(state, 'Solution') && isfield(state.Solution, 'C4')
        dic = reshape(value(state.Solution.C4), [], 1);
        return;
    end
    % Backends without an aqueous carbon model (original h2-biochem): sum the
    % transported carbonate tracer and the EOS-dissolved CO2.
    idxH2O = findComponentIndex(componentNames, {'H2O', 'Water'});
    xCO2 = componentColumn(state.x, idxCO2);
    xH2O = max(componentColumn(state.x, idxH2O), eps);
    co2Molality = xCO2 ./ xH2O .* 55.508;
    if isfield(state, 'tracerHCO3') && ~isempty(state.tracerHCO3)
        hco3Molality = reshape(value(state.tracerHCO3), [], 1) ./ 1000;
    else
        hco3Molality = zeros(model.G.cells.num, 1);
    end
    dic = hco3Molality + co2Molality;
end

function column = componentColumn(composition, componentIndex)
    if iscell(composition)
        column = value(composition{componentIndex});
    else
        column = value(composition(:, componentIndex));
    end
    column = column(:);
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

function printSummary(scenarios)
    names = [scenarios.name].';
    injected = [scenarios.injectedH2Moles].';
    consumed = [scenarios.consumedH2Moles].';
    loss = 100.*consumed./injected;
    runtime = [scenarios.runtimeSeconds].';
    summary = table(names, injected, consumed, loss, runtime, ...
        'VariableNames', {'Backend', 'InjectedH2_mol', 'ConsumedH2_mol', ...
        'FinalH2Loss_percent', 'Runtime_seconds'});
    disp(summary);
end

function createComparisonFigures(scenarios)
    colors = paperColors(numel(scenarios));
    labels = scenarios(1).snapshotLabels;

    fig = paperFigure([26, 16], ...
        'Aqueous H2 and dissolved inorganic carbon profiles');
    layout = tiledlayout(fig, 2, 3, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for component = 1:2
        for snapshot = 1:3
            ax = nexttile(layout);
            hold(ax, 'on');
            for backend = 1:numel(scenarios)
                if component == 1
                    profile = scenarios(backend).aqueousH2(:, snapshot);
                else
                    profile = scenarios(backend).aqueousDIC(:, snapshot);
                end
                plot(ax, scenarios(backend).xDimensionless, profile, ...
                    'LineWidth', 1.8, 'Color', colors(backend, :), ...
                    'DisplayName', char(scenarios(backend).name));
            end
            title(ax, labels(snapshot));
            if component == 1
                ylabel(ax, 'Aqueous H_2 mole fraction');
            else
                ylabel(ax, 'Dissolved inorganic carbon (mol kgw^{-1})');
            end
            xlabel(ax, 'Dimensionless distance');
            if component == 1 && snapshot == 3
                legend(ax, 'Location', 'northeast');
            end
            styleAxes(ax);
        end
    end
    paperExport(fig, 'three_backends_H2_DIC_profiles');

    fig = paperFigure([18, 11], 'H2 loss over time');
    ax = axes(fig);
    hold(ax, 'on');
    for backend = 1:numel(scenarios)
        plot(ax, scenarios(backend).timeDays, ...
            scenarios(backend).lossPercent, ...
            'LineWidth', 1.8, 'Color', colors(backend, :), ...
            'DisplayName', char(scenarios(backend).name));
    end
    xline(ax, 50, 'k:', 'End injection');
    xline(ax, 200, 'k:', 'End storage');
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Consumed injected H_2 (%)');
    legend(ax, 'Location', 'best');
    styleAxes(ax);
    paperExport(fig, 'three_backends_H2_loss_over_time');

    fig = paperFigure([18, 11], 'Spatial H2 consumption');
    ax = axes(fig);
    hold(ax, 'on');
    for backend = 1:numel(scenarios)
        plot(ax, scenarios(backend).xDimensionless, ...
            scenarios(backend).finalSpatialConsumptionMoles, ...
            'LineWidth', 1.8, 'Color', colors(backend, :), ...
            'DisplayName', char(scenarios(backend).name));
    end
    xlabel(ax, 'Dimensionless distance from injector');
    ylabel(ax, 'Cumulative H_2 consumed (mol cell^{-1})');
    legend(ax, 'Location', 'best');
    styleAxes(ax);
    paperExport(fig, 'three_backends_spatial_H2_consumption');
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
        error('runThreeBackendComparison:MissingPhreeqcDatabase', ...
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
%}
