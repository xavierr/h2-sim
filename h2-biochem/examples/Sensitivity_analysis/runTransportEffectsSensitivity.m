function results = runTransportEffectsSensitivity(varargin)
% Run a transport and bio-clogging sensitivity study and display its plots.
%
% The study is repeated for one or both reaction backends:
%   'original'  - MRST owns all reactions (BiochemistryModel), reduced
%                 fixed-pH CO2/HCO3 buffer, no PHREEQC. This is the first
%                 "Original h2-biochem" case from testAllDiffusionEffects.
%   'multirate' - Coarse global flow step followed by local reaction
%                 substeps, each closed by a PHREEQC equilibrium update
%                 (SequentialBiochemistryPhreeqcModel, built with
%                 convertToSequentialBiochemistryPhreeqcModel). MRST retains
%                 the microbial kinetics/biomass; PHREEQC owns the aqueous
%                 and mineral equilibrium. Requires PHREEQC_Modified.DAT.
%
% Each backend's cases are compared against that backend's own baseline.
% Cases that fail to converge are recorded and skipped rather than aborting
% the study (the coarse-flow/local-reaction split has only been exercised
% with the transport effects switched off, so some combinations may fail).
% Results remain in memory and figures are left open; no files are written.
%
% EXAMPLE:
%   % Original backend only (no PHREEQC needed):
%   results = runTransportEffectsSensitivity('backends', {'original'});
%
%   % Both backends:
%   db = '/path/to/PHREEQC_Modified.DAT';
%   results = runTransportEffectsSensitivity('phreeqcDatabaseFile', db);

    mrstModule add ad-props compositional deckformat h2-biochem

    opt = struct( ...
        'meaningfulLossThreshold', 0.1, ...
        'scheduleMode', 'complete', ...
        'rate', 'highrate', ...
        'backends', {{'original', 'multirate'}}, ...
        'phreeqcDatabaseFile', '', ...
        'flowTimestepMaxDt', 2*day, ...
        'reactionSubstepMaxDt', 0.4*day);
    opt = merge_options(opt, varargin{:});
    validateOptions(opt);

    backends = resolveBackends(opt);
    cases = transportCases();
    nCases = numel(cases);
    results = initializeResults(numel(backends)*nCases);

    commonOptions = { ...
        'rate', opt.rate, ...
        'scheduleMode', opt.scheduleMode, ...
        'injectionCO2', 0, ...
        'gridCells', 200, ...
        'carbonateBuffer', true, ...
        'carbonateBufferPH', 5.9, ...
        'initialHCO3', 4.4119e-2, ...
        'bacteriamodel', true};

    for b = 1:numel(backends)
        backend = backends(b);
        for caseNo = 1:nCases
            idx = (b - 1)*nCases + caseNo;
            fprintf('\n[%s] Case %d/%d: %s\n', backend.label, caseNo, ...
                nCases, cases(caseNo).name);
            setupOptions = [commonOptions, backend.setupOptions, ...
                cases(caseNo).options];
            results(idx).backend = string(backend.label);
            results(idx).caseName = string(cases(caseNo).name);
            results(idx).isBaseline = (caseNo == 1);
           % try
                [~, model, schedule, state0] = ...
                    setupH2StorageExampleWithSRB_benchmark(setupOptions{:});
                [model, schedule] = applyBackend(backend, model, schedule, ...
                    state0, opt);
                solver = NonLinearSolver();
                solver.maxTimestepCuts = 12;
                timer = tic();
                [~, states, report] = simulateScheduleAD( ...
                    state0, model, schedule, 'nonlinearSolver', solver);
                if isstruct(report) && isfield(report, 'Failure')
                    assert(~report.Failure, ...
                        'Simulation report indicates failure.');
                end
                metrics = computeMetrics(states, schedule, model);
                results(idx).completed = true;
                results(idx).runtimeSeconds = toc(timer);
                results(idx).timeDays = metrics.timeDays;
                results(idx).lossPercent = metrics.lossPercent;
                results(idx).finalLossPercent = metrics.finalLossPercent;
                results(idx).consumedH2Moles = metrics.consumedH2Moles;
                results(idx).injectedH2Moles = metrics.injectedH2Moles;
                results(idx).reactionTotalsMoles = metrics.reactionTotalsMoles;
                results(idx).finalSpatialConsumptionMoles = ...
                    metrics.finalSpatialConsumptionMoles;
                results(idx).xDimensionless = metrics.xDimensionless;
                results(idx).errorMessage = "";
            % catch err
            %     warning('runTransportEffectsSensitivity:caseFailed', ...
            %         '[%s] "%s" failed: %s', backend.label, ...
            %         cases(caseNo).name, err.message);
            %     results(idx).completed = false;
            %     results(idx).errorMessage = string(err.message);
            % end
        end
    end

    summary = buildSummary(results, opt);
    createPaperFigures(results, backends, cases, opt);
    disp(summary);
end

function validateOptions(opt)
    validateattributes(opt.meaningfulLossThreshold, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'}, ...
        mfilename, 'meaningfulLossThreshold');
    validateattributes(opt.flowTimestepMaxDt, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'flowTimestepMaxDt');
    validateattributes(opt.reactionSubstepMaxDt, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'reactionSubstepMaxDt');
end

function backends = resolveBackends(opt)
    keys = opt.backends;
    if ischar(keys) || isstring(keys)
        keys = cellstr(keys);
    end
    assert(iscellstr(keys) && ~isempty(keys), ...
        'backends must be a nonempty cell array of backend names.');
    keys = unique(lower(strtrim(keys)), 'stable');
    template = struct('key', '', 'label', '', 'setupOptions', {{}}, ...
        'coarseSplit', false);
    backends = repmat(template, numel(keys), 1);
    for i = 1:numel(keys)
        switch keys{i}
            case 'original'
                backends(i).key = 'original';
                backends(i).label = 'Original h2-biochem';
                backends(i).setupOptions = { ...
                    'paperBiomassKinetics', false, ...
                    'equilibrateInitialCO2', true, ...
                    'phreeqcTimestepCoupling', false};
                backends(i).coarseSplit = false;
            case 'multirate'
                db = resolvePhreeqcDatabaseFile(opt.phreeqcDatabaseFile);
                backends(i).key = 'multirate';
                backends(i).label = 'Multirate PHREEQC';
                backends(i).setupOptions = { ...
                    'paperBiomassKinetics', true, ...
                    'nbact0', 1, ...
                    'equilibrateInitialCO2', false, ...
                    'phreeqcTimestepCoupling', true, ...
                    'phreeqcBackend', 'sequential-h2biochem-phreeqc', ...
                    'phreeqcDatabaseFile', db};
                backends(i).coarseSplit = true;
            otherwise
                error(['Unknown backend "%s". Use ''original'' and/or ', ...
                    '''multirate''.'], keys{i});
        end
    end
end

function [model, schedule] = applyBackend(backend, model, schedule, state0, opt)
    if ~backend.coarseSplit
        assert(~model.phreeqcTimestepCoupling, ...
            'The original backend must not call PHREEQC.');
        assert(~model.isSequentialCompositionalPhreeqcBackend() && ...
            ~model.isSequentialH2BiochemPhreeqcBackend(), ...
            'The original backend must use MRST-owned reactions.');
        return;
    end
    assert(model.bacterialDecayOrder == 1, ...
        'The multirate backend requires first-order biomass decay.');
    assert(max(abs(value(state0.nbact) - 1), [], 'all') < 1e-12, ...
        'The multirate backend requires initial N/N0 = 1.');
    schedule = coarsenScheduleByControl(schedule, opt.flowTimestepMaxDt);
    model = convertToSequentialBiochemistryPhreeqcModel(model, ...
        'reactionSubstepMaxDt', opt.reactionSubstepMaxDt);
end

function cases = transportCases()
    % cases = [ ...
    %     makeCase("Baseline", false, false, false, false, false), ...
    %     makeCase("Microbial diffusion", true, false, false, false, false), ...
    %     makeCase("Chemotaxis", false, true, false, false, false), ...
    %     makeCase("Microbial diffusion + chemotaxis", ...
    %         true, true, false, false, false), ...
    %     makeCase("Molecular diffusion", false, false, true, false, false), ...
    %     makeCase("Mechanical dispersion", false, false, false, true, false), ...
    %     makeCase("Molecular diffusion + dispersion", ...
    %         false, false, true, true, false), ...
    %     makeCase("Bio-clogging", false, false, false, false, true), ...
    %     makeCase("All transport effects", true, true, true, true, false), ...
    %     makeCase("All transport effects + bio-clogging", ...
    %         true, true, true, true, true)];
    cases = [ ...
        makeCase("Baseline", false, false, false, false, false), ...
        makeCase("Microbial diffusion + chemotaxis", ...
        true, true, false, false, false), ...
        makeCase("Molecular diffusion + dispersion", ...
        false, false, true, true, false), ...
        makeCase("Bio-clogging", false, false, false, false, true), ...
        ];
end

function item = makeCase(name, bactDiffusion, chemotaxis, ...
        molecularDiffusion, mechanicalDispersion, bioClogging)
    item = struct( ...
        'name', char(name), ...
        'options', {{ ...
        'bactDiffusion', bactDiffusion, ...
        'chemotaxisEffect', chemotaxis, ...
        'molecularDiffusion', molecularDiffusion, ...
        'molecularDispersion', mechanicalDispersion, ...
        'bioClogging', bioClogging}});
end

function results = initializeResults(nEntries)
    template = struct( ...
        'backend', "", ...
        'caseName', "", ...
        'isBaseline', false, ...
        'completed', false, ...
        'runtimeSeconds', nan, ...
        'timeDays', [], ...
        'lossPercent', [], ...
        'finalLossPercent', nan, ...
        'consumedH2Moles', nan, ...
        'injectedH2Moles', nan, ...
        'reactionTotalsMoles', [], ...
        'finalSpatialConsumptionMoles', [], ...
        'xDimensionless', [], ...
        'errorMessage', "");
    results = repmat(template, nEntries, 1);
end

function metrics = computeMetrics(states, schedule, model)
    nReactions = model.biochemFluid.nbioreact;
    cumulative = zeros(numel(states), nReactions);
    finalSpatial = zeros(model.G.cells.num, nReactions);
    for reaction = 1:nReactions
        [~, reactionCumulative] = computeH2Consumption( ...
            states, schedule, model, reaction);
        cumulative(:, reaction) = sum(reactionCumulative, 1).';
        finalSpatial(:, reaction) = reactionCumulative(:, end);
    end
    total = sum(cumulative, 2);
    injected = prescribedInjectedH2(schedule, model);
    x = model.G.cells.centroids(:, 1);
    x = (x - min(x))./max(max(x) - min(x), eps);
    metrics = struct( ...
        'timeDays', cumsum(schedule.step.val(:))./day, ...
        'lossPercent', 100.*total./injected, ...
        'finalLossPercent', 100.*total(end)./injected, ...
        'consumedH2Moles', total(end), ...
        'injectedH2Moles', injected, ...
        'reactionTotalsMoles', cumulative(end, :), ...
        'finalSpatialConsumptionMoles', sum(finalSpatial, 2), ...
        'xDimensionless', x);
end

function injected = prescribedInjectedH2(schedule, model)
    names = model.EOSModel.CompositionalMixture.names;
    idxH2 = find(strcmpi(names, 'H2'), 1);
    assert(~isempty(idxH2), 'The compositional mixture does not contain H2.');
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
        error('runTransportEffectsSensitivity:MissingPhreeqcDatabase', ...
            ['PHREEQC_Modified.DAT was not found. Pass its absolute path ', ...
            'using ''phreeqcDatabaseFile'', databaseFile, or request only ', ...
            'the original backend with ''backends'', {''original''}.']);
    end
    [~, name, extension] = fileparts(databaseFile);
    assert(strcmpi([name, extension], 'PHREEQC_Modified.DAT'), ...
        'The database must be PHREEQC_Modified.DAT.');
end

function base = findBaseline(results, backendName)
    base = NaN;
    for i = 1:numel(results)
        if results(i).completed && results(i).isBaseline && ...
                results(i).backend == backendName
            base = results(i).finalLossPercent;
            return;
        end
    end
end

function summary = buildSummary(results, opt)
    n = numel(results);
    backend = strings(n, 1);
    caseName = strings(n, 1);
    status = strings(n, 1);
    finalLoss = nan(n, 1);
    delta = nan(n, 1);
    relative = nan(n, 1);
    consumed = nan(n, 1);
    runtime = nan(n, 1);
    direction = strings(n, 1);
    meaningful = strings(n, 1);
    for i = 1:n
        backend(i) = results(i).backend;
        caseName(i) = results(i).caseName;
        runtime(i) = results(i).runtimeSeconds;
        if ~results(i).completed
            status(i) = "FAILED";
            direction(i) = "n/a";
            meaningful(i) = "n/a";
            continue;
        end
        status(i) = "ok";
        finalLoss(i) = results(i).finalLossPercent;
        consumed(i) = results(i).consumedH2Moles;
        base = findBaseline(results, results(i).backend);
        if ~isnan(base)
            delta(i) = finalLoss(i) - base;
            relative(i) = 100.*delta(i)./max(abs(base), eps);
        end
        if isnan(delta(i))
            direction(i) = "n/a";
            meaningful(i) = "n/a";
        elseif abs(delta(i)) < opt.meaningfulLossThreshold
            direction(i) = "No material change";
            meaningful(i) = "No";
        elseif delta(i) > 0
            direction(i) = "Increased loss";
            meaningful(i) = "Yes";
        else
            direction(i) = "Decreased loss";
            meaningful(i) = "Yes";
        end
    end
    summary = table(backend, caseName, status, finalLoss, delta, relative, ...
        consumed, runtime, direction, meaningful, ...
        'VariableNames', {'Backend', 'Case', 'Status', 'FinalH2Loss_percent', ...
        'DeltaFromBaseline_percentagePoints', 'RelativeChange_percent', ...
        'ConsumedH2_mol', 'Runtime_seconds', 'Direction', 'Meaningful'});
end

function createPaperFigures(results, backends, cases, opt)
    nCases = numel(cases);
    caseNames = string({cases.name});
    for b = 1:numel(backends)
        sub = results((b - 1)*nCases + (1:nCases));
        done = [sub.completed];
        if ~any(done)
            continue;
        end
        label = string(backends(b).label);
        slug = matlab.lang.makeValidName(lower(char(label)));
        colors = paperColors(nCases);

        fig = paperFigure([22, 12], char(label + " -- H2 loss over time"));
        ax = axes(fig);
        hold(ax, 'on');
        for i = find(done)
            plot(ax, sub(i).timeDays, sub(i).lossPercent, ...
                'LineWidth', 1.6, 'Color', colors(i, :), ...
                'DisplayName', cases(i).name);
        end
        xlabel(ax, 'Time (days)');
        ylabel(ax, 'Consumed injected H_2 (%)');
        title(ax, label);
        legend(ax, 'Location', 'eastoutside');
        styleAxes(ax);
        paperExport(fig, sprintf('sensitivity_%s_H2_loss_over_time', slug));

        base = findBaseline(results, label);
        vals = nan(nCases, 1);
        for i = find(done)
            vals(i) = sub(i).finalLossPercent - base;
        end
        keep = ~isnan(vals);
        fig = paperFigure([20, 13], char(label + " -- delta H2 loss"));
        ax = axes(fig);
        names = categorical(caseNames(keep), caseNames(keep));
        v = vals(keep);
        bars = barh(ax, names, v, 'FaceColor', 'flat');
        bars.CData = repmat([0.16, 0.48, 0.70], numel(v), 1);
        bars.CData(v < 0, :) = repmat([0.78, 0.30, 0.22], nnz(v < 0), 1);
        xline(ax, 0, 'k-', 'LineWidth', 0.8);
        xline(ax, opt.meaningfulLossThreshold, 'k:', 'LineWidth', 1);
        xline(ax, -opt.meaningfulLossThreshold, 'k:', 'LineWidth', 1);
        xlabel(ax, '\Delta final H_2 loss (percentage points)');
        title(ax, label);
        styleAxes(ax);
        paperExport(fig, sprintf('sensitivity_%s_delta_H2_loss', slug));

        selected = [1, 2, 3, 4, 5, 6, 9, 10];
        selected = selected(ismember(selected, find(done)));
        fig = paperFigure([22, 12], char(label + " -- spatial H2 consumption"));
        ax = axes(fig);
        hold(ax, 'on');
        for i = selected
            plot(ax, sub(i).xDimensionless, ...
                sub(i).finalSpatialConsumptionMoles, ...
                'LineWidth', 1.6, 'Color', colors(i, :), ...
                'DisplayName', cases(i).name);
        end
        xlabel(ax, 'Dimensionless distance from injector');
        ylabel(ax, 'Cumulative H_2 consumption (mol cell^{-1})');
        title(ax, label);
        legend(ax, 'Location', 'eastoutside');
        styleAxes(ax);
        paperExport(fig, sprintf('sensitivity_%s_spatial_H2_consumption', slug));
    end
end


%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%}
