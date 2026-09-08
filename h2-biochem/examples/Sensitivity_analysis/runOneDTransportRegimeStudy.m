function study = runOneDTransportRegimeStudy(varargin)
% Regime analysis: when do transport and motility processes affect H2 loss?
%
% The UGFACT-style 1D benchmark is advection-dominated and its reactive
% front is H2-saturated, so molecular diffusion, mechanical dispersion,
% microbial diffusion and chemotaxis have almost no effect on the
% integrated biological H2 loss (see runTransportEffectsSensitivity). This
% script rebuilds the 1D case in a regime where those processes *can*
% matter -- a slower injection and a substrate-limited front -- and then
% isolates their effect. The changes are physically motivated (lower
% Peclet, higher Damkohler-to-Peclet ratio) rather than an arbitrary
% inflation of the transport coefficients.
%
% Cases (all share the sensitive-regime base state):
%   1. Baseline                          -- no transport/motility effects
%   2. Chemotaxis + microbial diffusion  -- bactDiffusion + chemotaxisEffect
%   3. Molecular diffusion + dispersion  -- molecularDiffusion + molecularDispersion
%   4. Bio-clogging                      -- porosity/permeability feedback
%   5. All effects                       -- 2 + 3 + 4 together
%
% All reactions are MRST-owned (BiochemistryModel, reduced fixed-pH
% carbonate buffer, no PHREEQC), so every transport switch is available and
% the run is fast.
%
% OUTPUTS:
%   study - struct array, one entry per case, with the loss history,
%           spatial profile and runtime. A results table, a dimensionless
%           regime table, and four figures (exported to EPS + PNG in the
%           'figures' sub-folder) are produced.
%
% OPTIONAL PARAMETERS (sensitive-regime knobs, all physically motivated):
%   'rateScale'        - Multiplier on the injection/production rate.
%                        Default 0.5 (halves the Peclet number).
%   'h2HalfSatScale'   - Multiplier on the H2 Monod half-saturation
%                        constant. Default 25 (moves the front from
%                        H2-saturated to H2-limited).
%   'growthScale'      - Multiplier on the maximum specific growth rate.
%                        Default 0.7 (slows the near-well biomass bloom so
%                        the front stays transport-limited).
%   'microbeDiffScale' - Multiplier on the microbial diffusion coefficient.
%                        Default 60.
%   'chemotaxisScale'  - Multiplier on the chemotaxis coefficient.
%                        Default 60.
%   'dispersivityScale'- Multiplier on the mechanical-dispersion
%                        dispersivities. Default 3.
%   'cloggingStrengthScale' - Strength of the bio-clogging feedback,
%                        relative to the module default. 0 removes it, 1 is
%                        the default. Sweep it (0, 0.25, 0.5, 1, 2) to
%                        verify the sign and magnitude of the bio-clogging
%                        effect are smooth rather than a model artifact.
%   'cloggingMode'     - Which part of the bio-clogging feedback to keep:
%                        'full'        (default) pore-volume + permeability,
%                        'porevolume'  only the pore-volume reduction
%                                      (concentrates H2/biomass -> raises
%                                      the rate in the substrate-limited
%                                      regime),
%                        'permeability' only the Kozeny-Carman permeability
%                                      reduction (alters the flow field).
%                        The 'full' effect is the sum of two opposing
%                        mechanisms and can be non-monotone; run the two
%                        parts separately to get a clean, reportable sign
%                        for each.
%   'meaningfulLossThreshold' - |delta| (percentage points) below which a
%                        case is reported as "no material change".
%                        Default 0.1.
%   'cloggingReference' - 'initial' (default) normalises the clogging
%                        feedback so the run starts from the same porosity
%                        and permeability as the no-clogging baseline
%                        (clogging then acts only on biomass grown past its
%                        initial value) -- the correct choice for a
%                        clogging-vs-baseline comparison. 'zero' is the
%                        legacy module behaviour, where the rock is already
%                        partly clogged at the initial biomass.
%   'gridCells'        - Number of cells in the 1D column. Default 50.
%                        Refine (100, 200) to check whether a non-monotone
%                        bio-clogging response is numerical or physical.
%   'cases'            - Which process cases to run, as a string array or
%                        cell array of names, or 'all' (default). The
%                        Baseline case is always included as the reference.
%                        Use e.g. 'cases', "Bio-clogging" to sweep only
%                        clogging without paying for the other cases.
%
% SEE ALSO:
%   runTransportEffectsSensitivity, setupH2StorageExampleWithSRB_benchmark,
%   paperFigure, paperExport

    mrstModule add ad-props compositional deckformat h2-biochem

    % Pull 'cases' out of varargin before merge_options, which type-checks
    % against the default and would reject a string/cellstr value.
    [wantCases, varargin] = extractOption(varargin, 'cases', 'all');

    opt = struct( ...
        'rateScale',               0.5, ...
        'h2HalfSatScale',          25, ...
        'growthScale',             0.7, ...
        'microbeDiffScale',        60, ...
        'chemotaxisScale',         60, ...
        'dispersivityScale',       3, ...
        'cloggingStrengthScale',   1, ...
        'cloggingMode',            'full', ...
        'cloggingReference',       'initial', ...
        'gridCells',               50, ...
        'meaningfulLossThreshold', 0.1);
    opt = merge_options(opt, varargin{:});
    opt.cases = wantCases;
    opt.cloggingReference = validatestring(opt.cloggingReference, ...
        {'initial', 'zero'}, mfilename, 'cloggingReference');
    checkPositive(opt, {'rateScale', 'h2HalfSatScale', 'growthScale', ...
        'microbeDiffScale', 'chemotaxisScale', 'dispersivityScale'});
    validateattributes(opt.gridCells, {'numeric'}, ...
        {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'gridCells');
    validateattributes(opt.cloggingStrengthScale, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, ...
        'cloggingStrengthScale');
    opt.cloggingMode = validatestring(opt.cloggingMode, ...
        {'full', 'porevolume', 'permeability'}, mfilename, 'cloggingMode');

    commonOptions = { ...
        'rate', 'medrate', ...
        'scheduleMode', 'complete', ...
        'gridCells', opt.gridCells, ...
        'injectionCO2', 0, ...
        'carbonateBuffer', true, ...
        'carbonateBufferPH', 5.9, ...
        'initialHCO3', 4.4119e-2, ...
        'equilibrateInitialCO2', true, ...
        'phreeqcTimestepCoupling', false, ...
        'bacteriamodel', true};

    cases = selectCases(transportCases(), opt.cases);
    study = initStudy(numel(cases));

    for c = 1:numel(cases)
        fprintf('\nCase %d/%d: %s\n', c, numel(cases), cases(c).name);
        study(c).name       = string(cases(c).name);
        study(c).isBaseline = (c == 1);
        setupOptions = [commonOptions, cases(c).options];
        [~, model, schedule, state0] = ...
            setupH2StorageExampleWithSRB_benchmark(setupOptions{:});
        [model, schedule, regime] = applySensitiveRegime( ...
            model, schedule, state0, opt, cases(c).wantClog);
        study(c).regime = regime;

        solver = NonLinearSolver();
        solver.maxTimestepCuts = 16;
        timer = tic();
        try
            [ws, states, report] = simulateScheduleAD( ...
                state0, model, schedule, 'nonlinearSolver', solver);
            if isstruct(report) && isfield(report, 'Failure')
                assert(~report.Failure, 'Simulation report indicates failure.');
            end
            m = computeMetrics(states, schedule, model);
            d = cloggingDiagnostics(ws, states, model);
            study(c).completed   = true;
            study(c).timeDays    = m.timeDays;
            study(c).lossPercent = m.lossPercent;
            study(c).finalLoss   = m.finalLoss;
            study(c).consumedMol = m.consumedMol;
            study(c).injectedMol = m.injectedMol;
            study(c).xDless      = m.xDless;
            study(c).finalSpatialMol = m.finalSpatialMol;
            study(c).minPoroRatio  = d.minPoroRatio;
            study(c).maxInjBHP_bar = d.maxInjBHP_bar;
        catch err
            warning('runOneDTransportRegimeStudy:caseFailed', ...
                'Case "%s" did not converge: %s', cases(c).name, err.message);
            study(c).completed = false;
        end
        study(c).runtimeSec = toc(timer);
    end

    summary   = buildSummary(study, opt);
    regimeTab = buildRegimeTable(study(1).regime);
    disp(summary);
    disp(regimeTab);
    makeFigures(study, regimeTab, opt);
end

% ------------------------------------------------------------------------
function cases = transportCases()
    cases = [ ...
        makeCase("Baseline", ...
            false, false, false, false, false), ...
        makeCase("Chemotaxis + microbial diffusion", ...
            true,  true,  false, false, false), ...
        makeCase("Molecular diffusion + dispersion", ...
            false, false, true,  true,  false), ...
        makeCase("Bio-clogging", ...
            false, false, false, false, true), ...
        makeCase("All effects", ...
            true,  true,  true,  true,  true)];
end

function [val, rest] = extractOption(args, name, default)
    val  = default;
    rest = args;
    for i = numel(args)-1:-1:1
        if (ischar(args{i}) || (isstring(args{i}) && isscalar(args{i}))) ...
                && strcmpi(char(args{i}), name)
            val  = args{i+1};
            rest = args([1:i-1, i+2:end]);
            return;
        end
    end
end

function cases = selectCases(cases, want)
    % Keep the Baseline (reference) plus any requested cases.
    if (ischar(want) || (isstring(want) && isscalar(want))) && ...
            strcmpi(want, 'all')
        return;
    end
    want = string(want);
    names = string({cases.name});
    keep  = strcmpi(names, "Baseline") | ismember(lower(names), lower(want));
    assert(any(keep & ~strcmpi(names, "Baseline")), ...
        'None of the requested case names matched: %s', strjoin(want, ', '));
    cases = cases(keep);
end

function item = makeCase(name, bactDiff, chemo, molDiff, molDisp, clog)
    % Bio-clogging is applied by applySensitiveRegime (so its strength can
    % be swept via 'cloggingStrengthScale'), never by the setup function.
    item = struct('name', char(name), 'wantClog', logical(clog), ...
        'options', {{ ...
        'bactDiffusion',       bactDiff, ...
        'chemotaxisEffect',    chemo, ...
        'molecularDiffusion',  molDiff, ...
        'molecularDispersion', molDisp, ...
        'bioClogging',         false}});
end

% ------------------------------------------------------------------------
function [model, schedule, regime] = applySensitiveRegime( ...
        model, schedule, state0, opt, wantClog)
    % 1. Slow the injection/production (lowers the Peclet number). Keep the
    %    injector rate-controlled even under heavy clogging by lifting its
    %    BHP limit, so the clogging effect on H2 loss is not confounded by
    %    the well silently switching to pressure control.
    for k = 1:numel(schedule.control)
        for w = 1:numel(schedule.control(k).W)
            W = schedule.control(k).W(w);
            W.val = W.val * opt.rateScale;
            if W.sign > 0
                if ~isfield(W, 'lims') || isempty(W.lims)
                    W.lims = struct();
                end
                W.lims.bhp = 1000 * barsa;
            end
            schedule.control(k).W(w) = W;
        end
    end

    % 2. Substrate-limited front + slower bloom + more mobile microbes.
    bf = model.biochemFluid;
    bf.alphaH2      = bf.alphaH2      .* opt.h2HalfSatScale;
    bf.alphasub     = bf.alphasub     .* opt.h2HalfSatScale;
    bf.Psigrowthmax = bf.Psigrowthmax .* opt.growthScale;
    bf.bbact        = bf.bbact        .* opt.growthScale;
    bf.bactdiff     = bf.bactdiff     .* opt.microbeDiffScale;
    bf.xch_seuil    = bf.xch_seuil    .* opt.chemotaxisScale;
    model.biochemFluid = bf;   % kinetic/motility params are read live

    % 3. Bio-clogging (rebuilds operators / state functions, so it must
    %    precede scaleDispersivity). Its strength is set by scaling the
    %    characteristic biomass nc: larger nc -> weaker clogging, so
    %    cloggingStrengthScale -> 0 removes the feedback and
    %    cloggingStrengthScale = 1 reproduces the module default.
    if wantClog
        ncBase = [180, 180, 180];
        nc     = ncBase ./ sqrt(opt.cloggingStrengthScale);
        cp     = [0.5, 0.5, 0.5];
        nbact0 = value(state0.nbact(1, :));
        [model, poro0, perm0] = setupBioCloggingModel(model, nbact0, nc, cp, ...
            true, opt.cloggingReference);

        % Isolate one mechanism by freezing the other, then rebuild the
        % operators/state functions that depend on the rock handles.
        switch opt.cloggingMode
            case 'porevolume'
                model.rock.perm = perm0;                 % flow field unchanged
            case 'permeability'
                model.rock.poro = poro0;                 % pore volume unchanged
                model.fluid.pvMultR = @(varargin) 1;
        end
        if ~strcmp(opt.cloggingMode, 'full')
            model = model.setupOperators();
            model.FlowDiscretization = BiochemicalFlowDiscretization(model);
            model = model.setupStateFunctionGroupings();
        end
    end

    % 4. Larger dispersivity (still O(cell size)).
    model = scaleDispersivity(model, opt.dispersivityScale);

    regime = regimeNumbers(model, schedule, opt);
end

function d = cloggingDiagnostics(ws, states, model)
    % Sanity checks for the bio-clogging cases: how far did porosity fall
    % and did the injector stay well-behaved.
    d = struct('minPoroRatio', 1, 'maxInjBHP_bar', nan);
    try
        phi0 = mean(value(baseValue(model.rock.poro)));
        r = 1;
        for i = 1:numel(states)
            pv = model.getProp(states{i}, 'PoreVolume');
            r  = min(r, min(value(pv) ./ model.G.cells.volumes) / phi0);
        end
        d.minPoroRatio = r;
    catch
        % leave default
    end
    try
        % Only look at the injector: during the shut-in and production
        % controls well 1 is the shut-in/producer, whose BHP constraint
        % would otherwise leak into the maximum.
        bhpMax = nan;
        for i = 1:numel(ws)
            if isempty(ws{i})
                continue;
            end
            for w = 1:numel(ws{i})
                isInj = (isfield(ws{i}(w), 'sign') && ws{i}(w).sign > 0) || ...
                    (isfield(ws{i}(w), 'name') && ...
                     strcmpi(ws{i}(w).name, 'Injector'));
                if isInj && isfinite(ws{i}(w).bhp)
                    bhpMax = max([bhpMax, ws{i}(w).bhp]);
                end
            end
        end
        d.maxInjBHP_bar = bhpMax / barsa;
    catch
        % leave default
    end
end

function model = scaleDispersivity(model, s)
    % DispersiveDiffusivity only exists when molecularDiffusion or
    % molecularDispersion is enabled; skip quietly otherwise.
    if s == 1
        return;
    end
    fd = model.FlowDiscretization;
    if ~fd.hasStateFunction('DispersiveDiffusivity')
        return;
    end
    sf = fd.getStateFunction('DispersiveDiffusivity');
    for f = {'alphaL_water', 'alphaT_water', 'alphaL_gas', ...
             'alphaT_gas', 'defaultLiquidDiffusivity'}
        sf.(f{1}) = sf.(f{1}) * s;
    end
    model.FlowDiscretization = ...
        fd.setStateFunction('DispersiveDiffusivity', sf);
end

% ------------------------------------------------------------------------
function regime = regimeNumbers(model, schedule, opt)
    G   = model.G;
    x   = G.cells.centroids(:, 1);
    L   = max(x) - min(x) + (x(2) - x(1));            % column length [m]
    dx  = L / G.cells.num;
    A   = G.cells.volumes(1) / dx;                    % cross-section [m^2]
    phi = mean(baseValue(model.rock.poro));

    % Peak reservoir injection rate over the schedule [m^3/s].
    q = 0;
    for k = 1:numel(schedule.control)
        for w = 1:numel(schedule.control(k).W)
            W = schedule.control(k).W(w);
            if W.sign > 0
                q = max(q, abs(W.val));
            end
        end
    end
    v = q / max(A * phi, eps);                        % pore velocity [m/s]

    Dmol  = 4.5e-9 * phi^(4/3);                       % effective H2 diffusivity
    alphaL = 5.0e-2 * opt.dispersivityScale;          % scaled longitudinal dispersivity [m]
    Ddisp = alphaL * v;
    Db    = max(model.biochemFluid.bactdiff);
    chi   = max(model.biochemFluid.xch_seuil);
    kmax  = max(model.biochemFluid.Psigrowthmax);
    tsim  = sum(schedule.step.val);

    regime = struct( ...
        'poreVelocity_m_s', v, ...
        'columnLength_m',   L, ...
        'Pe_molecular',     v * L / max(Dmol, eps), ...
        'Pe_dispersive',    v * L / max(Ddisp + Dmol, eps), ...
        'Pe_chemotaxis',    v * L / max(chi, eps), ...
        'Pe_microbial',     v * L / max(Db, eps), ...
        'Damkohler',        kmax * L / max(v, eps), ...
        'microbialSpread_L', sqrt(max(Db, 0) * tsim) / L, ...
        'h2HalfSatScale',   opt.h2HalfSatScale, ...
        'rateScale',        opt.rateScale);
end

function v = baseValue(x)
    if isa(x, 'function_handle')
        v = x(1);          % clogging porosity handle -> evaluate at unit args
    else
        v = x;
    end
    v = value(v);
end

% ------------------------------------------------------------------------
function m = computeMetrics(states, schedule, model)
    nReact = model.biochemFluid.nbioreact;
    cum    = zeros(numel(states), nReact);
    final  = zeros(model.G.cells.num, nReact);
    for r = 1:nReact
        [~, perCell] = computeH2Consumption(states, schedule, model, r);
        cum(:, r)   = sum(perCell, 1).';
        final(:, r) = perCell(:, end);
    end
    total    = sum(cum, 2);
    injected = prescribedInjectedH2(schedule, model);
    x = model.G.cells.centroids(:, 1);
    x = (x - min(x)) ./ max(max(x) - min(x), eps);

    m = struct( ...
        'timeDays',        cumsum(schedule.step.val(:)) / day, ...
        'lossPercent',     100 .* total ./ injected, ...
        'finalLoss',       100 .* total(end) ./ injected, ...
        'consumedMol',     total(end), ...
        'injectedMol',     injected, ...
        'xDless',          x, ...
        'finalSpatialMol', sum(final, 2));
end

function injected = prescribedInjectedH2(schedule, model)
    names   = model.EOSModel.CompositionalMixture.names;
    idxH2   = find(strcmpi(names, 'H2'), 1);
    assert(~isempty(idxH2), 'The compositional mixture does not contain H2.');
    gasIdx  = model.getVaporIndex();
    injected = 0;
    for s = 1:numel(schedule.step.val)
        W = schedule.control(schedule.step.control(s)).W;
        for w = 1:numel(W)
            if W(w).sign <= 0 || (isfield(W, 'status') && ~W(w).status)
                continue;
            end
            if strcmpi(W(w).type, 'grat')
                gasRate = W(w).val;
            elseif strcmpi(W(w).type, 'rate')
                gasRate = W(w).val * W(w).compi(gasIdx);
            else
                continue;
            end
            molarRate = gasRate * model.FacilityModel.pressure / ...
                (8.314462618 * model.FacilityModel.T);
            injected = injected + molarRate * W(w).components(idxH2) * ...
                schedule.step.val(s);
        end
    end
    assert(injected > 0, 'No prescribed H2 injection was found.');
end

% ------------------------------------------------------------------------
function study = initStudy(n)
    template = struct('name', "", 'isBaseline', false, 'completed', false, ...
        'timeDays', [], 'lossPercent', [], 'finalLoss', nan, ...
        'consumedMol', nan, 'injectedMol', nan, 'xDless', [], ...
        'finalSpatialMol', [], 'runtimeSec', nan, 'regime', struct(), ...
        'minPoroRatio', nan, 'maxInjBHP_bar', nan);
    study = repmat(template, n, 1);
end

function summary = buildSummary(study, opt)
    n     = numel(study);
    name  = strings(n, 1);
    final = nan(n, 1);
    delta = nan(n, 1);
    rel   = nan(n, 1);
    cons  = nan(n, 1);
    rt    = nan(n, 1);
    poro  = nan(n, 1);
    bhp   = nan(n, 1);
    verd  = strings(n, 1);
    if study(1).completed
        base = study(1).finalLoss;
    else
        base = NaN;
    end
    for i = 1:n
        name(i) = study(i).name;
        rt(i)   = study(i).runtimeSec;
        if ~study(i).completed
            verd(i) = "FAILED (no convergence)";
            continue;
        end
        final(i) = study(i).finalLoss;
        cons(i)  = study(i).consumedMol;
        poro(i)  = study(i).minPoroRatio;
        bhp(i)   = study(i).maxInjBHP_bar;
        if isnan(base)
            verd(i) = "ok (no baseline)";
            continue;
        end
        delta(i) = study(i).finalLoss - base;
        rel(i)   = 100 * delta(i) / max(abs(base), eps);
        if abs(delta(i)) < opt.meaningfulLossThreshold
            verd(i) = "no material change";
        elseif delta(i) > 0
            verd(i) = "increased loss";
        else
            verd(i) = "decreased loss";
        end
    end
    summary = table(name, final, delta, rel, cons, rt, poro, bhp, verd, ...
        'VariableNames', {'Case', 'FinalH2Loss_percent', ...
        'DeltaFromBaseline_pp', 'RelativeChange_percent', ...
        'ConsumedH2_mol', 'Runtime_s', 'MinPoroRatio', 'MaxInjBHP_bar', ...
        'Verdict'});
end

function regimeTab = buildRegimeTable(r)
    quantity = ["Pore velocity (m/s)"; "Column length (m)"; ...
        "Peclet, molecular"; "Peclet, dispersive"; ...
        "Peclet, chemotaxis"; "Peclet, microbial diffusion"; ...
        "Damkohler (k L / v)"; "Microbial spread sqrt(D t)/L"; ...
        "H2 half-saturation scale"; "Injection-rate scale"];
    value_ = [r.poreVelocity_m_s; r.columnLength_m; ...
        r.Pe_molecular; r.Pe_dispersive; r.Pe_chemotaxis; r.Pe_microbial; ...
        r.Damkohler; r.microbialSpread_L; r.h2HalfSatScale; r.rateScale];
    regimeTab = table(quantity, value_, ...
        'VariableNames', {'Quantity', 'Value'});
end

% ------------------------------------------------------------------------
function makeFigures(study, regimeTab, opt)
    n      = numel(study);
    colors = distinctColors(n);
    names  = arrayfun(@(s) char(s.name), study, 'UniformOutput', false);
    done   = [study.completed];

    % (1) H2 loss over time
    fig = paperFigure([20, 12], 'Regime study -- H2 loss over time');
    ax = axes(fig); hold(ax, 'on');
    for i = find(done)
        plot(ax, study(i).timeDays, study(i).lossPercent, ...
            'LineWidth', 1.8, 'Color', colors(i, :), 'DisplayName', names{i});
    end
    xline(ax, 50,  'k:', 'end injection', 'HandleVisibility', 'off');
    xline(ax, 200, 'k:', 'end storage',   'HandleVisibility', 'off');
    xlabel(ax, 'Time (days)');
    ylabel(ax, 'Cumulative biological H_2 loss (%)');
    title(ax, 'Sensitive-regime 1D case');
    legend(ax, 'Location', 'northwest');
    styleAxes(ax);
    paperExport(fig, 'regime_1D_H2_loss_over_time');

    % (2) Final spatial consumption
    fig = paperFigure([20, 12], 'Regime study -- spatial H2 consumption');
    ax = axes(fig); hold(ax, 'on');
    for i = find(done)
        plot(ax, study(i).xDless, study(i).finalSpatialMol, ...
            'LineWidth', 1.8, 'Color', colors(i, :), 'DisplayName', names{i});
    end
    xlabel(ax, 'Dimensionless distance from injector, x/L');
    ylabel(ax, 'Cumulative H_2 consumed (mol cell^{-1})');
    title(ax, 'Final spatial H_2 consumption');
    legend(ax, 'Location', 'northeast');
    styleAxes(ax);
    paperExport(fig, 'regime_1D_spatial_H2_consumption');

    % (3) Delta from baseline
    if done(1)
        fig = paperFigure([20, 11], 'Regime study -- delta H2 loss');
        ax = axes(fig);
        base   = study(1).finalLoss;
        idx    = find(done); idx = idx(idx ~= 1);
        dd     = arrayfun(@(k) study(k).finalLoss - base, idx(:));
        lbls   = names(idx);
        cats   = categorical(lbls, lbls);
        b = barh(ax, cats, dd, 'FaceColor', 'flat');
        b.CData = repmat([0.16, 0.48, 0.70], numel(dd), 1);
        b.CData(dd < 0, :) = repmat([0.78, 0.30, 0.22], nnz(dd < 0), 1);
        xline(ax, 0, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        xline(ax,  opt.meaningfulLossThreshold, 'k:', 'HandleVisibility', 'off');
        xline(ax, -opt.meaningfulLossThreshold, 'k:', 'HandleVisibility', 'off');
        xlabel(ax, '\Delta final H_2 loss vs baseline (percentage points)');
        title(ax, 'Effect of each process in the sensitive regime');
        styleAxes(ax);
        paperExport(fig, 'regime_1D_delta_H2_loss');
    else
        warning('runOneDTransportRegimeStudy:noBaseline', ...
            'Baseline case failed; skipping the delta-from-baseline figure.');
    end

    % (4) Dimensionless regime numbers
    fig = paperFigure([18, 11], 'Regime study -- dimensionless numbers');
    ax = axes(fig);
    pick = ["Peclet, molecular", "Peclet, dispersive", ...
            "Peclet, chemotaxis", "Peclet, microbial diffusion", ...
            "Damkohler (k L / v)"];
    lbl  = ["Pe_{mol}", "Pe_{disp}", "Pe_{chemo}", "Pe_{micro}", "Da"];
    vals = zeros(numel(pick), 1);
    for i = 1:numel(pick)
        vals(i) = regimeTab.Value(regimeTab.Quantity == pick(i));
    end
    barh(ax, categorical(lbl, lbl), vals);
    set(ax, 'XScale', 'log');
    xline(ax, 1, 'k--', 'transport = advection', 'HandleVisibility', 'off');
    xlabel(ax, 'Dimensionless group (log scale)');
    title(ax, 'Where the sensitive-regime case sits');
    styleAxes(ax);
    paperExport(fig, 'regime_1D_dimensionless_numbers');
end

function c = distinctColors(n)
    base = [ ...
        0.15 0.15 0.15;   % baseline   - near-black
        0.00 0.45 0.74;   % blue
        0.85 0.33 0.10;   % orange-red
        0.47 0.67 0.19;   % green
        0.49 0.18 0.56];  % purple
    if n <= size(base, 1)
        c = base(1:n, :);
    else
        c = [base; lines(n - size(base, 1))];
    end
end

% ------------------------------------------------------------------------
function checkPositive(opt, fields)
    for i = 1:numel(fields)
        validateattributes(opt.(fields{i}), {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'}, mfilename, fields{i});
    end
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
%}
