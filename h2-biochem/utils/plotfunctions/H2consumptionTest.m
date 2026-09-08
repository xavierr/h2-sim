% H2consumptionTest - interactive post-processing diagnostic: reports and
% plots per-reaction H2 consumption (MET/ACE/SRB) for a completed
% simulation, and the fraction of injected H2 consumed overall.
%
% Requires `states`, `ws`, `model`, `schedule` (the outputs of a
% simulateScheduleAD/simulateSequentialH2BiochemPhreeqc run) already in
% the workspace. Not a function; not called from anywhere else in the
% module. Paste-run after any of the module's example drivers to sanity
% check H2 loss against the expected benchmark range.

for aa = 1:1
    assert(numel(states) == numel(schedule.step.val), ...
        'states and schedule must contain the same number of timesteps.');
    assert(numel(ws) == numel(schedule.step.val), ...
        'well solutions and schedule must contain the same number of timesteps.');

    eosNames = model.EOSModel.CompositionalMixture.names;
    idxH2 = find(strcmp(eosNames, 'H2'), 1);
    assert(~isempty(idxH2), 'The compositional mixture does not contain H2.');
    % The imposed rate and injector composition define the amount supplied
    % to the model. ComponentTotalFlux is not suitable for this accounting:
    % with molecular diffusion it includes the local diffusive contribution
    % and can therefore report several times the prescribed H2 input.
    gasIndex = model.getVaporIndex();
    gasConstant = 8.314462618; % J/(mol K)
    surfacePressure = model.FacilityModel.pressure;
    surfaceTemperature = model.FacilityModel.T;
    totalInjectedH2 = 0;

    for t = 1:numel(schedule.step.val)
        control = schedule.step.control(t);
        wells = schedule.control(control).W;

        for w = 1:numel(wells)
            if wells(w).sign <= 0 || ...
                    (~strcmpi(wells(w).type, 'bhp') && wells(w).val == 0) || ...
                    (isfield(wells, 'status') && ~wells(w).status)
                continue;
            end

            if ~isfield(wells, 'components') || numel(wells(w).components) < idxH2
                error('Injection control at timestep %d has no H2 component composition.', t);
            end
            if strcmpi(wells(w).type, 'grat')
                gasRate = wells(w).val;
            elseif strcmpi(wells(w).type, 'rate')
                gasRate = wells(w).val*wells(w).compi(gasIndex);
            else
                error(['Cannot compute prescribed H2 input for injection control type "%s" ' ...
                    'at timestep %d.'], wells(w).type, t);
            end

            h2MolarRate = gasRate*surfacePressure/(gasConstant*surfaceTemperature) * ...
                wells(w).components(idxH2);
            totalInjectedH2 = totalInjectedH2 + h2MolarRate*schedule.step.val(t);
        end
    end
    fprintf('Total injected H2: %.2f mol\n', totalInjectedH2);
    totalInjectedH2_all(aa) = totalInjectedH2;

    if model.bacteriamodel

        nReactions = model.biochemFluid.nbioreact;
        H2cum = cell(nReactions, 1);
        totalCum = zeros(numel(states), nReactions);
        finalCum = zeros(model.G.cells.num, nReactions);

        for reaction = 1:nReactions
            [~, H2cum{reaction}] = computeH2Consumption( ...
                states, schedule, model, reaction);
            totalCum(:, reaction) = sum(H2cum{reaction}, 1)';
            finalCum(:, reaction) = H2cum{reaction}(:, end);
        end
        assert(all(isfinite(finalCum(:)) & finalCum(:) >= 0), ...
            'Computed H2 consumption contains invalid values.');
        reactionTotals = sum(finalCum, 1);
        totalConsumedH2 = sum(reactionTotals);

        timeDays = cumsum(schedule.step.val)./day;
        reactionNames = cellstr(model.biochemFluid.metabolicReaction);
        reactionNames = reactionNames(:);

        figure;
        plot(timeDays, totalCum, 'LineWidth', 2);
        hold on;
        plot(timeDays, sum(totalCum, 2), 'k--', 'LineWidth', 2);
        xlabel('Time (days)');
        ylabel('Cumulative H2 consumed (mol)');
        legend([reactionNames; {'Total'}], 'Location', 'best');
        title('Total H2 Consumption over Time');
        grid on;

        xCoords = model.G.cells.centroids(:, 1);
        xNorm = (xCoords - min(xCoords))./(max(xCoords) - min(xCoords));

        figure;
        plot(xNorm, finalCum, 'LineWidth', 2);
        hold on;
        plot(xNorm, sum(finalCum, 2), 'k--', 'LineWidth', 2);
        xlabel('Dimensionless length (distance from injector)');
        ylabel('Cumulative H2 consumed per cell (mol)');
        legend([reactionNames; {'Total'}], 'Location', 'best');
        title('Spatial Distribution of H2 Consumption');
        grid on;

        for reaction = 1:nReactions
            fprintf('%s consumed H2: %.6f mol\n', ...
                reactionNames{reaction}, reactionTotals(reaction));
        end
        fprintf('Total consumed H2: %.6f mol\n', totalConsumedH2);
        fprintf('Consumed injected H2: %.3f %%\n', ...
            100*totalConsumedH2./totalInjectedH2);
        % ---- Spatial distribution as stacked bar chart ----
        xCoords = model.G.cells.centroids(:, 1);
        xNorm = (xCoords - min(xCoords))./(max(xCoords) - min(xCoords));

        figure;
        % Stacked bars: each cell's total height = sum of all reactions
        hb = bar(xNorm, finalCum, 'stacked');
        xlabel('Dimensionless length (distance from injector)');
        ylabel('Cumulative H₂ consumed per cell (mol)');
        title('Spatial Distribution of H₂ Consumption');
        legend(reactionNames, 'Location', 'best');
        grid on;
        ylim([0, 250]);
    end
end