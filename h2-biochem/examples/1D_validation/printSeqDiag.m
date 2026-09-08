% printSeqDiag - interactive post-processing diagnostics for the last
% exampleSequentialBiochemistryPhreeqc1D run (Workflow 3, the
% coarse-flow/local-reaction PHREEQC split). Paste-run after that example
% to check H2 loss and catch resource-starved/runaway kinetics.
%
% Requires `summary` (the struct returned by exampleSequentialBiochemistryPhreeqc1D)
% already in the workspace. Prints, per reaction: cumulative H2 consumed,
% final biomass relative to its cap, how hard the carbon/H2 resource
% limiter clipped the kinetics each step, and the aqueous carbon pool
% trend. Not a function; not called from anywhere else in the module.

bf = summary.model.biochemFluid;

%% 1) Per-reaction cumulative H2 at end of run
% Compositional reference @ day 50: MET 512.15, ACE 2.2231, SRB 41.887 mol
fprintf('\nPer-reaction H2 consumed @ end of run:\n');
for r = 1:bf.nbioreact
    [~, cum] = computeH2Consumption(summary.states, summary.schedule, summary.model, r);
    fprintf('  %-25s %10.4g mol\n', bf.metabolicReaction{r}, sum(cum(:, end)));
end

%% 2) Biomass growth: cap is N/N0 = 1e4
nb = value(summary.states{end}.nbact);
fprintf('\nFinal N/N0 per reaction (max over cells):\n');
for r = 1:bf.nbioreact
    fprintf('  %-25s %10.4g\n', bf.metabolicReaction{r}, max(nb(:, r)));
end

%% 3) How hard the resource limiter clipped the kinetics each step
reports = summary.report.ControlstepReports;
extMin = nan(numel(reports), 1);
for k = 1:numel(reports)
    try
        nr = reports{k}.StepReports{1}.NonlinearReport{1};
        st = nr.SequentialBiochemistryPhreeqcStages;
        extMin(k) = min(cellfun(@(s) s.ExtentScaleMinimum, ...
            st.ReactionSubstepReports));
    catch
    end
end
fprintf(['\nExtentScaleMinimum (1 = unlimited kinetics, ', ...
    '<<1 = resource-starved):\n']);
fprintf('  min %.3g | median %.3g | steps below 0.5: %d of %d\n', ...
    min(extMin), median(extMin, 'omitnan'), ...
    nnz(extMin < 0.5), nnz(~isnan(extMin)));

%% 4) Aqueous carbon pool over time
hco3 = cellfun(@(s) max(s.tracerHCO3), summary.states);
fprintf(['\ntracerHCO3 max over cells: start %.3g | min %.3g | ', ...
    'end %.3g\n'], hco3(1), min(hco3), hco3(end));
