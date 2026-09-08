function plotH2LossVsBaseline(model, schedule, statesNoBact, statesBact)
% plotH2LossVsBaseline  H2 loss attributable to bacteria, via direct diff
%
%   loss(t)        = mass_noBacteria(t) - mass_withBacteria(t)
%   lossPercent(t)  = 100 * loss(t) / mass_noBacteria(t)   (instantaneous)
%
%   Prints total lost H2 in kg and as % of initial H2 mass.
%
%   INPUTS:
%     model         - fluid model (must contain compFluid.names)
%     schedule      - shared schedule (for time axis)
%     statesNoBact  - states from the abiotic (no‑bacteria) run
%     statesBact    - states from the biotic (with‑bacteria) run

compNames = model.compFluid.names;
iH2 = find(strcmpi(compNames, 'Hydrogen'));
if isempty(iH2), iH2 = 2; end

nSteps = min(numel(statesNoBact), numel(statesBact));
if numel(statesNoBact) ~= numel(statesBact)
    warning('plotH2LossVsBaseline:StepMismatch', ...
        'statesNoBact (%d) and statesBact (%d) have different lengths - truncating to %d.', ...
        numel(statesNoBact), numel(statesBact), nSteps);
end

if iscell(schedule.step.val)
    dt = vertcat(schedule.step.val{:});
else
    dt = schedule.step.val(:);
end
dt = dt(1:nSteps);
timeDays = cumsum(dt) / day;

massNoBact = zeros(nSteps, 1);
massBact   = zeros(nSteps, 1);
for k = 1:nSteps
    massNoBact(k) = sum(statesNoBact{k}.FlowProps.ComponentTotalMass{iH2});
    massBact(k)   = sum(statesBact{k}.FlowProps.ComponentTotalMass{iH2});
end

lossMass    = massNoBact - massBact;                       % positive = consumed by bacteria
lossPercent = 100 * lossMass ./ max(massNoBact, 1e-12);   % relative to instantaneous no‑bacteria mass

% Compute total loss as % of initial H2 (constant denominator)
initialMass = massNoBact(1);
totalLossPct = 100 * lossMass(end) / max(initialMass, 1e-12);

fprintf('\n=== H2 Loss Due to Bacteria (baseline-differenced) ===\n');
fprintf('Initial H2 mass (both runs):   %.6f kg\n', initialMass);
fprintf('Final H2 mass (no bacteria):   %.6f kg\n', massNoBact(end));
fprintf('Final H2 mass (with bacteria): %.6f kg\n', massBact(end));
fprintf('H2 consumed by bacteria:       %.6f kg  (%.4f%% of instantaneous no‑bacteria mass)\n', ...
    lossMass(end), lossPercent(end));
fprintf('Total H2 lost:                 %.6f kg  (%.4f%% of initial H2 mass)\n', ...
    lossMass(end), totalLossPct);

figure('Name', 'H2 loss due to bacteria (baseline-differenced)');
subplot(2,1,1);
plot(timeDays, massNoBact, 'b', timeDays, massBact, 'r', 'LineWidth', 1.5);
legend('No bacteria', 'With bacteria', 'Location', 'best');
xlabel('Time (days)'); ylabel('H_2 mass in reservoir (kg)');
title('H_2 mass: bacteria vs. no-bacteria run'); grid on;

subplot(2,1,2);
plot(timeDays, lossPercent, 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('H_2 loss due to bacteria (%)');
title(sprintf('Bacteria-attributable H_2 loss: %.4f%%', lossPercent(end)));
grid on;
end