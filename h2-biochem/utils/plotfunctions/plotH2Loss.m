function plotH2Loss(model, schedule, states, ws)
% plotH2Loss  Compute and plot cumulative H2 loss (% of total introduced)
%
%   Mass balance, computed purely from states + ws (no internal kinetics):
%       loss = (initialMass + cumNetWellH2) - currentMass
%       lossPercent = 100 * loss / (initialMass + cumInjected)
%
%   well.H2 (kg/s) is used AS-IS, signed: positive = mass added to
%   reservoir (injection), negative = mass removed (production). We do
%   NOT re-derive direction from well.sign/well.val and re-apply it via
%   abs() - well.H2's sign already encodes direction, and reconstructing
%   it separately can silently disagree with the actual flow (crossflow,
%   BHP-controlled wells, etc), corrupting the balance.

% Hydrogen component index
compNames = model.compFluid.names;
iH2 = find(strcmpi(compNames, 'Hydrogen'));
if isempty(iH2), iH2 = 2; end

nSteps = numel(states);

% Time step (handles cell arrays)
if iscell(schedule.step.val)
    dt = vertcat(schedule.step.val{:});
else
    dt = schedule.step.val(:);
end
timeDays = cumsum(dt) / day;

% Total H2 mass in reservoir at each step
currentMass = zeros(nSteps, 1);
for k = 1:nSteps
    currentMass(k) = sum(states{k}.FlowProps.ComponentTotalMass{iH2});
end
initialMass = currentMass(1);

% ----- Net H2 rate from wells, taken signed and as-is -----
netH2RateStep      = zeros(nSteps, 1);   % signed: +injection, -production
injectedH2RateStep = zeros(nSteps, 1);   % positive part only (for % denominator)

for k = 1:nSteps
    if isempty(ws{k}), continue; end
    for w = 1:numel(ws{k})
        well = ws{k}(w);
        if isempty(well) || ~isfield(well, 'H2') || isempty(well.H2)
            continue;
        end
        % Ignore shut/inactive wells - .H2 is meaningless when the well
        % isn't actually flowing (val==0 and/or status==false).
        if isfield(well, 'val') && well.val == 0
            continue;
        end
        if isfield(well, 'status') && ~well.status
            continue;
        end
        h2Rate = well.H2;   % kg/s, signed
        netH2RateStep(k) = netH2RateStep(k) + h2Rate;
        if h2Rate > 0
            injectedH2RateStep(k) = injectedH2RateStep(k) + h2Rate;
        end
    end
end

cumNetWellH2 = cumsum(netH2RateStep .* dt);
cumInjected  = cumsum(injectedH2RateStep .* dt);

% ----- Mass balance -----
expectedMass = initialMass + cumNetWellH2;
lossMass     = expectedMass - currentMass;                 % positive = lost
lossPercent  = 100 * lossMass ./ max(initialMass + cumInjected, 1e-12);

% ----- Output -----
fprintf('\n=== H2 Mass Balance Summary ===\n');
fprintf('Initial H2 mass:        %.6f kg\n', initialMass);
fprintf('Total H2 injected:      %.6f kg\n', cumInjected(end));
fprintf('Expected final mass:    %.6f kg\n', expectedMass(end));
fprintf('Measured final mass:    %.6f kg\n', currentMass(end));
fprintf('H2 lost:                %.6e kg (%.4f%%)\n', lossMass(end), lossPercent(end));

% ----- Plots -----
figure('Name', 'H2 mass balance');
subplot(2,1,1);
plot(timeDays, currentMass, 'b', timeDays, expectedMass, 'r--', 'LineWidth', 1.5);
legend('Current mass', 'Expected (init + net well H_2)', 'Location', 'best');
xlabel('Time (days)'); ylabel('H_2 mass (kg)');
title('Mass balance check'); grid on;

subplot(2,1,2);
plot(timeDays, lossPercent, 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('H_2 loss (%)');
title(sprintf('Cumulative H_2 loss: %.4f%%', lossPercent(end)));
grid on;

% Optional bacteria growth plot
if isfield(states{end}, 'nbact') && ~isempty(states{end}.nbact)
    figure('Name', 'Bacterial growth');
    totalBact = zeros(nSteps, 1);
    for k = 1:nSteps
        if isfield(states{k}, 'nbact') && ~isempty(states{k}.nbact)
            totalBact(k) = sum(states{k}.nbact(:,1));
        end
    end
    plot(timeDays, totalBact, 'g', 'LineWidth', 1.5);
    xlabel('Time (days)'); ylabel('Total bacteria');
    title('Bacterial Growth'); grid on;
end
end
