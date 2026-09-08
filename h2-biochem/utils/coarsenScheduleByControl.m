function schedule = coarsenScheduleByControl(schedule, maxDt)
% Partition every operating period (injection/storage/production, i.e.
% every distinct schedule.control) into equal-size flow steps no larger
% than maxDt, preserving each period's total duration exactly.
%
% Used to give the coarse-flow/local-reaction split
% (SequentialBiochemistryPhreeqcModel) a coarse global flow-step
% resolution while it independently subcycles local reactions at a finer
% scale (see reactionSubstepMaxDt).
%
% SYNOPSIS:
%   schedule = coarsenScheduleByControl(schedule, 2*day)
%
% SEE ALSO:
%   SequentialBiochemistryPhreeqcModel, convertToSequentialBiochemistryPhreeqcModel

    control = schedule.step.control(:);
    val = schedule.step.val(:);
    newVal = [];
    newControl = [];
    controls = unique(control, 'stable');
    for i = 1:numel(controls)
        mask = control == controls(i);
        total = sum(val(mask));
        nSteps = ceil(total/maxDt);
        coarseVal = repmat(total/nSteps, nSteps, 1);
        newVal = [newVal; coarseVal]; %#ok<AGROW>
        newControl = [newControl; repmat(controls(i), nSteps, 1)]; %#ok<AGROW>
    end
    schedule.step.val = newVal;
    schedule.step.control = newControl;
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
