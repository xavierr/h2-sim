function c = paperColors(n)
% Return n RGB color rows for paper/presentation figures with multiple
% series. The first 4 rows are the module's fixed palette (used by
% runThreeBackendComparison.m's 3-4 backend comparisons); beyond that,
% additional distinguishable colors are drawn from MATLAB's 'lines'
% colormap (needed by runTransportEffectsSensitivity.m's 10-case sweep).
% Shared so every 1D_validation/ and Sensitivity_analysis/ figure with
% more than one series uses the same colors for the same number of series.
%
% SYNOPSIS:
%   c = paperColors(4)
%   c = paperColors(10)
%
% SEE ALSO:
%   paperFigure, styleAxes

    fixed = [ ...
        0.15, 0.15, 0.15; ...
        0.00, 0.45, 0.74; ...
        0.85, 0.33, 0.10; ...
        0.47, 0.67, 0.19];
    validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive'}, ...
        mfilename, 'n');
    if n <= size(fixed, 1)
        c = fixed(1:n, :);
    else
        c = [fixed; lines(n - size(fixed, 1))];
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
