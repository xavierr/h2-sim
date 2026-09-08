function fig = paperFigure(sizeCm, name)
% Create a white, centimeter-sized figure window for paper/presentation
% figures. Shared by every 1D_validation/ and Sensitivity_analysis/
% script so all module figures share one physical sizing convention.
%
% SYNOPSIS:
%   fig = paperFigure([18, 11])
%   fig = paperFigure([24, 15], 'Aqueous H2 and CO2 profiles')
%
% PARAMETERS:
%   sizeCm - [width, height] figure size in centimeters.
%   name   - Optional figure name (default: '').
%
% SEE ALSO:
%   styleAxes, paperColors

    if nargin < 2
        name = '';
    end
    fig = figure('Name', name, 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, sizeCm], ...
        'DefaultAxesFontName', 'Times New Roman', ...
        'DefaultAxesFontSize', 12, ...
        'DefaultTextFontName', 'Times New Roman', ...
        'DefaultTextFontSize', 12, ...
        'PaperPositionMode', 'auto');
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
