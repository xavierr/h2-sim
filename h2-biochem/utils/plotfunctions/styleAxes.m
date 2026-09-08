function styleAxes(ax, fontSize)
% Apply the module's shared paper-figure axis style: gridded, boxed,
% Times New Roman, outward ticks, 12 pt. Shared by every example script that
% produces analysis figures so the whole module has one visual style.
%
% SYNOPSIS:
%   styleAxes(ax)
%   styleAxes(ax, 14)   % override the default 12 pt font
%
% The font size is applied to the axes, its title, its x/y/z labels, and any
% legend attached to it, so a single call is enough to standardise a plot.
%
% SEE ALSO:
%   paperFigure, paperColors, paperExport

    if nargin < 1 || isempty(ax)
        ax = gca;
    end
    if nargin < 2 || isempty(fontSize)
        fontSize = 12;
    end
    grid(ax, 'on');
    box(ax, 'on');
    ax.FontName  = 'Times New Roman';
    ax.FontSize  = fontSize;
    ax.LineWidth = 0.8;
    ax.TickDir   = 'out';
    ax.Layer     = 'top';

    handles = [ax.Title, ax.XLabel, ax.YLabel, ax.ZLabel];
    set(handles, 'FontName', 'Times New Roman', 'FontSize', fontSize);

    lg = get(ax, 'Legend');
    if ~isempty(lg)
        set(lg, 'FontName', 'Times New Roman', 'FontSize', fontSize, ...
            'Box', 'off');
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
