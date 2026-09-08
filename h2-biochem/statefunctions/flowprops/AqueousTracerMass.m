classdef AqueousTracerMass < StateFunction
    % Moles of mobile non-volatile aqueous tracers per grid cell.
    %
    % SYNOPSIS:
    %   tm = AqueousTracerMass(model)
    %
    % DESCRIPTION:
    %   Computes the amount of each active aqueous tracer in each grid
    %   cell, accounting for pore volume and liquid saturation. These
    %   tracers are non-volatile: unlike EOS components, they never enter
    %   the flash and are transported only by advection with the liquid
    %   phase (see BiochemicalFlowDiscretization.tracerConservationEquation).
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with mobile aqueous tracers enabled
    %
    % SEE ALSO:
    %   BacterialMass, BiochemistryModel, SRBTracerConvRate

    methods
        function tm = AqueousTracerMass(model, varargin)
            tm@StateFunction(model, varargin{:});
            if ismethod(model, 'getAqueousTracerNames')
                tracerNames = model.getAqueousTracerNames();
            else
                tracerNames = {'SO4', 'HS'};
            end
            tm = tm.dependsOn(cellfun(@lower, tracerNames, 'UniformOutput', false), 'state');
            tm = tm.dependsOn('s', 'state');
            tm = tm.dependsOn('PoreVolume', 'PVTPropertyFunctions');
            tm.label = 'M_{tracer}';
        end

        function m = evaluateOnDomain(prop, model, state)
            s   = model.getProp(state, 's');
            pv  = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
            L_ix = model.getLiquidIndex();
            if iscell(s)
                sL = max(s{L_ix}, 1.0e-8);
            else
                sL = max(s(:, L_ix), 1.0e-8);
            end

            if ismethod(model, 'getAqueousTracerNames')
                tracerNames = model.getAqueousTracerNames();
            else
                tracerNames = {'SO4', 'HS'};
            end
            m = cell(1, numel(tracerNames));
            for i = 1:numel(tracerNames)
                concentration = model.getProp(state, lower(tracerNames{i}));
                m{i} = pv .* sL .* concentration;
            end
        end
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
