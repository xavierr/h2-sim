classdef DiffusiveTracerFlux < StateFunction
    % Diffusive/dispersive flux of mobile aqueous tracers.
    %
    % SYNOPSIS:
    %   flux = DiffusiveTracerFlux(model)
    %
    % DESCRIPTION:
    %   Computes the face-wise diffusive/dispersive flux of the active
    %   tracers in the liquid phase as
    %
    %       J_i = - T_i .* Grad(c_i)
    %
    %   where T_i is the tracer transmissibility from
    %   TracerTransmissibility and c_i is the tracer concentration
    %   [mol/m^3 liquid]. Unlike DiffusiveBactFlux, no face-averaged
    %   liquid density factor is needed here: the tracers are already
    %   carried as a volumetric molar concentration (see
    %   AqueousTracerMass), not a mass fraction.
    %
    % REQUIRED PARAMETERS:
    %   model - Model with `TracerTransmissibility` state function
    %           registered and sulfateReduction enabled.
    %
    % SEE ALSO:
    %   StateFunction, TracerTransmissibility, DiffusiveBactFlux

    methods
        function df = DiffusiveTracerFlux(model, varargin)
            df@StateFunction(model, varargin{:});
            df = merge_options(df, varargin{:});

            df = df.dependsOn('TracerTransmissibility');
            if ismethod(model, 'getAqueousTracerNames')
                tracerNames = model.getAqueousTracerNames();
            else
                tracerNames = {'SO4', 'HS'};
            end
            df = df.dependsOn(cellfun(@lower, tracerNames, 'UniformOutput', false), 'state');

            df.label = 'J_{tracer}^{diff}';
        end

        function J = evaluateOnDomain(prop, model, state)
            op = model.operators;
            T = prop.getEvaluatedDependencies(state, 'TracerTransmissibility');

            if ismethod(model, 'getAqueousTracerNames')
                tracerNames = model.getAqueousTracerNames();
            else
                tracerNames = {'SO4', 'HS'};
            end
            J = cell(1, numel(tracerNames));
            for i = 1:numel(tracerNames)
                c = model.getProp(state, lower(tracerNames{i}));
                J{i} = - T{i} .* op.Grad(c);
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
