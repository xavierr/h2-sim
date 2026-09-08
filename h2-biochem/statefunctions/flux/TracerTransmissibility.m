classdef TracerTransmissibility < StateFunction
    % Computes face-based transmissibility for mobile aqueous tracers
    % from their cell-based effective diffusivity.
    %
    % SYNOPSIS:
    %   T = TracerTransmissibility(model)
    %
    % DESCRIPTION:
    %   Converts cell-based effective diffusivity D_eff{i} (from
    %   TracerDiffusivity) to face-based transmissibility T{i} using a
    %   two-point harmonic average (via DynamicFlowTransmissibility),
    %   mirroring MicrobialTransmissibility for the bacterial species.
    %
    % REQUIRED PARAMETERS:
    %   model - BiochemistryModel with TracerDiffusivity state function
    %
    % RETURNS:
    %   T_face - Cell array of face-based transmissibilities
    %            (nInternalFaces x 1), ordered by getAqueousTracerNames.
    %
    % SEE ALSO:
    %   TracerDiffusivity, DynamicFlowTransmissibility, DiffusiveTracerFlux

    properties
        transmissibilityComputer  % DynamicFlowTransmissibility instance
    end

    methods
        function tf = TracerTransmissibility(model, varargin)
            tf@StateFunction(model, varargin{:});
            tf = merge_options(tf, varargin{:});

            tf = tf.dependsOn('TracerDiffusivity');
            tf.transmissibilityComputer = DynamicFlowTransmissibility(model, 'dummy');

            tf.label = 'T_{tracer}';
        end

        function T_face = evaluateOnDomain(tf, model, state)
            D_eff_cell = tf.getEvaluatedDependencies(state, 'TracerDiffusivity');

            op = model.operators;
            nInternal = numel(op.internalConn);

            if ismethod(model, 'getAqueousTracerNames')
                tracerNames = model.getAqueousTracerNames();
            else
                tracerNames = {'SO4', 'HS'};
            end
            T_face = cell(1, numel(tracerNames));
            for i = 1:numel(tracerNames)
                D = D_eff_cell{i};
                if isnumeric(D) && all(D(:) == 0)
                    T_face{i} = zeros(nInternal, 1);
                    continue;
                end
                T_all = tf.transmissibilityComputer.getTransmissibility(D);
                T_face{i} = T_all(op.internalConn);
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
