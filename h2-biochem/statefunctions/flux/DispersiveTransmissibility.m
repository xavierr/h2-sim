classdef DispersiveTransmissibility < StateFunction
    % Computes face-based dispersive/diffusive transmissibility from diffusivity.
    %
    % SYNOPSIS:
    %   T = DispersiveTransmissibility(model)
    %
    % DESCRIPTION:
    %   Converts cell-based effective dispersive/diffusive diffusivity
    %   D_eff{c,ph} (from DispersiveDiffusivity) to face-based transmissibility
    %   T{c,ph} using DynamicFlowTransmissibility with two-point flux approximation:
    %
    %       T_{c,α,f} = twoPoint(D_{c,α}) harmonically averaged to faces
    %
    %   Each component/phase diffusivity is treated as a conductivity field.
    %
    % REQUIRED PARAMETERS:
    %   model - Model with DispersiveDiffusivity state function
    %
    % RETURNS:
    %   T - Cell array T{c, ph} of face-based transmissibilities
    %       (nfac_internal × 1) for each component c and phase ph
    %
    % SEE ALSO:
    %   DispersiveDiffusivity, DynamicFlowTransmissibility, ComponentPhaseDispFlux

    properties
        transmissibilityComputer  % DynamicFlowTransmissibility instance
    end

    methods
        function tf = DispersiveTransmissibility(model, varargin)
            % Constructor.
            %
            % PARAMETERS:
            %   model - MRST model with dispersive physics enabled
            %
            % RETURNS:
            %   tf - DispersiveTransmissibility instance

            tf@StateFunction(model, varargin{:});
            tf = merge_options(tf, varargin{:});

            % Depend on DispersiveDiffusivity for cell-based effective diffusivity
            tf = tf.dependsOn('DispersiveDiffusivity');

            % Create a DynamicFlowTransmissibility computer for computing face transmissibilities
            tf.transmissibilityComputer = DynamicFlowTransmissibility(model, 'dummy');

            tf.label = 'T_{c,\\alpha}';
        end

        function T_face = evaluateOnDomain(tf, model, state)
            % Convert cell-based diffusivity to face-based transmissibility.
            % Uses cached two-point and harmonic-averaging operators for efficiency.
            %
            % PARAMETERS:
            %   model - MRST model
            %   state - Reservoir state structure
            %
            % RETURNS:
            %   T_face - Cell array T_face{c, ph} of face-based transmissibilities
            %            (nfac_internal × 1) for each component/phase pair

            % Get cell-based effective diffusivity
            D_eff_cell = tf.getEvaluatedDependencies(state, 'DispersiveDiffusivity');

            % Get number of components and phases
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();

            % Initialize output cell array
            T_face = cell(ncomp, nph);

            % Get operators (cached in transmissibilityComputer)
            op = model.operators;
            nfac_internal = numel(op.internalConn);

            % Convert each component/phase pair efficiently
            for c = 1:ncomp
                for ph = 1:nph
                    D_c_ph = D_eff_cell{c, ph};

                    % Handle zero/empty diffusivity
                    if isnumeric(D_c_ph) && all(D_c_ph(:) == 0)
                        T_face{c, ph} = zeros(nfac_internal, 1);
                        continue;
                    end

                    % Use cached getTransmissibility which avoids recreating operators
                    T_all = tf.transmissibilityComputer.getTransmissibility(D_c_ph);

                    % Extract internal connections only
                    T_face{c, ph} = T_all(op.internalConn);
                end
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
