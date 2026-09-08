classdef ChemotaxisTransmissibility < StateFunction
    % Computes face-based transmissibility for bacterial chemotaxis from chemotactic diffusivity.
    %
    % SYNOPSIS:
    %   T = ChemotaxisTransmissibility(model)
    %
    % DESCRIPTION:
    %   Converts cell-based chemotactic diffusivity D_chemo{spec} (from
    %   MicrobialChemotaxis) to face-based transmissibility T{spec}
    %   using a two-point harmonic average (via DynamicFlowTransmissibility):
    %
    %       T_{spec,f} = twoPoint(D_chemo{spec}) harmonically averaged to faces
    %
    %   Each bacterial species is treated independently.
    %
    % REQUIRED PARAMETERS:
    %   model - BiochemistryModel with MicrobialChemotaxis state function
    %
    % RETURNS:
    %   T_face - Cell array T_face{spec} of face-based transmissibilities
    %            (nInternalFaces × 1) for each bacterial species 'spec'.
    %
    % SEE ALSO:
    %   MicrobialChemotaxis, DynamicFlowTransmissibility, ChemoBactFlux

    properties
        transmissibilityComputer  % DynamicFlowTransmissibility instance
    end

    methods
        function tf = ChemotaxisTransmissibility(model, varargin)
            % Constructor.
            %
            % PARAMETERS:
            %   model - MRST biochemical model with chemotaxis enabled
            %
            % RETURNS:
            %   tf - ChemotaxisTransmissibility instance

            tf@StateFunction(model, varargin{:});
            tf = merge_options(tf, varargin{:});

            % Depend on MicrobialChemotaxis for cell-based chemotactic diffusivity
            tf = tf.dependsOn('MicrobialChemotaxis');

            % Create a DynamicFlowTransmissibility computer for harmonic averaging
            tf.transmissibilityComputer = DynamicFlowTransmissibility(model, 'dummy');

            tf.label = 'T_{chemo}';
        end

        function T_face = evaluateOnDomain(tf, model, state)
            % Convert cell-based chemotactic diffusivity to face-based transmissibility.
            %
            % PARAMETERS:
            %   model - MRST model
            %   state - Reservoir state structure
            %
            % RETURNS:
            %   T_face - Cell array T_face{spec} of face transmissibilities,
            %            one element per bacterial species.

            % Get cell-based chemotactic diffusivity (cell array per species)
            D_chemo_cell = tf.getEvaluatedDependencies(state, 'MicrobialChemotaxis');

            % Number of bacterial species
            nspec = model.biochemFluid.nbioreact;

            % Internal face count
            op = model.operators;
            nInternal = numel(op.internalConn);

            % Initialize output cell array
            T_face = cell(nspec, 1);

            % Convert each species diffusivity to face transmissibility
            for spec = 1:nspec
                D = D_chemo_cell{spec};

                % Skip if diffusivity is zero everywhere
                if isnumeric(D) && all(D(:) == 0)
                    T_face{spec} = zeros(nInternal, 1);
                    continue;
                end

                % Compute harmonic transmissibility on all faces
                T_all = tf.transmissibilityComputer.getTransmissibility(D);

                % Retain only internal faces
                T_face{spec} = T_all(op.internalConn);
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