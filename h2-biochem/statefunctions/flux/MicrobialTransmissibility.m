classdef MicrobialTransmissibility < StateFunction
    % Computes face-based transmissibility for microbial species from diffusivity.
    %
    % SYNOPSIS:
    %   T = MicrobialTransmissibility(model)
    %
    % DESCRIPTION:
    %   Converts cell-based microbial diffusivity D_eff{spec} (from
    %   MicrobialDiffusivity) to face-based transmissibility T{spec}
    %   using a two-point harmonic average (via DynamicFlowTransmissibility).
    %
    %       T_{spec,f} = twoPoint(D_eff{spec}) harmonically averaged to faces
    %
    %   Each bacterial species is treated as an independent diffusing
    %   quantity.
    %
    % REQUIRED PARAMETERS:
    %   model - BiochemistryModel with MicrobialDiffusivity state function
    %
    % RETURNS:
    %   T_face - Cell array T_face{spec} of face-based transmissibilities
    %            (nInternalFaces × 1) for each bacterial species 'spec'.
    %
    % SEE ALSO:
    %   MicrobialDiffusivity, DynamicFlowTransmissibility, BactFlux

    properties
        transmissibilityComputer  % DynamicFlowTransmissibility instance
    end

    methods
        function tf = MicrobialTransmissibility(model, varargin)
            % Constructor.
            %
            % PARAMETERS:
            %   model - MRST biochemical model with bacterial diffusion enabled
            %
            % RETURNS:
            %   tf - MicrobialTransmissibility instance

            tf@StateFunction(model, varargin{:});
            tf = merge_options(tf, varargin{:});

            % Depend on MicrobialDiffusivity for cell-based effective diffusivity
            tf = tf.dependsOn('MicrobialDiffusivity');

            % Create a DynamicFlowTransmissibility computer for harmonic averaging
            tf.transmissibilityComputer = DynamicFlowTransmissibility(model, 'dummy');

            tf.label = 'T_{bact}';
        end

        function T_face = evaluateOnDomain(tf, model, state)
            % Convert cell-based diffusivity to face-based transmissibility.
            %
            % PARAMETERS:
            %   model - MRST model
            %   state - Reservoir state structure
            %
            % RETURNS:
            %   T_face - Cell array T_face{spec} of face transmissibilities,
            %            one element per bacterial species.

            % Get cell-based microbial diffusivity (cell array per species)
            D_eff_cell = tf.getEvaluatedDependencies(state, 'MicrobialDiffusivity');

            % Number of bacterial species
            nspec = model.biochemFluid.nbioreact;

            % Internal face count
            op = model.operators;
            nInternal = numel(op.internalConn);

            % Initialize output cell array
            T_face = cell(nspec, 1);

            % Convert each species diffusivity to face transmissibility
            for spec = 1:nspec
                D = D_eff_cell{spec};

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