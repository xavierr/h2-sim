classdef DiffusiveBactFlux < StateFunction
    % Diffusive flux of bacteria (mass‑fraction based).
    %
    % SYNOPSIS:
    %   flux = DiffusiveBactFlux(model)
    %
    % DESCRIPTION:
    %   Computes the face-wise diffusive flux of bacteria in the liquid
    %   phase as
    %       J = - rho_l^f * T * Grad(nbact)
    %   where T is the microbial transmissibility provided by
    %   `MicrobialTransmissibility` and ρ_l^f is the face‑averaged liquid
    %   density.
    %
    % REQUIRED PARAMETERS:
    %   model - Model with `MicrobialTransmissibility` state function
    %           registered and bacterial modelling enabled.
    %
    % SEE ALSO:
    %   StateFunction, MicrobialTransmissibility

    properties
    end

    methods
        %-----------------------------------------------------------------%
        function df = DiffusiveBactFlux(model, varargin)
            df@StateFunction(model, varargin{:});
            df = merge_options(df, varargin{:});

            % Dependencies
            df = df.dependsOn('MicrobialTransmissibility');
            df = df.dependsOn('nbact', 'state');
            df = df.dependsOn('Density', 'PVTPropertyFunctions');

            df.label = 'J_{bact}^{diff}';
        end

        %-----------------------------------------------------------------%
        function DN = evaluateOnDomain(prop, model, state)
            op = model.operators;
            bcrm=model.biochemFluid;
            nbioreact=bcrm.nbioreact;

            % Retrieve microbial transmissibility
            T = prop.getEvaluatedDependencies(state, 'MicrobialTransmissibility');

            % State variables and properties
            [nbact, rho] = model.getProps(state, 'nbact', 'Density');

            % Initialize with zero growth rate
            DN=cell(1,nbioreact);
            [DN{:}] = deal(0);

            % Liquid density (face averaged)
            for i=1:nbioreact
                L_ix = model.getLiquidIndex();
                if iscell(rho)
                    rhoL = rho{L_ix};
                else
                    rhoL = rho(:, L_ix);
                end
                if iscell(nbact)
                    nbacti=nbact{i};
                else
                    nbacti=nbact(:,i);
                end
                rhoLf = model.operators.faceAvg(rhoL);

                % Diffusive flux (all faces)
                DN{i} = - rhoLf .* T{i} .* op.Grad(nbacti);
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