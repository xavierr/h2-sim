classdef ChemotaxisBactFlux < StateFunction
    % Diffusive flux of bacteria (mass‑fraction based).
    %
    % SYNOPSIS:
    %   flux = ChemotaxisBactFlux(model)
    %
    % DESCRIPTION:
    %   Computes the face-wise diffusive flux of bacteria in the liquid
    %   phase as
    %Jchemotaxis = Kchemo * nbact *(1 - nbact) * Grad(xH2)
    % (volume-filling effect, generalized Keller-Segel model)
    %   where T is the chemotaxis transmissibility provided by
    %   `chemotaxisTransmissibility` and rho_l^f is the face-averaged liquid
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
        function df = ChemotaxisBactFlux(model, varargin)
            df@StateFunction(model, varargin{:});
            df = merge_options(df, varargin{:});

            % Dependencies
            df = df.dependsOn('ChemotaxisTransmissibility');
            df = df.dependsOn('x', 'state');
            df = df.dependsOn('Density', 'PVTPropertyFunctions');

            df.label = 'J_{bact}^{chemo}';
        end

        %-----------------------------------------------------------------%
        function DN = evaluateOnDomain(prop, model, state)
            op = model.operators;
            bcrm=model.biochemFluid;
            nbioreact=bcrm.nbioreact;


            % Retrieve microbial transmissibility
            T = prop.getEvaluatedDependencies(state, 'ChemotaxisTransmissibility');

            % State variables and properties
            [x, rho] = model.getProps(state, 'x', 'Density');

            % Initialize with zero growth rate
            DN=cell(1,nbioreact);
            [DN{:}] = deal(0);

            % Liquid density (face averaged)
            namecp = model.getComponentNames();
            L_ix = model.getLiquidIndex();
            if iscell(rho)
                rhoL = rho{L_ix};
            else
                rhoL = rho(:, L_ix);
            end
            rhoLf = model.operators.faceAvg(rhoL);

            for i=1:nbioreact
                indH2= find(strcmpi(namecp, bcrm.rH2(i)), 1);
                if iscell(x)
                    xH2=x{indH2};
                else
                    xH2=x(:,indH2);
                end

                % Diffusive flux (all faces)
                DN{i} = + rhoLf .* T{i} .* op.Grad(xH2);
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