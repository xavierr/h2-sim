classdef BacterialMass < StateFunction & ComponentProperty
    % Bacterial mass per cell accounting for pore volume and phase saturation
    %
    % SYNOPSIS:
    %   bm = BacterialMass(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the mass of bacterial cells in each grid cell, accounting for:
    %   - Bacterial concentration (nbact)
    %   - Liquid phase saturation
    %   - Pore volume
    %   - Phase densities
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with bacterial modeling enabled
    %
    % OPTIONAL PARAMETERS:
    %   'minimumDerivatives' - Minimum absolute values for derivatives
    %
    % RETURNS:
    %   Class instance ready for use in simulation
    %
    % SEE ALSO:
    %   CompositionalModel, EquationsCompositional

    properties (Access = protected)
        minimumDerivatives = 1.0e-3; % Minimum derivative cutoff value
    end

    methods
        function gp = BacterialMass(model, varargin)
            % Constructor for bacterial mass calculator

            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('nbact', 'state');          % Bacterial concentration
            gp = gp.dependsOn('s', 'state');              % Phase saturations
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            gp.label = 'M_{bio}';
        end

        function mb = evaluateOnDomain(prop, model, state)
            % Compute bacterial mass in each grid cell
            %
            % Theory: Mass = pore_volume * saturation * liquid_density *
            % bacterial mass fraction.
            % Units: [m³] * [dimensionless] * [kg/m³] * [dimensionless] = [kg]
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct with fields
            %
            % RETURNS:
            %   mb    - Bacterial mass per cell [kg]

            % Get required properties from state
            s = model.getProp(state, 's');
            nbact = model.getProp(state, 'nbact');

            % Get pore volume
            pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
            rho = model.PVTPropertyFunctions.get(model, state, 'Density');
            bcm=model.biochemFluid;
            nbioreact=bcm.nbioreact;
            mb=cell(1,nbioreact);
            [mb{:}] = deal(0);
            % Get liquid saturation
            L_ix = model.getLiquidIndex();
            if iscell(s)
                sL = max(s{L_ix}, 1.0e-8);
            else
                sL = max(s(:, L_ix), 1.0e-8);
            end

            if iscell(rho)
                rhoL = rho{L_ix};
            else
                rhoL = rho(:, L_ix);
            end

            for i=1:nbioreact
                % Compute bacterial mass: pv * S_l * rho_l * nbact
                % where nbact is the model's dimensionless bacterial mass
                % fraction and rhoL converts it to a mass concentration.
                if iscell(nbact)
                    nbacti=nbact{i};
                else
                    nbacti=nbact(:,i);
                end
                mb{i} = pv .* sL .* rhoL .* nbacti;

                % Ensure minimum derivatives for stability
                mb{i} = ensureMinimumDerivatives(prop, model, mb{i});
            end
        end

        function prop = setMinimumDerivatives(prop, der)
            if nargin < 2
                % Just set to a very small value, since nothing was
                % specified.
                der = 1e-10;
            end
            prop.minimumDerivatives = der;
        end

        function der = getMinimumDerivatives(prop)
            der = prop.minimumDerivatives;
        end

        function mass = ensureMinimumDerivatives(prop, model, mass)
            der = prop.minimumDerivatives;
            if isempty(der)
                return;
            end
            nc = numel(mass);
            if numel(der) == 1
                der = repmat(der, 1, nc);
            end
            isDiag = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
            if isDiag
                rowMajor = model.AutoDiffBackend.rowMajor;
                for c = 1:nc
                    m = mass;
                    if isnumeric(m) || size(m.jac{1}.diagonal, 2 - rowMajor) < c
                        continue
                    end
                    d = der(c);
                    if rowMajor
                        bad = abs(m.jac{1}.diagonal(c, :)) < d;
                        m.jac{1}.diagonal(c, bad) = d;
                    else
                        bad = abs(m.jac{1}.diagonal(:, c)) < d;
                        m.jac{1}.diagonal(bad, c) = d;
                    end
                    mass = m;
                end
            else
                for c = 1:nc
                    m = mass;
                    d = der(c);
                    if isnumeric(m)
                        continue;
                    end
                    if numel(m.jac) < c
                        % Not matching number of derivatives - stop early
                        break;
                    end
                    J = m.jac{c};
                    [n, l] = size(J);
                    if n ~= l
                        continue;
                    end
                    diagonal = diag(J);
                    bad = abs(diagonal) < d;
                    if any(bad)
                        m.jac{c} = m.jac{c} + sparse(1:n, 1:n, bad.*d, n, n);
                    end
                    mass = m;
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