classdef ComponentPhaseDispFlux < StateFunction
    % Computes face-based dispersive + diffusive flux for components.
    %
    % SYNOPSIS:
    %   flux = ComponentPhaseDispFlux(model)
    %   flux = ComponentPhaseDispFlux(model, 'pn', pv, ...)
    %
    % DESCRIPTION:
    %   Computes the face-wise dispersive/diffusive flux for each component
    %   in each phase:
    %
    %       j*_{i,α} = -ρ_α^f · D_{eff,i,α}^f · ∇w_{i,α}
    %       J_{i,α}  = j*_{i,α} - w_{i,α}^f Σ_j j*_{j,α}
    %
    %   where w is component mass fraction and D_{eff,i,α}^f is the
    %   face-averaged effective dispersive/diffusive diffusivity (from
    %   DispersiveDiffusivity cell values). The correction velocity in the
    %   second expression ensures that diffusion does not create a net
    %   phase mass flux.
    %
    %   D_{eff,i,α} = φ · S_α · (D_{disp,α} + D_{diff,i,α})
    %
    %   - D_disp: Mechanical dispersion (phase-dependent, component-independent)
    %   - D_diff: Molecular diffusion (component and phase-dependent)
    %
    % REQUIRED PARAMETERS:
    %   model - Model with `DispersiveDiffusivity` state function registered.
    %
    % SEE ALSO:
    %   StateFunction, DispersiveDiffusivity

    properties
    end

    methods
        %-----------------------------------------------------------------%
        function df = ComponentPhaseDispFlux(model, varargin)
            % Constructor.
            %
            % PARAMETERS:
            %   model - MRST model with dispersive physics enabled
            %   varargin - Optional property/value pairs
            %
            % RETURNS:
            %   df - ComponentPhaseDispFlux instance

            df@StateFunction(model, varargin{:});
            df = merge_options(df, varargin{:});

            % Dependencies
            df = df.dependsOn('DispersiveTransmissibility');
            df = df.dependsOn('ComponentPhaseMassFractions', ...
                'PVTPropertyFunctions');
            df = df.dependsOn('Density', 'PVTPropertyFunctions');

            df.label = 'J_{i,\\alpha}^{disp+diff}';
        end

        %-----------------------------------------------------------------%
        function J = evaluateOnDomain(df, model, state)
            % Evaluate the dispersive/diffusive flux on the domain.
            %
            % PARAMETERS:
            %   model - MRST model
            %   state - Reservoir state structure
            %
            % RETURNS:
            %   J - Cell array J{i, ph} containing face-based flux values
            %       for component i in phase ph (nfac × 1)

            op = model.operators;

            % Get number of components and phases
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            % Retrieve dispersive transmissibility (component-phase-face)
            % Expected: T_disp{i, ph} (nfac × 1)
            T_disp = df.getEvaluatedDependencies(state, 'DispersiveTransmissibility');

            % ComponentPhaseFlux is a mass flux, so diffusion must use
            % mass fractions rather than the EOS mole fractions.
            massFraction = model.getProp(state, 'ComponentPhaseMassFractions');

            % Get density
            rho = model.getProp(state, 'Density');

            % Initialize output
            J = cell(ncomp, nph);
            [J{:}] = deal(0);

            % Phase loop
            for ph = 1:nph
                % Face-average density for this phase
                rho_ph = ifcell(rho, ph);
                rho_f = op.faceAvg(rho_ph);

                rawFlux = cell(ncomp, 1);
                totalRawFlux = 0;

                % Compute uncorrected mass-fraction Fickian fluxes.
                for c = 1:ncomp
                    T_c_ph = T_disp{c, ph};
                    w = massFraction{c, ph};
                    if isempty(w)
                        rawFlux{c} = 0;
                    else
                        rawFlux{c} = -rho_f .* T_c_ph .* op.Grad(w);
                    end
                    totalRawFlux = totalRawFlux + rawFlux{c};
                end

                % Diffusion is internal redistribution. Correct the
                % component fluxes to ensure sum_c J_{c,ph} = 0.
                for c = 1:ncomp
                    w = massFraction{c, ph};
                    if isempty(w)
                        J{c, ph} = rawFlux{c};
                    else
                        J{c, ph} = rawFlux{c} - ...
                            op.faceAvg(w).*totalRawFlux;
                    end
                end
            end
        end
    end
end

% ------------------------------------------------------------------------
% Helper functions
% ------------------------------------------------------------------------
function val = ifcell(field, ph)
    % Extract phase from cell array or matrix.
    %
    % PARAMETERS:
    %   field - Cell array or matrix
    %   ph    - Phase index
    %
    % RETURNS:
    %   val - Phase-specific values

    if iscell(field)
        val = field{ph};
    else
        val = field(:, ph);
    end
end

%{
------------------------------------------------------------------------------
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
------------------------------------------------------------------------------
%}