classdef TracerDiffusivity < StateFunction
    % Computes cell-based effective dispersive/diffusive diffusivity for
    % mobile non-volatile aqueous tracers.
    %
    % SYNOPSIS:
    %   d = TracerDiffusivity(model)
    %
    % DESCRIPTION:
    %   Mirrors DispersiveDiffusivity (used for the volatile EOS
    %   components), but for liquid-only aqueous tracers, which are
    %   never part of the EOS mole-fraction vectors (see
    %   BiochemistryModel.sulfateReduction / AqueousTracerMass):
    %
    %       D_eff{i} = phi * S_l * (D_disp + D_diff_i)
    %
    %   D_disp is the same isotropic mechanical-dispersion approximation
    %   as DispersiveDiffusivity (liquid dispersivities only, since the
    %   tracers only exist in the liquid phase). D_diff_i uses fixed
    %   reference aqueous molecular diffusivities for SO4^2- and HS-,
    %   scaled by Millington-Quirk tortuosity.
    %
    % REQUIRED PARAMETERS:
    %   model - BiochemistryModel with sulfateReduction enabled and
    %           molecularDiffusion/molecularDispersion flags set.
    %
    % SEE ALSO:
    %   DispersiveDiffusivity, TracerTransmissibility, DiffusiveTracerFlux

    properties
        minPorosity = 1e-12;

        % --- Mechanical dispersion (liquid phase only) ---
        alphaL_water = 5.0e-2;   % [m]
        alphaT_water = 5.0e-3;   % [m]

        % --- Molecular diffusion ---
        tortuosityExponent = 7/3;   % Millington-Quirk exponent
        minDiffusivity = 1e-15;    % reasonable floor for diffusivities

        % Reference aqueous molecular diffusivities at infinite dilution
        % [m^2/s]. These are configurable approximate values; PHREEQC
        % supplies speciation but MRST transports the analytical totals.
        D_SO4 = 1.07e-9;
        D_HS  = 1.73e-9;
        D_HCO3 = 1.18e-9;
        D_CA = 0.79e-9;
        D_MG = 0.71e-9;
    end

    methods
        %-----------------------------------------------------------------%
        function d = TracerDiffusivity(model, varargin)
            d@StateFunction(model, varargin{:});
            d = merge_options(d, varargin{:});

            d = d.dependsOn('s', 'state');

            if model.molecularDispersion
                d = d.dependsOn('PhaseFlux');
            end

            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                d = d.dependsOn('nbact', 'state');
            end

            d.label = 'D_{tracer}^{eff}';
        end

        %-----------------------------------------------------------------%
        function D_eff = evaluateOnDomain(d, model, state)
            G = model.G;
            L_ix = model.getLiquidIndex();

            % --- Porosity (cell-based) ---
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                nbactArray = model.extractBactValues(nbact);
                phi = model.rock.poro(p, nbactArray{:});
            else
                phi = model.rock.poro;
            end
            phi = max(phi, d.minPorosity);

            % --- Liquid saturation ---
            s = model.getProp(state, 's');
            if iscell(s)
                sL = max(s{L_ix}, 0);
            else
                sL = max(s(:, L_ix), 0);
            end

            % --- Mechanical dispersion (liquid phase) ---
            D_disp = zeros(G.cells.num, 1);
            if model.molecularDispersion
                op = model.operators;
                phase_flux = model.getProp(state, 'PhaseFlux');
                u_face = ifcell(phase_flux, L_ix);
                internalFaces = find(op.internalConn);
                internalFaceMap = sparse(internalFaces, ...
                    1:numel(internalFaces), 1, G.faces.num, ...
                    numel(internalFaces));
                u_face_all = internalFaceMap*u_face;

                v_mag = faceFlux2cellSpeed(G, u_face_all);

                D_disp = ((d.alphaL_water + 2*d.alphaT_water)/3) .* v_mag;
                D_disp = max(D_disp, 0);
            end

            % --- Molecular diffusion (Millington-Quirk tortuosity) ---
            if ismethod(model, 'getAqueousTracerNames')
                tracerNames = model.getAqueousTracerNames();
            else
                tracerNames = {'SO4', 'HS'};
            end
            D_eff = cell(1, numel(tracerNames));
            if model.molecularDiffusion
                phiS = phi .* sL;
                tau_MQ = max((phiS).^(d.tortuosityExponent) .* (phi.^(-2)), 0);
            end
            for i = 1:numel(tracerNames)
                Dref = d.(['D_', upper(tracerNames{i})]);
                if model.molecularDiffusion
                    D_eff{i} = phi .* sL .* (D_disp + ...
                        tau_MQ .* max(Dref, d.minDiffusivity));
                else
                    D_eff{i} = phi .* sL .* D_disp;
                end
            end
        end
    end
end

% ------------------------------------------------------------------------
function val = ifcell(field, ph)
if iscell(field)
    val = field{ph};
else
    val = field(:, ph);
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
