classdef DispersiveDiffusivity < StateFunction
    % Computes cell-based effective dispersive diffusivity.
    %
    % SYNOPSIS:
    %   d = DispersiveDiffusivity(model)
    %   d = DispersiveDiffusivity(model, 'pn', pv, ...)
    %
    % DESCRIPTION:
    %   Evaluates the cell-centred effective dispersive diffusivity for
    %   each component and phase:
    %
    %       D_eff{i, ph} = phi * S_ph * (D_disp + D_diff)
    %
    %   where:
    %       D_disp = Mechanical dispersion coefficient (phase-dependent)
    %       D_diff = Molecular diffusion coefficient (component-dependent)
    %
    %   The class checks the model flags:
    %       model.molecularDispersion  (true/false)
    %       model.molecularDiffusion   (true/false)
    %   and activates the corresponding physics automatically.
    %
    % REQUIRED PARAMETERS:
    %   model - Model instance with dispersion/diffusion flags and
    %           parameters (alphaL_water, alphaT_water, etc.)
    %
    % PROPERTIES (optional, set via property/value pairs):
    %   minPorosity = 1e-12
    %   alphaL_water = 5.0e-2   [m]
    %   alphaT_water = 5.0e-3   [m]
    %   alphaL_gas   = 1.5e-1   [m]
    %   alphaT_gas   = 1.5e-2   [m]
    %   Tref = 273.15 + 40      [K]
    %   pref = 101325           [Pa]
    %   tortuosityExponent = 7/3
    %   gasDiffExponent = 1.5
    %   defaultLiquidDiffusivity = 1e-9  [m^2/s]
    %   minDiffusivity = 1e-15           [m^2/s]
    %   defaultSigma = 3.5               [Angstrom]
    %   defaultEpsilon = 150.0           [K]
    %
    % SEE ALSO:
    %   StateFunction, DispersiveTransmissibility, ComponentPhaseDispFlux

    properties
        % --- General ---
        minPorosity = 1e-12;

        % --- Mechanical dispersion parameters ---
        alphaL_water = 5.0e-2;   % [m]
        alphaT_water = 5.0e-3;
        alphaL_gas   = 1.5e-1;
        alphaT_gas   = 1.5e-2;

        % --- Molecular diffusion parameters ---
        Tref = 273.15 + 40;               % [K]
        pref = 101325;                     % [Pa] (1 atm)
        tortuosityExponent = 7/3;          % Millington-Quirk exponent
        gasDiffExponent = 1.5;            % Chapman-Enskog temperature exponent
        defaultLiquidDiffusivity = 1e-9;   % [m^2/s]
        minDiffusivity = 1e-15;            % reasonable floor for diffusivities

        % Lennard-Jones defaults for unknown components
        defaultSigma   = 3.5;              % [Angstrom]
        defaultEpsilon = 150.0;            % [K]
    end

    methods
        %-----------------------------------------------------------------%
        function d = DispersiveDiffusivity(model, varargin)
            % Constructor.
            %
            % PARAMETERS:
            %   model - MRST model with dispersive physics enabled
            %   varargin - Optional property/value pairs
            %
            % RETURNS:
            %   d - DispersiveDiffusivity instance

            d@StateFunction(model, varargin{:});
            d = merge_options(d, varargin{:});

            % Dependencies always needed
            d = d.dependsOn({'s', 'x', 'y', 'pressure'}, 'state');

            % Mechanical dispersion depends on velocity
            if model.molecularDispersion
                d = d.dependsOn('PhaseFlux');
            end

            % Molecular diffusion depends on temperature
            if model.molecularDiffusion
                d = d.dependsOn('temperature', 'state');
            end

            % Dynamic porosity (e.g., bio-clogging)
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                d = d.dependsOn('nbact', 'state');
            end

            d.label = 'D_{eff}';
        end

        %-----------------------------------------------------------------%
        function D_eff = evaluateOnDomain(d, model, state)
            % Evaluate cell-based effective dispersive diffusivity.
            % Separates mechanical dispersion (phase-only) from diffusion (component+phase)
            % to avoid redundant computations.
            %
            % PARAMETERS:
            %   model - MRST model
            %   state - Reservoir state structure
            %
            % RETURNS:
            %   D_eff - Cell array D_eff{i, ph} containing cell-based
            %           effective diffusivity for component i in phase ph
            %           (ncells x 1), units m^2/s.

            if ~isfield(state, 'x')
                ncomp = model.getNumberOfComponents();
                nph = model.getNumberOfPhases();
                D_eff = cell(ncomp, nph);
                [D_eff{:}] = deal(0);
                return;
            end

            op = model.operators;
            G = model.G;

            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            nm = model.getPhaseNames();
            ncells = G.cells.num;

            % --- Porosity (cell-based) ----------------------------------------
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                nbactArray = model.extractBactValues(nbact);
                phi = model.rock.poro(p, nbactArray{:}); % Apply both modifications
            else
                phi = model.rock.poro;
            end
            phi = max(phi, d.minPorosity);

            % --- Pre-compute mechanical dispersion (phase-only) ---------------
            % D_disp{ph} computed ONCE per phase, not per component
            D_disp_ph = cell(nph, 1);
            [D_disp_ph{:}] = deal(zeros(ncells, 1));

            if model.molecularDispersion
                phase_flux = model.getProp(state, 'PhaseFlux');

                % Dispersivity coefficients
                aL = zeros(1, nph);
                aT = zeros(1, nph);
                aL(L_ix) = d.alphaL_water;
                aT(L_ix) = d.alphaT_water;
                if V_ix ~= L_ix
                    aL(V_ix) = d.alphaL_gas;
                    aT(V_ix) = d.alphaT_gas;
                end

                % Compute dispersion for each phase (ONCE)
                for ph = 1:nph
                    u_face = ifcell(phase_flux, ph);
                    internalFaces = find(op.internalConn);
                    internalFaceMap = sparse(internalFaces, ...
                        1:numel(internalFaces), 1, G.faces.num, ...
                        numel(internalFaces));
                    u_face_all = internalFaceMap*u_face;

                    v_mag = faceFlux2cellSpeed(G, u_face_all);

                    % Isotropic dispersive diffusivity (Scheidegger)
                    % D_disp = aL*|v| for longitudinal
                    % Full tensor: D = aL*(vv'/|v|^2) + aT*(I - vv'/|v|^2)
                    % Isotropic approximation: D = ((aL + 2*aT)/3)*|v|
                    D_disp_ph{ph} = ((aL(ph) + 2*aT(ph))/3) .* v_mag;
                    D_disp_ph{ph} = max(D_disp_ph{ph}, 0);
                end
            end

            % --- Pre-compute molecular diffusion parameters -------------------
            if model.molecularDiffusion
                [p, T] = model.getProps(state, 'pressure', 'temperature');
                p_safe = max(p, 1e-8*barsa);
                gasScale = (T./d.Tref).^d.gasDiffExponent .* (d.pref./p_safe);
                [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(d, model);
                Dij_ref = localBinaryDiffusionReference(d, model, paramLJ);
            end

            % --- Initialize output -------------------------------------------
            D_eff = cell(ncomp, nph);
            [D_eff{:}] = deal(0);

            % --- Phase loop --------------------------------------------------
            for ph = 1:nph
                s = model.getProp(state, ['s', nm(ph)]);
                s = max(s, 0);

                % Get dispersion for this phase (already computed)
                D_disp = D_disp_ph{ph};

                % >>> Molecular Diffusion (component-dependent) <<<
                if model.molecularDiffusion && (ph == L_ix || ph == V_ix)
                    % Millington-Quirk tortuosity (same for all components in this phase)
                    phiS = phi .* s;
                    tau_MQ = (phiS).^(d.tortuosityExponent) .* (phi.^(-2));
                    tau_MQ = max(tau_MQ, 0);

                    % Component loop for diffusion only
                    for c = 1:ncomp
                        if ph == L_ix
                            % Liquid diffusion
                            D_i = Dliq_ref(c);
                        elseif ph == V_ix
                            % Gas diffusion (Wilke multi-component)
                            % D_i,mix = (1 - y_i) / sum(y_j/D_ij) for j!=i
                            yAll = model.getProp(state, 'y');
                            yc_curr = localGetComponentVector(yAll, c);
                            numerator = max(1 - yc_curr, 1e-6);

                            invDi = 0;
                            for j = 1:ncomp
                                if j == c, continue; end
                                yj = localGetComponentVector(yAll, j);
                                Dij = Dij_ref(c, j) .* gasScale;
                                Dij = max(Dij, d.minDiffusivity);
                                invDi = invDi + yj ./ Dij;
                            end
                            Di_mix = numerator ./ max(invDi, d.minDiffusivity);
                            D_i = Di_mix;
                        end

                        D_i = max(D_i, d.minDiffusivity);

                        % Apply tortuosity
                        D_diff = tau_MQ .* D_i;

                        % Combine: phi*s*(D_disp + D_diff) per component
                        D_eff{c, ph} = phi .* s .* (D_disp + D_diff);
                    end
                else
                    % Dispersion only (no diffusion) - reuse D_disp for all components
                    for c = 1:ncomp
                        D_eff{c, ph} = phi .* s .* D_disp;
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
if iscell(field)
    val = field{ph};
else
    val = field(:, ph);
end
end

function v = localGetComponentVector(z, j)
% Extract component j from cell array or matrix.
if iscell(z)
    v = z{j};
else
    v = z(:, j);
end
end

function [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(d, model)
% Load liquid diffusion coefficients and Lennard-Jones parameters.
% Returns diffusivity in m^2/s and LJ parameters [sigma (Angstrom), epsilon (K)]
namecp = model.compFluid.names();
ncomp  = numel(namecp);
Dliq_ref = d.defaultLiquidDiffusivity * ones(ncomp, 1);
paramLJ  = repmat([d.defaultSigma, d.defaultEpsilon], ncomp, 1);

db.Dliq = struct('h2',6.44e-9,'c1',2.15e-9,'methane',2.15e-9,'co2',2.72e-9, ...
    'h2o',3.29e-9,'water',3.29e-9,'n2',2.86e-9,'c2',1.72e-9, ...
    'ethane',1.72e-9,'c3',1.43e-9,'propane',1.43e-9,'h2s',2.15e-9, ...
    'nc4',1.15e-9,'butane',1.15e-9);
db.LJ = struct('h2',[2.92,59.7],'c1',[3.758,148.6],'methane',[3.758,148.6], ...
    'co2',[3.996,195.2],'h2o',[2.641,809.1],'water',[2.641,809.1], ...
    'n2',[3.798,71.4],'c2',[4.443,215.7],'ethane',[4.443,215.7], ...
    'c3',[5.118,237.1],'propane',[5.118,237.1],'h2s',[3.60,301.0], ...
    'nc4',[5.206,289.5],'butane',[5.206,289.5]);

for i = 1:ncomp
    key = lower(namecp{i});
    if isfield(db.Dliq, key)
        Dliq_ref(i) = db.Dliq.(key);
    end
    if isfield(db.LJ, key)
        paramLJ(i,:) = db.LJ.(key);
    end
end
end

function Dij_ref = localBinaryDiffusionReference(d, model, paramLJ)
% Computes binary gas diffusion coefficients using Chapman-Enskog theory
% with Lennard-Jones collision integral (SI units, m^2/s at reference T/P).
%
% Formula: D_ij = (1.858e-7 * T^1.5) / (sigma_ij^2 * sqrt(Mij*) * Omega_D)
% where Omega_D is the collision integral (Neufeld fit)

ncomp = size(paramLJ, 1);
Molmass = 1e3 .* model.compFluid.molarMass;   % g/mol (converted from kg/mol)
epsilon = paramLJ(:, 2);                       % Lennard-Jones energy [K]
sigma = paramLJ(:, 1);                         % collision diameter [Angstrom]

const = 1.858e-7;  % Chapman-Enskog constant (m^2/s)
fTref = d.Tref^1.5;   % temperature factor at reference conditions
Dij_ref = zeros(ncomp);

for i = 1:ncomp
    for j = 1:ncomp
        if i == j
            Dij_ref(i,j) = inf;
            continue;
        end

        % Combining rules (Lorentz-Berthelot)
        sigma_ij = 0.5 * (sigma(i) + sigma(j));
        epsilon_ij = sqrt(epsilon(i) * epsilon(j));

        % Reduced temperature for collision integral
        T_star = d.Tref / epsilon_ij;

        % Neufeld collision integral approximation
        % Omega_D(T*) = A/(T*)^B + C/exp(D*T*) + E/exp(F*T*) + G/exp(H*T*)
        A = 1.06036;
        B = 0.15610;
        C = 0.19300;
        D = 0.47635;
        E = 1.03587;
        F = 1.52996;
        G = 1.76474;
        H = 3.89411;

        Omega_D = A * (T_star^(-B)) + C * exp(-D*T_star) + ...
            E * exp(-F*T_star) + G * exp(-H*T_star);

        % Molecular mass term
        sqrtMass = sqrt((1.0 / Molmass(i)) + (1.0 / Molmass(j)));
        sigma_ij2 = sigma_ij^2;

        % Chapman-Enskog formula with collision integral
        Dij_ref(i,j) = (const * fTref * sqrtMass) / (sigma_ij2 * Omega_D);
    end
end

% Symmetrize matrix: D_ij = D_ji
Dij_ref = 0.5 * (Dij_ref + Dij_ref.');
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
------------------------------------------------------------------------------
%}