classdef BactConvertionRate < StateFunction
    % Bacterial conversion rate model for compositional simulations.
    %
    % SYNOPSIS:
    %   bcr = BactConvertionRate(model)
    %   bcr = BactConvertionRate(model, 'pn', pv, ...)
    %
    % DESCRIPTION:
    %   Computes component conversion rates due to bacterial growth for
    %   methanogenic archaea. The conversion rate follows the stoichiometry
    %   of the metabolic reaction:
    %
    %       4 H2 + CO2 -> CH4 + 2 H2O
    %
    %   The conversion rate for each component c is computed as:
    %
    %       qbiot(c) = gamma_norm(c) * CarbonLimitedGrowthRate * bmass / Y_H2
    %
    %   where:
    %       gamma_norm = mass-weighted stoichiometric coefficient,
    %                    normalized by hydrogen and nbactMax
    %       bmass      = Bacterial mass in the liquid pore volume [kg]
    %       Y_H2       = Reaction yield scale
    %
    %   The stoichiometric coefficients are normalized by the H2 coefficient
    %   and scaled by the bacterial carrying capacity (nbactMax).
    %
    % REQUIRED PARAMETERS:
    %   model - Model instance with:
    %           .bacteriamodel (true)
    %           .biochemFluid.metabolicReaction ('MethanogenicArchae')
    %           .biochemFluid.rH2 (reactant component name, e.g., 'H2')
    %           .biochemFluid.rsub (substrate component name, e.g., 'CO2')
    %           .biochemFluid.Y_H2 (yield coefficient)
    %           .biochemFluid.nbactMax (carrying capacity)
    %           .gammak (stoichiometric coefficients)
    %
    % RETURNS:
    %   qbiot - Cell array of conversion rates for each component [kg/s]
    %
    % SEE ALSO:
    %   StateFunction, BiochemistryModel, ComponentSource

    properties
        % No additional properties needed - all parameters come from model
    end

    methods
        %-----------------------------------------------------------------%
        function bcr = BactConvertionRate(model, varargin)
            % Constructor for bacterial conversion rate calculator.
            %
            % PARAMETERS:
            %   model - Compositional model with bacterial physics
            %   varargin - Optional property/value pairs
            %
            % RETURNS:
            %   bcr - BactConvertionRate instance

            bcr@StateFunction(model, varargin{:});

            % This state function is also registered on the facility model.
            rm = model;
            if isprop(model, 'ReservoirModel') && ~isempty(model.ReservoirModel)
                rm = model.ReservoirModel;
            end
            if isa(rm, 'BiochemistryPhreeqcModel')
                growthRateName = 'CarbonLimitedGrowthRate';
            else
                growthRateName = 'PsiGrowthRate';
            end
            bcr = bcr.dependsOn(growthRateName, 'state');
            if isprop(model, 'phreeqcTimestepCoupling') && model.phreeqcTimestepCoupling
                bcr = bcr.dependsOn('phreeqcPH', 'state');
            end

            bcr.label = 'Q_biot';
        end

        %-----------------------------------------------------------------%
        function qbiot = evaluateOnDomain(bcr, model, state)
            % Compute bacterial conversion rate for each component.
            %
            % PARAMETERS:
            %   model - Compositional model instance
            %   state - Reservoir state structure
            %
            % RETURNS:
            %   qbiot - Cell array qbiot{c} containing conversion rates
            %           for component c (ncells × 1), units [kg/s]

            % Get reservoir model and number of components
            if isprop(model, 'ReservoirModel') && ~isempty(model.ReservoirModel)
                rm = model.ReservoirModel;
            else
                rm = model;
            end
            ncomp = rm.EOSModel.getNumberOfComponents();

            % Initialize output with zeros
            qbiot = cell(ncomp, 1);
            [qbiot{:}] = deal(0);

            % Early return if bacterial model is not active
            if ~(rm.bacteriamodel && rm.liquidPhase)
                return;
            end
            if ismethod(rm, 'isSequentialCompositionalPhreeqcBackend') && ...
                    rm.isSequentialCompositionalPhreeqcBackend()
                % The compositional PHREEQC COM split owns MET/ACE/SRB chemistry after
                % convergence. Keep the MRST transport equations, but do
                % not apply their microbial component sources as well.
                return;
            end

            % Get biochemical fluid properties
            bcrm = rm.biochemFluid;
            nbioreact = bcrm.nbioreact;

            % Validate reactor types
            for i = 1:nbioreact
                rxn = bcrm.metabolicReaction{i};
                if ~ismember(rxn, {'MethanogenicArchae', 'AcetogenicBacteria', 'SulfateReducingBacteria'})
                    error('BactConvertionRate:UnsupportedReaction', ...
                        'Unsupported metabolic reaction: %s', rxn);
                end
            end

            % Get component names
            compNames = rm.EOSModel.getComponentNames();

            % Get primary state variables directly (avoids redundant Density calculations)
            bmass = rm.PVTPropertyFunctions.get(rm, state, 'BacterialMass');
            if isa(rm, 'BiochemistryPhreeqcModel')
                growthRateName = 'CarbonLimitedGrowthRate';
            else
                growthRateName = 'PsiGrowthRate';
            end
            psigrowth = model.getProps(state, growthRateName);
            % Get model parameters
            Y_H2 = bcrm.Y_H2;                      % Reaction yield scales
            gamma = rm.gammak;                     % Stoichiometric coefficients (nbioreact x ncomp)
            nbactMax = bcrm.nbactMax;              % Population scales
            mc = rm.EOSModel.CompositionalMixture.molarMass;  % Component molar masses [kg/mol]

            % Calculate conversion rate for each reactor and component
            for i = 1:nbioreact
                rxn = bcrm.metabolicReaction{i};
                isSRB = strcmp(rxn, 'SulfateReducingBacteria');

                % Find indices for this reactor's H2 and substrate
                idxH2   = find(strcmpi(compNames, bcrm.rH2{i}), 1);
                idxsub  = find(strcmpi(compNames, bcrm.rsub{i}), 1);
                idxprod = find(strcmpi(compNames, bcrm.p2{i}), 1);

                % Validate required components exist. For SRB, the
                % substrate (SO4) and one of the products (HS-) are
                % non-volatile aqueous tracers, not EOS components (see
                % BiochemistryModel.sulfateReduction /
                % SRBTracerConvRate), so they are legitimately absent
                % from compNames -- only require H2 to be present.
                if isempty(idxH2)
                    error('BactConvertionRate:MissingH2', ...
                        'Required component H2 (%s) not found for reactor %d', bcrm.rH2{i}, i);
                end
                if isempty(idxsub) && ~isSRB
                    error('BactConvertionRate:MissingSubstrate', ...
                        'Required component substrate (%s) not found for reactor %d', bcrm.rsub{i}, i);
                end

                % Get bacterial mass and growth rate for this reactor
                if iscell(bmass)
                    bmass_i = bmass{i};
                else
                    bmass_i = bmass(:, i);
                end
                if iscell(psigrowth)
                    psigrowth_i = psigrowth{i};
                else
                    psigrowth_i = psigrowth(:, i);
                end

                % Normalize stoichiometric coefficients by H2 coefficient for this reactor
                % γ_c^H2 = (γ_c * m_c) / |γ_H2|
                gamma_row = gamma(i, :);
                gamma_norm = nbactMax(i) .* gamma_row .* mc ./ abs(gamma_row(idxH2));

                % Calculate base conversion rate
                % The source is based on microbial mass and the
                % reaction-specific yield scale.
                qbase = psigrowth_i.*bmass_i./Y_H2(i);

                % For SRB, only a pH/salinity-dependent fraction of the
                % total sulfide produced is volatile H2S; the rest stays
                % as the HS- aqueous tracer (see SRBTracerConvRate,
                % which accounts for the complementary (1-f_H2S) part).
                % Splitting here, at the source term, avoids overwriting
                % flash mole fractions directly (see the note in
                % SoreideWhitsonEos.getMixtureFugacityCoefficients).
                if isSRB && ~isempty(idxprod)
                    if isfield(state, 'phreeqcPH')
                        pH = state.phreeqcPH;
                    else
                        pH = [];
                    end

                    fH2S = rm.EOSModel.fractionH2SVolatile(state.T, pH);
                end

                % Apply normalized stoichiometric coefficients to each component
                for c = 1:ncomp
                    contrib = gamma_norm(c) .* qbase;
                    if isSRB && c == idxprod
                        contrib = contrib .* fH2S;
                    end
                    qbiot{c} = qbiot{c} + contrib;
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