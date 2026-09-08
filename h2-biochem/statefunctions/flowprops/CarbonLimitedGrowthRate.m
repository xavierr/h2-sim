classdef CarbonLimitedGrowthRate < StateFunction
    % Growth rate with a finite aqueous inorganic-carbon constraint.
    %
    % For carbonate-buffered runs, the kinetic CO2 proxy used by
    % GrowthBactRateSRC can exceed the carbon actually present in a cell.
    % This property limits methanogenic and acetogenic growth so their
    % combined CO2 demand over the implicit timestep cannot exceed the
    % aqueous EOS CO2 plus HCO3- inventory at the timestep start.

    methods
        function clgr = CarbonLimitedGrowthRate(model, varargin)
            clgr@StateFunction(model, varargin{:});
            clgr = clgr.dependsOn('PsiGrowthRate', 'state');

            if isprop(model, 'ReservoirModel') && ~isempty(model.ReservoirModel)
                rm = model.ReservoirModel;
            else
                rm = model;
            end
            if isprop(rm, 'carbonateBuffer') && rm.carbonateBuffer
                clgr = clgr.dependsOn('tracerHCO3', 'state');
                clgr = clgr.dependsOn('carbonSubstrateDt', 'state');
                clgr = clgr.dependsOn('carbonSubstrateMoles', 'state');
            end
            clgr.label = '\Psi_{growth,carbon-limited}';
        end

        function psigrowth = evaluateOnDomain(prop, model, state)
            if isprop(model, 'ReservoirModel') && ~isempty(model.ReservoirModel)
                rm = model.ReservoirModel;
            else
                rm = model;
            end

            psigrowth = model.getProps(state, 'PsiGrowthRate');
            if ~(isprop(rm, 'carbonateBuffer') && rm.carbonateBuffer)
                return;
            end

            assert(isfield(state, 'carbonSubstrateDt') && isfield(state, 'carbonSubstrateMoles') && ...
                isscalar(state.carbonSubstrateDt) && state.carbonSubstrateDt > 0, ...
                ['CarbonLimitedGrowthRate requires the positive current timestep and ', ...
                 'the start-of-step aqueous carbon inventory.']);

            bcrm = rm.biochemFluid;
            cnames = rm.EOSModel.getComponentNames();
            idxCO2 = find(strcmpi(cnames, 'CO2'), 1);
            assert(~isempty(idxCO2), ...
                'carbonateBuffer requires an EOS CO2 component.');

            L_ix = rm.getLiquidIndex();
            % [mol/s]: finite liquid CO2 plus HCO3- inventory at the
            % timestep start, divided by the current timestep duration.
            carbonAvailable = state.carbonSubstrateMoles./state.carbonSubstrateDt;

            bmass = rm.PVTPropertyFunctions.get(rm, state, 'BacterialMass');
            rho = rm.PVTPropertyFunctions.get(rm, state, 'Density');
            if iscell(rho)
                rhoL = rho{L_ix};
            else
                rhoL = rho(:, L_ix);
            end

            gamma = rm.gammak;
            molarMass = rm.EOSModel.CompositionalMixture.molarMass;
            co2Demand = 0;
            carbonReactions = false(1, bcrm.nbioreact);
            for i = 1:bcrm.nbioreact
                reaction = bcrm.metabolicReaction(i);
                carbonReactions(i) = strcmpi(reaction, 'MethanogenicArchae') || ...
                    strcmpi(reaction, 'AcetogenicBacteria');
                if ~carbonReactions(i)
                    continue;
                end

                assert(strcmpi(bcrm.rsub(i), 'CO2'), ...
                    'Carbon-limited growth requires CO2 as the MET/ACE substrate.');
                idxH2 = find(strcmpi(cnames, bcrm.rH2(i)), 1);
                assert(~isempty(idxH2) && gamma(i, idxH2) < 0 && ...
                    gamma(i, idxCO2) < 0, ...
                    'Carbon-limited growth requires consuming H2 and CO2 stoichiometries.');

                if iscell(psigrowth)
                    psigrowth_i = psigrowth{i};
                else
                    psigrowth_i = psigrowth(:, i);
                end
                if iscell(bmass)
                    bmass_i = bmass{i};
                else
                    bmass_i = bmass(:, i);
                end

                qbase = psigrowth_i.*bmass_i./bcrm.Y_H2(i);
                gammaCO2 = bcrm.nbactMax(i).*gamma(i, idxCO2).* ...
                    molarMass(idxCO2)./abs(gamma(i, idxH2));

                % This is exactly the CO2 molar sink represented by
                % BactConvertionRate after BiochemistryModel divides its
                % mass source by rhoL.
                co2Demand = co2Demand + max( ...
                    -gammaCO2.*qbase./rhoL./molarMass(idxCO2), 0);
            end

            % Smoothly approximate min(1, available/demand). A hard min
            % places the implicit reaction solve exactly on the
            % zero-carbon boundary and creates a nonsmooth Jacobian that
            % triggers severe timestep cutting. The p-norm form remains
            % conservative while retaining useful derivatives.
            carbonRatio = carbonAvailable./max(co2Demand, 1e-30);
            smoothnessOrder = 8;
            carbonScale = carbonRatio./ ...
                (1 + carbonRatio.^smoothnessOrder).^(1/smoothnessOrder);
            for i = find(carbonReactions)
                if iscell(psigrowth)
                    psigrowth{i} = psigrowth{i}.*carbonScale;
                else
                    psigrowth(:, i) = psigrowth(:, i).*carbonScale;
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
