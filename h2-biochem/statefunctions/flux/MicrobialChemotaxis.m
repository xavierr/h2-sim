classdef MicrobialChemotaxis < StateFunction
    % Compute the effective microbial diffusivity in the liquid phase.
    %
    % SYNOPSIS:
    %   d = MicrobialChemotaxis(model)
    %   d = MicrobialChemotaxis(model, 'pn', pv, ...)
    %
    % DESCRIPTION:
    %   Evaluates the effective cell-centered microbial Chemotaxis [m^2/s]
    %
    % REQUIRED PARAMETERS:
    %   model - Model instance, must have fields `bactDiffusion` (logical)
    %           and `biochemFluid.bactdiff` (positive numeric).
    %
    % OPTIONAL PARAMETERS (property/value pairs):
    %   (none defined in this base class; can be extended via `merge_options`)
    %
    % RETURNS:
    %   dN - Cell-centered effective diffusivity (ncells x 1), units m^2/s.
    %
    % SEE ALSO:
    %   StateFunction, DynamicFlowTransmissibility

    properties
        % No additional properties needed
    end

    methods
        %-----------------------------------------------------------------%
        function md = MicrobialChemotaxis(model, varargin)
            md@StateFunction(model, varargin{:});
            md = merge_options(md, varargin{:});

            % Dependencies: pressure, saturation, and nbact (for dynamic porosity)
            md = md.dependsOn({'s','pressure','nbact'}, 'state');
            md.label = 'D_{chemot}^{eff}';
        end

        %-----------------------------------------------------------------%
        function dN = evaluateOnDomain(prop, model, state)
            bcrm=model.biochemFluid;
            nbioreact=bcrm.nbioreact;
            % Fetch primary state variables
            [ s, nbact, p] = model.getProps(state, 's','nbact', 'pressure');

            % Porosity - dynamic (bio-clogging) or static
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                nbactArray = model.extractBactValues(nbact);
                phi = model.rock.poro(p, nbactArray{:}); % Apply both modifications
            else
                phi = model.rock.poro;
            end


            % Liquid saturation
            L_ix = model.getLiquidIndex();
            if iscell(s)
                sL = s{L_ix};
            else
                sL = s(:, L_ix);
            end

            Voln = max(sL, 0);

            % Initialize with zero growth rate
            dN=cell(1,nbioreact);
            [dN{:}] = deal(0);

            % Effective microbial chemotaxis
            for i=1:nbioreact
                if iscell(nbact)
                    nbacti=nbact{i};
                else
                    nbacti=nbact(:,i);
                end
                % Effective microbial chemotaxis (volume-filling term)
                if model.chemotaxisEffect
                    if model.biochemFluid.xch_seuil(i) <= 0
                        warning(['MicrobialChemotaxis: chemotaxisEffect flag enabled but ' ...
                                 'xch_seuil coefficient is zero or negative for species %d. ' ...
                                 'Chemotaxis will be inactive.'], i);
                    end
                    populationFraction = max(0, min(nbacti ./ model.bact_maxProp, 1));
                    dN{i} = model.biochemFluid.xch_seuil(i) .* phi .* Voln .* ...
                        populationFraction .* (1.0 - populationFraction);
                else
                    dN{i} = 0;
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