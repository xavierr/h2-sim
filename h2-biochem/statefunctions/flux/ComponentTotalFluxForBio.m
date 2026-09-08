classdef ComponentTotalFluxForBio < ComponentTotalFlux
    % Total component flux for bio‑chemical models.
    % Inherits from ComponentTotalFlux and optionally adds molecular
    % dispersion and diffusion fluxes when the corresponding model flags
    % are enabled:
    %   model.molecularDispersion  (logical)
    %   model.molecularDiffusion   (logical)
    %
    % SEE ALSO: ComponentTotalFlux, ComponentTotalMolecularDispFlux

    properties (Access = protected)
        hasDispersion_  % Cache flag for performance
        hasDiffusion_   % Cache flag for performance
    end

    methods
        function sf = ComponentTotalFluxForBio(model, varargin)
            % Construct a ComponentTotalFluxForBio object.
            %
            % PARAMETERS:
            %   model - MRST model with optional flags:
            %           .molecularDispersion (logical)
            %           .molecularDiffusion  (logical)
            %   varargin - Additional arguments passed to parent constructor.
            %
            % RETURNS:
            %   sf - ComponentTotalFluxForBio instance.

            % Call parent constructor first
            sf = sf@ComponentTotalFlux(model, varargin{:});

            % Cache model flags for performance
            sf.hasDispersion_ = isprop(model, 'molecularDispersion') && ...
                                model.molecularDispersion;

            sf.hasDiffusion_ = isprop(model, 'molecularDiffusion') && ...
                               model.molecularDiffusion;

            % Add dependencies for the extra physics
            if sf.hasDispersion_ || sf.hasDiffusion_
                sf = sf.dependsOn('ComponentTotalDispFlux');
            end

            % Set the label for this flux (used in equation display)
            sf.label = 'V_i';
        end

        function v = evaluateOnDomain(sf, model, state)
            % Evaluate the total component flux on the domain.
            % 1. Get the standard advective total flux (cell array per component)
            v = evaluateOnDomain@ComponentTotalFlux(sf, model, state);

            % 2. Add molecular diffusion/dispersion flux if enabled
            if sf.hasDispersion_ || sf.hasDiffusion_
                % Get the disp flux from dependencies
                Jdiff = sf.getEvaluatedDependencies(state, 'ComponentTotalDispFlux');

                % Ensure we have a cell array and add to each component
                if iscell(Jdiff)
                    ncomp = min(numel(v), numel(Jdiff));
                    for c = 1:ncomp
                        if ~isempty(Jdiff{c})
                            v{c} = v{c} + Jdiff{c};
                        end
                    end
                elseif ~isempty(Jdiff)
                    % Single flux value (same for all components)
                    for c = 1:numel(v)
                        v{c} = v{c} + Jdiff;
                    end
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