classdef BiochemicalFlowDiscretization < FlowDiscretization
    % BiochemicalFlowDiscretization
    % Discretization and state function grouping for bio-chemical flow
    % within a compositional model with microbial growth.

    properties
        % No additional properties needed
    end

    methods
        %-----------------------------------------------------------------%
        function props = BiochemicalFlowDiscretization(model)
            % Constructor: inherit base FlowDiscretization properties
            props = props@FlowDiscretization(model);
            props = props.setStateFunction('ComponentTotalFlux', ...
                ComponentTotalFluxForBio(model));

            % Add dispersion state functions
            if (isprop(model, 'molecularDispersion') && model.molecularDispersion) || ...
                    (isprop(model, 'molecularDiffusion') && model.molecularDiffusion)
                props = props.setStateFunction('DispersiveDiffusivity', DispersiveDiffusivity(model));
                props = props.setStateFunction('DispersiveTransmissibility', ...
                    DispersiveTransmissibility(model));
                props = props.setStateFunction('ComponentPhaseDispFlux', ...
                    ComponentPhaseDispFlux(model));
                props = props.setStateFunction('ComponentTotalDispFlux', ...
                    ComponentTotalDispFlux(model));
            end

            if model.bacteriamodel
                % Ensure Porosity exists in FlowDiscretization for diffusion (also when no clogging)
                props = props.setStateFunction('Porosity', PorosityFromRock(model));

                % Set up transmissibility and porosity functions
                if model.dynamicFlowTrans
                    % Dynamic transmissibility computed each nonlinear iteration
                    props = props.setStateFunction('Permeability', BactPermeability(model));
                    props = props.setStateFunction('Transmissibility', ...
                        DynamicFlowTransmissibility(model, 'Permeability'));
                end

                if model.dynamicFlowPv
                    % Dynamic porosity / pore volume
                    props = props.setStateFunction('Porosity', BactPorosity(model));
                    props = props.setStateFunction('PoreVolume', ...
                        DynamicFlowPoreVolume(model, 'Porosity'));
                end

                if model.bactDiffusion
                    props = props.setStateFunction('MicrobialDiffusivity', MicrobialDiffusivity(model));
                    props = props.setStateFunction('MicrobialTransmissibility', ...
                        MicrobialTransmissibility(model));
                    props = props.setStateFunction('BactFlux', DiffusiveBactFlux(model));
                end
                 if model.chemotaxisEffect
                    props = props.setStateFunction('MicrobialChemotaxis', MicrobialChemotaxis(model));
                    props = props.setStateFunction('ChemotaxisTransmissibility', ...
                        ChemotaxisTransmissibility(model));
                    props = props.setStateFunction('ChemoBactFlux', ChemotaxisBactFlux(model));
                end

                if ismethod(model, 'hasMobileAqueousTracers')
                    hasTracers = model.hasMobileAqueousTracers();
                else
                    hasTracers = isprop(model, 'sulfateReduction') && model.sulfateReduction;
                end
                if hasTracers
                    if model.sulfateReduction
                    props = props.setStateFunction('SRBTracerConvRate', SRBTracerConvRate(model));
                    end

                    if (isprop(model, 'molecularDispersion') && model.molecularDispersion) || ...
                            (isprop(model, 'molecularDiffusion') && model.molecularDiffusion)
                        props = props.setStateFunction('TracerDiffusivity', TracerDiffusivity(model));
                        props = props.setStateFunction('TracerTransmissibility', ...
                            TracerTransmissibility(model));
                        props = props.setStateFunction('DiffusiveTracerFlux', DiffusiveTracerFlux(model));
                    end
                end
            end
        end

        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            % Call parent method for standard component conservation
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
        end

        %-----------------------------------------------------------------%
        function [acc, bflux, name, type] = bacteriaConservationEquation(fd, model, state, state0, dt)
            % Computes bacterial mass conservation equation
            bactmass  = model.getProp(state, 'BacterialMass');
            bactmass0 = model.getProp(state0, 'BacterialMass');

            % Accumulation term: d(bactmass)/dt
            nbioreact=model.biochemFluid.nbioreact;
            acc=cell(1,nbioreact);
            [acc{:}] = deal(0);
            for i=1:nbioreact
                acc{i} = (bactmass{i} - bactmass0{i}) ./ dt;
            end
            % Add microbial diffusion contributions if present
             bflux =cell(1,nbioreact);
            [bflux{:}] = deal(0);
            if (model.bactDiffusion && ~model.chemotaxisEffect)
                flowState = fd.buildFlowState(model, state, state0, dt);
                bflux = model.getProp(flowState, 'BactFlux');
            elseif (~model.bactDiffusion && model.chemotaxisEffect)
                flowState = fd.buildFlowState(model, state, state0, dt);
                bflux = model.getProp(flowState, 'ChemoBactFlux');
            elseif (model.bactDiffusion && model.chemotaxisEffect)
                flowState = fd.buildFlowState(model, state, state0, dt);
                BactFlux = model.getProp(flowState, 'BactFlux');
                ChemoFlux = model.getProp(flowState, 'ChemoBactFlux');
                bflux = cell(size(BactFlux));
                for i = 1:numel(BactFlux)
                    bflux{i} = BactFlux{i} + ChemoFlux{i};
                end
            end

            % Output variable names and types
            name = model.biochemFluid.bactnames;
            type=cell(1,nbioreact);
            [type{:}] = deal('cell');
        end

        %-----------------------------------------------------------------%
        function [acc, tflux, name, type] = tracerConservationEquation(fd, model, state, state0, dt)
            % Mass conservation for mobile non-EOS aqueous tracers.
            %
            % Unlike EOS components, these tracers never enter the flash:
            % they are advected with the liquid phase, using the same
            % interior Darcy flux/upstream weighting as any other
            % dissolved species (see PhaseFlux, PhaseUpwindFlag), with no
            % vapor-phase contribution. When molecularDiffusion and/or
            % molecularDispersion are enabled, a Fickian flux (see
            % TracerDiffusivity/DiffusiveTracerFlux) is added on top of
            % the advective flux, analogous to how those flags affect the
            % volatile EOS components via ComponentTotalFluxForBio.
            tracermass  = model.getProp(state, 'AqueousTracerMass');
            tracermass0 = model.getProp(state0, 'AqueousTracerMass');
            if ismethod(model, 'getAqueousTracerNames')
                name = model.getAqueousTracerNames();
            else
                name = {'SO4', 'HS'};
            end
            ntracer = numel(name);

            acc = cell(1, ntracer);
            for i = 1:ntracer
                acc{i} = (tracermass{i} - tracermass0{i}) ./ dt;
            end

            q    = model.getProp(state, 'PhaseFlux');
            flag = model.getProp(state, 'PhaseUpwindFlag');
            L_ix = model.getLiquidIndex();
            qL    = q{L_ix};
            flagL = flag{L_ix};

            c = cell(1, ntracer);
            for i = 1:ntracer
                c{i} = model.getProp(state, lower(name{i}));
            end

            tflux = cell(1, ntracer);
            for i = 1:ntracer
                tflux{i} = model.operators.faceUpstr(flagL, c{i}) .* qL;
            end

            if model.molecularDispersion || model.molecularDiffusion
                Jdiff = model.getProp(state, 'DiffusiveTracerFlux');
                for i = 1:ntracer
                    tflux{i} = tflux{i} + Jdiff{i};
                end
            end

            type = repmat({'cell'}, 1, ntracer);
        end

        %-----------------------------------------------------------------%
        function c = getPopulationMass(model, state, extra)
            % Compute microbial mass in each phase
            c = component.getPhaseComposition(model, state);

            if nargin < 4
                rho = model.getProp(state, 'Density');
            else
                rho = extra.rho;
            end

            for ph = 1:numel(c)
                if ~isempty(c{ph})
                    c{ph} = rho{ph} .* c{ph};
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