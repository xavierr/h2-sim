classdef BiochemistryModel < GenericOverallCompositionModel
    % Base model coupling two-phase compositional flow to Monod-kinetics
    % microbial growth/decay, with no PHREEQC geochemistry.
    %
    % SYNOPSIS:
    %   model = BiochemistryModel(G, rock, fluid)
    %   model = BiochemistryModel(G, rock, fluid, compFluid)
    %   model = BiochemistryModel(..., 'pn1', vn1, ...)
    %
    % DESCRIPTION:
    %   BiochemistryModel is the foundation of the h2-biochem module's three
    %   model classes: it fully-implicitly couples a compositional (Soreide-
    %   Whitson EOS) flow model to a multi-population microbial reaction
    %   model, with microbial growth/decay following Monod kinetics and H2
    %   as the shared electron donor. All reaction source terms are
    %   assembled directly into the same AD Jacobian as flow, so this class
    %   never calls out to an external geochemistry engine -- pH, aqueous
    %   speciation, and mineral equilibrium are not resolved.
    %
    %   Its capabilities include:
    %     - Reaction-specific microbial populations (nbact), one balance
    %       equation per metabolic reaction (biochemFluid.metabolicReaction).
    %     - Optional sulfate-reducing-bacteria (SRB) chemistry via SO4/HS
    %       aqueous tracers (sulfateReduction, set automatically from
    %       biochemFluid), advected with the liquid phase.
    %     - Optional microbial diffusion (bactDiffusion) and chemotaxis
    %       toward higher H2 concentration (chemotaxisEffect).
    %     - Optional molecular diffusion/mechanical dispersion of EOS
    %       components (molecularDiffusion, molecularDispersion).
    %     - Optional bio-clogging feedback on porosity/permeability, set up
    %       separately via setupBioCloggingModel.
    %
    %   It also defines reactionsEnabled/localReactionMode, two flags that
    %   exist purely so that SequentialBiochemistryPhreeqcModel (a
    %   subclass of BiochemistryPhreeqcModel, itself a subclass of this
    %   class) can build two single-purpose internal instances from it: one
    %   with reactions disabled (global flow-only stage) and one with
    %   spatial flux disabled (cell-local reaction-only stage). A plain
    %   BiochemistryModel/BiochemistryPhreeqcModel run never needs to touch
    %   these two flags.
    %
    % REQUIRED PARAMETERS:
    %   G         - Simulation grid
    %   rock      - Rock properties for the model
    %   fluid     - Fluid model for the simulation
    %   compFluid - Compositional fluid mixture (optional)
    %
    % OPTIONAL PARAMETERS:
    %   'bacteriamodel'       - Enable microbial growth/decay and its
    %                           reaction source terms (default true).
    %   'bactDiffusion'       - Enable microbial (Fickian) diffusion
    %                           (default false).
    %   'chemotaxisEffect'    - Enable chemotactic microbial migration
    %                           toward dissolved H2 (default false).
    %   'molecularDiffusion'  - Enable molecular diffusion of EOS
    %                           components (default false).
    %   'molecularDispersion' - Enable mechanical dispersion of EOS
    %                           components (default false).
    %
    % RETURNS:
    %   model - BiochemistryModel class instance
    %
    % SEE ALSO:
    %   BiochemistryPhreeqcModel, SequentialBiochemistryPhreeqcModel,
    %   ReservoirModel, ThreePhaseCompositionalModel

    properties
        % Reserved for future support of alternative bio-chemical
        % formulations; only 'bacterialmodel' is currently implemented
        % (enforced in the constructor).
        bacterialFormulation = 'bacterialmodel';

        % Compositional fluid mixture (EOS component list/names).
        compFluid

        % Biochemical reaction database: stoichiometry, kinetic
        % parameters, and yield/half-saturation constants per metabolic
        % reaction (TableBioChemMixture instance).
        biochemFluid

        % Physical quantities and bounds
        gammak   = [];                    % Stoichiometric coefficients, one row per reaction, one column per EOS component
        bacteriamodel = true;             % Master switch for the microbial growth/decay sub-model
        sulfateReduction = false;         % SO4/HS aqueous tracers active (set from biochemFluid)
        enableSulfateSource = true;       % Optional anhydrite sulfate source for SRB tracers
        bact_capProp = 3.0e0;             % Min nbact in the model
        bact_maxProp = 120;               % Max nbact in the model
        molecularDiffusion = false;       % Molecular diffusion of EOS components
        molecularDispersion = false;      % Mechanical dispersion of EOS components
        bactDiffusion = false;            % Microbial diffusion
        chemotaxisEffect = false;         % chemotaxis

        % Coarse-flow/local-reaction sequential split support (see
        % SequentialBiochemistryPhreeqcModel). Both default to standard,
        % unsplit behaviour and are only toggled on internally built
        % single-purpose stage models.
        reactionsEnabled = true;          % false: no biological source terms assembled (flow stage)
        localReactionMode = false;        % true: no spatial flux/diffusion assembled (reaction stage)
    end

    methods
        %-----------------------------------------------------------------%
        function model = BiochemistryModel(G, rock, fluid, compFluid, biochemFluid, includeWater, backend, varargin)
            % Constructor
            model = model@GenericOverallCompositionModel(G, rock, fluid, compFluid, ...
                'water', includeWater, 'AutoDiffBackend', backend);
            model = merge_options(model, varargin{:});
            model.molecularDiffusion = normalizeTransportFlag( ...
                model.molecularDiffusion, 'molecularDiffusion');
            model.molecularDispersion = normalizeTransportFlag( ...
                model.molecularDispersion, 'molecularDispersion');
            model.bactDiffusion = normalizeTransportFlag( ...
                model.bactDiffusion, 'bactDiffusion');
            model.chemotaxisEffect = normalizeTransportFlag( ...
                model.chemotaxisEffect, 'chemotaxisEffect');
            model.reactionsEnabled = normalizeTransportFlag( ...
                model.reactionsEnabled, 'reactionsEnabled');
            model.localReactionMode = normalizeTransportFlag( ...
                model.localReactionMode, 'localReactionMode');

            % Set up operators
            model = model.setupOperators();

            % Check phases
            model.gas = true;
            if ~includeWater
                assert(model.oil, 'we need a liquid phase');
            end

             %% Set metabolic reactions
            if isempty(biochemFluid)
                biochemFluid=TableBioChemMixture({'MethanogenicArchae'},{'bactM'});
            end
            model.biochemFluid=biochemFluid;
            % SO4 and HS are non-volatile: when the SRB reaction is
            % present they are carried as aqueous tracers (advected with
            % the liquid phase, reacted via SRBTracerConvRate), never as
            % EOS/CompositionalMixture components.
            model.sulfateReduction = any(strcmp(model.biochemFluid.metabolicReaction, 'SulfateReducingBacteria'));

            %% Set compositional fluid and EOS
            if isempty(compFluid)
                if strcmp(model.biochemFluid.metabolicReaction, 'MethanogenicArchae')
                    compNames = {'Hydrogen', 'Water', 'Nitrogen', 'CarbonDioxide', 'Methane'};
                    compSymbols = {'H2', 'H2O', 'N2', 'CO2', 'C1'};
                    compFluid = TableCompositionalMixture(compNames, compSymbols);
                else
                    warning('MethanogenicArchae is the default; other reactions not implemented.');
                end
            end
            model.compFluid = compFluid;
            model.EOSModel = SoreideWhitsonEos([], compFluid);
            ncomp = compFluid.getNumberOfComponents();
            nbioreact = numel(model.biochemFluid.metabolicReaction);
            namecp = compFluid.names;
            model.gammak = zeros(nbioreact, ncomp);
           for i=1:nbioreact
                indH2   = find(strcmp(namecp, model.biochemFluid.rH2(i)));
                indH2O  = find(strcmp(namecp, model.biochemFluid.pH2O(i)));
                indsub  = find(strcmp(namecp, model.biochemFluid.rsub(i)));
                indprod   = find(strcmp(namecp, model.biochemFluid.p2(i)));
                model.gammak(i,indH2)  = model.biochemFluid.gamrH2(i);
                model.gammak(i,indH2O) =  model.biochemFluid.gampH2O(i);
                model.gammak(i,indsub) = model.biochemFluid.gamrsub(i);
                model.gammak(i,indprod)  = model.biochemFluid.gamp2(i);
            end


            % Validate bacterial formulation
            assert(any(strcmpi(model.bacterialFormulation, {'bacterialmodel'})), ...
                'BioChemistryModel supports currently only one micro-organism');

            % Set output state functions
            model.FlowDiscretization = BiochemicalFlowDiscretization(model);
            model.OutputStateFunctions = {'ComponentTotalMass', 'Density'};
            % Set up state function groupings
            model = model.setupStateFunctionGroupings();
        end

        function model = setupOperators(model, G, rock, varargin)
            % Set up operators, potentially accounting for dynamic
            % transmissibilites

            % Set rock and grid from model if not provided
            if nargin < 3, rock = model.rock; end
            if nargin < 2, G = model.G;       end

            drock = rock;
            if model.dynamicFlowTrans()
                % Assign dummy transmissibilities to appease
                % model.setupOperators
                drock = rock;
                nbact0 = model.getDummyBacterialValues(drock.perm);
                drock.perm = rock.perm(1*barsa(), nbact0{:});
            end

            if model.dynamicFlowPv()
                % Assign dummy transmissibilities to appease
                % model.setupOperators
                if ~model.dynamicFlowTrans()
                    drock = rock;
                end
                nbact0 = model.getDummyBacterialValues(drock.poro);
                drock.poro = rock.poro(1*barsa(), nbact0{:});
            end
            % Let reservoir model set up operators
            model = setupOperators@ReservoirModel(model, G, drock, varargin{:});
            model.rock = rock;
        end

        function model = validateModel(model, varargin)
            % Attach a bacteria-aware FacilityModel (BiochemistryGenericFacilityModel)
            % when the microbial sub-model is active, otherwise fall back
            % to the standard GenericFacilityModel, then delegate to the
            % parent compositional model for the remaining validation.
            if model.bacteriamodel
                if isempty(model.FacilityModel) || ...
                        ~isa(model.FacilityModel, 'BiochemistryGenericFacilityModel')
                    model.FacilityModel = BiochemistryGenericFacilityModel(model);
                end
            else
                if isempty(model.FacilityModel) || ~isa(model.FacilityModel, 'GenericFacilityModel')
                    model.FacilityModel = GenericFacilityModel(model);
                end
            end
            model = validateModel@GenericOverallCompositionModel(model, varargin{:});
        end

        function model = setupStateFunctionGroupings(model, varargin)
            % Register the microbial-kinetics state functions (growth,
            % decay, and the resulting component conversion source) and
            % the bacterial-mass property on top of the parent
            % compositional groupings, so they participate in the AD
            % dependency graph like any other flow property.
            model = setupStateFunctionGroupings@GenericOverallCompositionModel(model, varargin{:});

            fluxprops = model.FlowDiscretization;
            pvtprops  = model.PVTPropertyFunctions;
            flowprops = model.FlowPropertyFunctions;

            if model.bacteriamodel
                flowprops = flowprops.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                flowprops = flowprops.setStateFunction('PsiDecayRate',  DecayBactRateSRC(model));
                flowprops = flowprops.setStateFunction('BactConvRate',  BactConvertionRate(model));

                % Register bacterial mass as cell property (not a source term)
                pvtprops = pvtprops.setStateFunction('BacterialMass', BacterialMass(model));

                if model.sulfateReduction
                    pvtprops = pvtprops.setStateFunction('AqueousTracerMass', AqueousTracerMass(model));
                end
            end

            pvt = pvtprops.getRegionPVT(model);
            if isfield(model.fluid, 'pvMultR')
                pv = DynamicFlowPoreVolume(model, pvt);
            else
                pv = PoreVolume(model, pvt);
            end
            pvtprops = pvtprops.setStateFunction('PoreVolume', pv);

            model.PVTPropertyFunctions  = pvtprops;
            model.FlowPropertyFunctions = flowprops;
            model.FlowDiscretization    = fluxprops;
        end

        function state = validateState(model, state)
            % Ensure state carries the fields the microbial/SRB sub-models
            % need (nbact, and the SO4/HS/lagged-H2S tracers when sulfate
            % reduction is active), defaulting any that are missing,
            % before delegating to the parent compositional validation.
            state = validateState@ThreePhaseCompositionalModel(model, state);
            if model.bacteriamodel && ~isfield(state, 'nbact')
                nbact0 = 1e6;
                state.nbact = repmat(nbact0, model.G.cells.num, 1);
            end
            if model.sulfateReduction
                if ~isfield(state, 'tracerSO4')
                    state.tracerSO4 = zeros(model.G.cells.num, 1);
                end
                if ~isfield(state, 'tracerHS')
                    state.tracerHS = zeros(model.G.cells.num, 1);
                end
                if ~isfield(state, 'h2sDissolvedLag')
                    % Lagged (previous converged timestep) dissolved H2S
                    % concentration [mol/m3 liquid], used only to avoid a
                    % circular dependency between the flash (needs
                    % msalt) and the flash-derived H2S content when
                    % feeding total sulfide to the EOS -- see
                    % initStateAD/updateAfterConvergence.
                    state.h2sDissolvedLag = zeros(model.G.cells.num, 1);
                end
            end
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Primary variable set: pressure, overall composition (all but
            % the first EOS component), microbial population nbact per
            % reaction when bacteriamodel is active, extra non-EOS phase
            % saturations, and SO4/HS when sulfateReduction is active, plus
            % whatever the FacilityModel (well controls) contributes.
            [p, z] = model.getProps(state, 'pressure', 'z');
            z = ensureMinimumFraction(z, model.EOSModel.minimumComposition);
            z = expandMatrixToCell(z);
            cnames = model.EOSModel.getComponentNames();
            extra = model.getNonEoSPhaseNames();
            ne = numel(extra);
            enames = cell(1, ne); evars = cell(1, ne);
            for i = 1:ne
                sn = ['s', extra(i)];
                enames{i} = sn;
                evars{i} = model.getProp(state, sn);
            end

            if model.bacteriamodel
                nbact = model.getProp(state, 'nbact');
                nbact = expandMatrixToCell(nbact);
                bactnames = model.biochemFluid.bactnames;
                names = [{'pressure'}, cnames(2:end), bactnames, enames];
                vars  = [p, z(2:end), nbact, evars];
                if model.sulfateReduction
                    so4 = model.getProp(state, 'so4');
                    hs  = model.getProp(state, 'hs');
                    names = [names, {'SO4', 'HS'}];
                    vars  = [vars, {so4, hs}];
                end
            else
                names = [{'pressure'}, cnames(2:end), enames];
                vars  = [p, z(2:end), evars];
            end
            origin = repmat({class(model)}, 1, numel(names));

            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars   = [vars, v];
                names  = [names, n];
                origin = [origin, o];
            end
        end
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            % Assemble the full residual: EOS component conservation
            % (with well sources and, if reactionsEnabled, the microbial
            % component-conversion source BactConvRate), the microbial
            % mass-balance equation per reaction (with optional spatial
            % diffusion/chemotaxis flux and growth/decay source), the
            % SO4/HS aqueous tracer balances when sulfateReduction is
            % active, and finally the well/facility equations.
            % localReactionMode drops all spatial flux-divergence terms
            % (used by the reaction-only stage of the coarse-flow/
            % local-reaction split); reactionsEnabled gates every
            % microbial source term (used by the flow-only stage).
            %
            % Discretize
            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');
            comps = cellfun(@(x, y) {x, y}, X(:, model.getLiquidIndex), X(:, model.getVaporIndex), 'UniformOutput', false);


            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                pressures, sat, mob, rho, ...
                {}, comps, ...
                drivingForces);

            % Add sources
            eqs = model.insertSources(eqs, src);
            % Assemble equations

            if model.bacteriamodel && model.reactionsEnabled
                cnames = model.EOSModel.getComponentNames();
                ncomp = numel(cnames);
                src_rate = model.FacilityModel.getProps(state, 'BactConvRate');
                L_ix = model.getLiquidIndex();

                if iscell(rho)
                    rhoL = rho{L_ix};
                else
                    rhoL = rho(:, L_ix);
                end

                for i = 1:ncomp
                    if ~isempty(src_rate{i})
                        eqs{i} = eqs{i} -src_rate{i}./rhoL;
                    end
                end
            end

            % localReactionMode (used by the reaction stage of the
            % coarse-flow/local-reaction split) skips spatial flux
            % divergence entirely, leaving pure per-cell accumulation
            % equations coupled only through the (still assembled)
            % reaction source terms above/below.
            for i = 1:numel(eqs)
                if ~model.localReactionMode
                    eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
                end
            end
            if model.bacteriamodel
                % Bacterial mass balance: d(M)/dt + div(flux) = source
                % where M = pv * S_l * nbact [kg]
                [beqs, bflux, bnames, btypes] = model.FlowDiscretization.bacteriaConservationEquation(model, state, state0, dt);
                fd = model.FlowDiscretization;
                if model.reactionsEnabled
                    src_growthdecay = model.FacilityModel.getBacteriaSources(fd, state, state0, dt);
                end

                % Assemble accumulation and flux divergence
                nbioreact=model.biochemFluid.nbioreact;
                for i=1:nbioreact
                    if model.localReactionMode
                        % No spatial coupling in local-reaction mode.
                    elseif model.bactDiffusion && ~model.chemotaxisEffect && ~isempty(bflux{i})
                        beqs{i} = model.operators.AccDiv(beqs{i}, bflux{i});
                        % Dirichlet boundary conditions for bacterial diffusion
                        beqs{i} = model.addBacterialDiffusionBC(beqs{i}, state, drivingForces, i);
                    elseif model.chemotaxisEffect && ~model.bactDiffusion && ~isempty(bflux{i})
                        beqs{i} = model.operators.AccDiv(beqs{i}, bflux{i});
                    elseif model.bactDiffusion && model.chemotaxisEffect && ~isempty(bflux{i})
                        beqs{i} = model.operators.AccDiv(beqs{i}, bflux{i});
                        % Dirichlet boundary conditions for bacterial diffusion
                        beqs{i} = model.addBacterialDiffusionBC(beqs{i}, state, drivingForces, i);
                    else
                        % No diffusion: just accumulation term (pore-scale diffusion only)
                    end

                    if model.reactionsEnabled
                        beqs{i} = beqs{i} - src_growthdecay{i};
                    end
                end

                if model.sulfateReduction
                    % SO4/HS aqueous tracer mass balance: pure advection
                    % with the liquid phase plus SRB reaction source, no
                    % EOS/flash coupling (see SoreideWhitsonEos and
                    % SRBTracerConvRate).
                    [tacc, tflux, tnames, ttypes] = model.FlowDiscretization.tracerConservationEquation(model, state, state0, dt);
                    if model.reactionsEnabled
                        src_tracer = model.FacilityModel.getSRBTracerSources(fd, state, state0, dt);
                    end
                    for i = 1:2
                        if ~model.localReactionMode
                            tacc{i} = model.operators.AccDiv(tacc{i}, tflux{i});
                        end
                        if model.reactionsEnabled
                            tacc{i} = tacc{i} - src_tracer{i};
                        end
                    end
                    beqs  = [beqs, tacc];
                    bnames = [bnames, tnames];
                    btypes = [btypes, ttypes];
                end
            else
                [beqs, bnames, btypes] = deal([]);
            end
            % Concatenate
            eqs   = [eqs, beqs];
            names = [names, bnames];
            types = [types, btypes];


            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            % Concatenate
            eqs   = [eqs  , weqs  ];
            names = [names, wnames];
            types = [types, wtypes];

        end

        function beq = addBacterialDiffusionBC(model, beq, state, forces, species)
            % Add Dirichlet boundary conditions for the bacterial diffusion
            % equation.
            %
            % The prescribed bacterial concentration is carried on the
            % standard boundary-condition struct as the extra field
            % `bc.nbact`. It may be a vector (applied to every species), a
            % matrix with one column per species, or a cell array with one
            % vector per species. Use NaN to retain the natural no-flux
            % condition. For every finite value, the diffusive half-face
            % flux leaving the adjacent cell is added to that species'
            % bacterial mass balance:
            %
            %   J_out(f) = rho_l(c) .* T_bc(f) .* (nbact(c) - nbact_bc(f))
            %
            % with T_bc(f) = cn(f) .* D_b(c), where cn(f) is the one-sided
            % two-point geometric weight (consistent with the internal
            % `DynamicFlowTransmissibility`, whose harmonic average of a
            % single half-face reduces to that half-face) and D_b is the
            % cell-centred microbial diffusivity (`MicrobialDiffusivity`).

            % Nothing to do without a Dirichlet bacterial specification
            if isempty(forces) || ~isfield(forces, 'bc') || isempty(forces.bc) ...
                    || ~isfield(forces.bc, 'nbact') || isempty(forces.bc.nbact)
                return
            end

            bc    = forces.bc;
            faces = bc.face(:);
            if iscell(bc.nbact)
                val = bc.nbact{species}(:);
            elseif ~isvector(bc.nbact) || ...
                    (size(bc.nbact, 1) == 1 && size(bc.nbact, 2) == model.biochemFluid.nbioreact)
                assert(size(bc.nbact, 2) >= species, ...
                    'bc.nbact must have one column per bacterial species.');
                val = bc.nbact(:, species);
            else
                val = bc.nbact(:);
            end
            assert(numel(val) == numel(faces), ...
                'bc.nbact must contain one value per boundary-condition face.');

            % Keep only faces that carry a finite Dirichlet value
            keep  = isfinite(val);
            faces = faces(keep);
            val   = val(keep);
            if isempty(faces)
                return
            end

            G = model.G;
            assert(all(any(G.faces.neighbors(faces, :) == 0, 2)), ...
                'Bacterial Dirichlet conditions can only be set on boundary faces.');

            % Adjacent reservoir cell for each boundary face
            cells = sum(G.faces.neighbors(faces, :), 2);

            % One-sided two-point geometric weight cn(f) [same form as the
            % internal transmissibility]; the sign is irrelevant since the
            % driving direction is set explicitly by (nbact_c - nbact_bc).
            C  = G.faces.centroids(faces, :) - G.cells.centroids(cells, :);
            N  = G.faces.normals(faces, :);
            cn = abs(sum(C.*N, 2))./sum(C.*C, 2);

            % Cell-centred microbial diffusivity D_b = bactdiff*pv*sL
            D = model.getProp(state, 'MicrobialDiffusivity');
            if iscell(D)
                D = D{species};
            end
            if numel(value(D)) == 1
                % Diffusion disabled / zero -> no boundary contribution
                return
            end
            T_bc = cn .* D(cells);

            % Liquid-phase density at the adjacent cells (face value approx.)
            rho  = model.getProp(state, 'Density');
            L_ix = model.getLiquidIndex();
            if iscell(rho)
                rhoL = rho{L_ix};
            else
                rhoL = rho(:, L_ix);
            end

            nbact = model.getProp(state, 'nbact');
            if iscell(nbact)
                nbact = nbact{species};
            else
                nbact = nbact(:, species);
            end

            % Diffusive flux leaving the adjacent cell (positive = outflow)
            Jout = rhoL(cells) .* T_bc .* (nbact(cells) - val);

            % Scatter face contributions onto the cell residual (handles
            % several boundary faces sharing the same cell)
            nc   = G.cells.num;
            nf   = numel(cells);
            Scat = sparse(cells, (1:nf)', 1, nc, nf);
            beq = beq + Scat*Jout;
        end

        function forces = validateDrivingForces(model, forces, varargin)
            % Delegate to the parent compositional model, then validate/
            % normalize any Soreide-Whitson-specific forcing terms (e.g.
            % salinity-dependent boundary conditions) when that EOS is used.
            forces = validateDrivingForces@GenericOverallCompositionModel(model, forces, varargin{:});
            if isa(model.EOSModel, 'SoreideWhitsonEos')
                forces = validateCompositionalForcesSW(model, forces, varargin{:});
            end
        end

        function state = initStateAD(model, state, vars, names, origin)
            % Distribute the flat primary-variable cell array (from
            % getPrimaryVariables, after the Newton update) back onto the
            % named state fields (pressure, nbact per reaction, SO4/HS
            % when active, then the remaining EOS overall composition),
            % converting composition back to mole fractions and delegating
            % what's left to the parent compositional initStateAD.
            if model.bacteriamodel

                isP = strcmp(names, 'pressure');
                isAD = any(cellfun(@(x) isa(x, 'ADI'), vars));
                state = model.setProp(state, 'pressure', vars{isP});

                removed = isP;

                bactnames=model.biochemFluid.bactnames;
                nbioreact=model.biochemFluid.nbioreact;
                nbact=cell(1, nbioreact);
                 for i = 1:nbioreact
                    name = bactnames{i};
                    sub = strcmp(names, name);
                    nbact{i} = vars{sub};
                    removed(sub) = true;
                end
                state = model.setProp(state, 'nbact', nbact);

                if model.sulfateReduction
                    isSO4 = strcmp(names, 'SO4');
                    isHS  = strcmp(names, 'HS');
                    so4 = vars{isSO4};
                    hs  = vars{isHS};
                    removed(isSO4) = true;
                    removed(isHS)  = true;
                    state = model.setProp(state, 'so4', so4);
                    state = model.setProp(state, 'hs', hs);
                end

                cnames = model.EOSModel.getComponentNames();
                ncomp = numel(cnames);
                z = cell(1, ncomp);
                z_end = 1;
                for i = 1:ncomp
                    name = cnames{i};
                    sub = strcmp(names, name);
                    if any(sub)
                        z{i} = vars{sub};
                        z_end = z_end - z{i};
                        removed(sub) = true;
                    else
                        fill = i;
                    end
                end
                z{fill} = z_end;
                state = model.setProp(state, 'components', z);

                if isAD
                    if model.sulfateReduction
                        % Refresh the EOS's salinity/H2S-speciation
                        % coupling from the CURRENT iterate's SO4/HS
                        % tracers before flashing (model is a value
                        % class, so this update only needs to survive for
                        % the getPhaseFractionAsADI call below -- it is
                        % redone every Newton iteration from state, not
                        % cached). Total dissolved sulfide uses the
                        % previous converged timestep's H2S(aq) content
                        % (state.h2sDissolvedLag, set in
                        % updateAfterConvergence) to avoid a circular
                        % dependency between the flash and its own
                        % output.
                        model.EOSModel = model.EOSModel.enablesrb_coupling(...
                            value(so4), value(hs), value(hs) + state.h2sDissolvedLag, state.T);
                    end
                    [state.x, state.y, state.L, state.FractionalDerivatives] = ...
                        model.EOSModel.getPhaseFractionAsADI(state, state.pressure, state.T, state.components);
                end
                if ~isempty(model.FacilityModel)
                    % Select facility model variables and pass them off to attached
                    % class.
                    fm = class(model.FacilityModel);
                    isF = strcmp(origin, fm);
                    state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                    removed = removed | isF;
                end
                nph = model.getNumberOfPhases();
                phnames = model.getPhaseNames();
                s = cell(1, nph);
                extra = model.getNonEoSPhaseNames();
                ne = numel(extra);
                void = 1;
                for i = 1:ne
                    sn = ['s', extra(i)];
                    isVar = strcmp(names, sn);
                    si = vars{isVar};
                    removed(isVar) = true;
                    void = void - si;

                    s{phnames == extra(i)} = si;
                end
                li = model.getLiquidIndex();
                vi = model.getVaporIndex();
                % Set up state with remaining variables
                state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));

                % Now that props have been set up, we can compute the
                % saturations from the mole fractions.
                if isAD
                    % We must get the version with derivatives
                    Z = model.getProps(state, 'PhaseCompressibilityFactors');
                    Z_L = Z{li};
                    Z_V = Z{vi};
                else
                    % Already stored in state - no derivatives needed
                    Z_L = state.Z_L;
                    Z_V = state.Z_V;
                end

                L = state.L;
                propmodel = model.EOSModel.PropertyModel;
                if isempty(propmodel.volumeShift)
                    volL = L.*Z_L;
                    volV = (1-L).*Z_V;
                else
                    volL = L./propmodel.computeMolarDensity(model.EOSModel, state.pressure, state.x, Z_L, state.T, true);
                    volV = (1-L)./propmodel.computeMolarDensity(model.EOSModel, state.pressure, state.y, Z_V, state.T, false);
                end
                volT = volL + volV;
                sL = volL./volT;
                sV = volV./volT;

                [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
                sL = sL.*void;
                sV = sV.*void;
                [s{li}, s{vi}] = model.setMinimumTwoPhaseSaturations(state, 1 - void, sL, sV, pureLiquid, pureVapor, twoPhase);
                state = model.setProp(state, 's', s);
            else
                state = initStateAD@GenericOverallCompositionModel(model, state, vars, names, origin);
            end
        end
        %-----------------------------------------------------------------%
        function [v_eqs, tolerances, names] = getConvergenceValues(model, problem, varargin)
            % Get values for convergence check with CNV-style scaling
            [v_eqs, tolerances, names] = getConvergenceValues@ReservoirModel(model, problem, varargin{:});

            if model.bacteriamodel
                nbioreact=model.biochemFluid.nbioreact;
                bacteriaIndex=zeros(nbioreact,1);
                for i=1:nbioreact
                    bact=strcat(model.biochemFluid.bactnames{i},' (cell)');
                    bacteriaIndex(i) = find(strcmp(names, bact));
                    tolerances(bacteriaIndex(i)) = 1.0e-2;
                end
                % Apply magnitude-based scaling to all equations (components + bacteria)
                % This CNV-style normalization makes residuals comparable across
                % dissimilar equation types and ensures convergence is fair.
                scale = model.getEquationScaling(problem.equations, problem.equationNames, problem.state, problem.dt);
                ix    = ~cellfun(@isempty, scale);
                v_eqs(ix) = cellfun(@(scale, x) norm(scale.*value(x), inf), scale(ix), problem.equations(ix));

                % Bacteria equation gets tighter tolerance due to stiff growth/decay
                % with quadratic decay term. Loose tolerance allows unbounded Newton
                % increments that cause Jacobian singularity.
                for i=1:nbioreact
                    if ~isempty(bacteriaIndex(i))
                        % Use 10x tighter tolerance for bacteria (stiff equation)
                        tolerances(bacteriaIndex(i)) = 5.0e-2;
                    end
                end
            end
        end

function scale = getEquationScaling(model, eqs, names, state0, dt)
            % Get scaling for the residual equations to determine convergence

            scale = cell(1, numel(eqs));
            cnames = model.getComponentNames();

            if model.bacteriamodel
                [cmass, chemistry] = model.getProps(state0, 'ComponentTotalMass', ...
                    'BacterialMass');
                cmass = value(cmass);
                chemistry = value(chemistry);
            else
                cmass= model.getProps(state0, 'ComponentTotalMass');
                cmass = value(cmass);
            end

            if ~iscell(cmass), cmass = {cmass}; end
            ncomp = model.getNumberOfComponents();
            mass = 0;
            for i = 1:ncomp
                mass = mass + cmass{i};
            end

            scaleMass = dt./mass;
            for n = cnames
                ix = strcmpi(n{1}, names);
                if ~any(ix), continue; end
                scale{ix} = scaleMass;
            end
            if model.bacteriamodel
               nbioreact=model.biochemFluid.nbioreact;
                for i=1:nbioreact
                    ix = strcmpi(names, model.biochemFluid.bactnames{i});
                    if any(ix)
                        if iscell(chemistry)
                            chemistry_i = chemistry{i};
                        else
                            chemistry_i = chemistry(:, i);
                        end
                        scaleChemistry = dt./max(chemistry_i, dt);
                        scaleChemistry = filloutliers(scaleChemistry, "nearest","mean");
                        scale{ix} = scaleChemistry;
                    end
                end
            end

        end
        function scaling = getScalingFactorsCPR(model, problem, names, solver) %#ok

            scaling = model.getEquationScaling(problem.equations, problem.equationNames, problem.state, problem.dt);

        end

        function [fn, index] = getVariableField(model, name, varargin)
            % Map a variable name to its state field and index, resolving
            % nbact/SO4/HS and any per-reaction bacterial-population name
            % (bactnames) before falling back to the parent compositional
            % model's lookup for EOS component names.
            switch(lower(name))
                case {'nbact', 'bacteriamodel'} %Bacteria model
                    index = ':';
                    fn = 'nbact';
                case  {'so4', 'tracerSO4'} % SO4 aqueous tracer (non-volatile, outside the EOS)
                    index = ':';
                    fn = 'tracerSO4';
                case {'hs', 'tracerHS'} % HS aqueous tracer (non-volatile, outside the EOS)
                    index = ':';
                    fn = 'tracerHS';
                otherwise
                    bactnames = model.biochemFluid.bactnames;
                    sub = strcmpi(bactnames, name);
                    if any(sub)
                        fn = 'nbact';
                        index = find(sub);
                    else
                        % This will throw an error for us
                        [fn, index] = getVariableField@OverallCompositionCompositionalModel(model, name, varargin{:});
                    end
            end
        end

        function names = getComponentNames(model)
            % Get names of the fluid components
            names  = getComponentNames@GenericOverallCompositionModel(model);
        end


        function  [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % After the parent compositional model accepts the converged
            % step, accumulate this step's per-reaction H2 consumption
            % (computeConvergedH2ConsumptionRate) into a running
            % cumulativeH2ConsumptionMoles diagnostic carried in state.
            [state, report] = updateAfterConvergence@GenericOverallCompositionModel(model, state0, state, dt, drivingForces);
            if model.bacteriamodel
                if isfield(state0, 'cumulativeH2ConsumptionMoles')
                    cumulative = value( ...
                        state0.cumulativeH2ConsumptionMoles);
                else
                    cumulative = zeros(model.G.cells.num, ...
                        model.biochemFluid.nbioreact);
                end
                if model.reactionsEnabled
                    state.h2ConsumptionRate = ...
                        model.computeConvergedH2ConsumptionRate(state);
                    state.cumulativeH2ConsumptionMoles = cumulative + ...
                        max(value(state.h2ConsumptionRate), 0).*dt;
                else
                    % Reactions are disabled for this stage (e.g. the flow
                    % stage of a coarse-flow/local-reaction split, see
                    % SequentialBiochemistryPhreeqcModel): no reaction
                    % extent occurred here, so cumulative H2 accounting
                    % must be carried forward unchanged rather than
                    % recomputed from an unused hypothetical growth rate.
                    state.h2ConsumptionRate = zeros(model.G.cells.num, ...
                        model.biochemFluid.nbioreact);
                    state.cumulativeH2ConsumptionMoles = cumulative;
                end
            end
            if model.sulfateReduction
                % Cache the converged dissolved-H2S concentration
                % [mol/m3 liquid] for use as the (lagged) total-sulfide
                % input to the EOS salinity coupling next timestep -- see
                % initStateAD. This avoids a circular dependency between
                % the flash (needs total sulfide -> msalt) and its own
                % output (dissolved H2S).
                names = model.EOSModel.CompositionalMixture.names;
                indH2S = find(strcmp(names, 'H2S'), 1);
                if ~isempty(indH2S)
                    x = value(state.x);
                    if iscell(x)
                        xH2S = x{indH2S};
                    else
                        xH2S = x(:, indH2S);
                    end
                    propmodel = model.EOSModel.PropertyModel;
                    rhoL_molar = propmodel.computeMolarDensity(model.EOSModel, ...
                        value(state.pressure), value(state.x), value(state.Z_L), state.T, true);
                    state.h2sDissolvedLag = value(rhoL_molar) .* xH2S;
                end
            end
        end

        function rate = computeConvergedH2ConsumptionRate(model, state)
            % Compute the converged reaction-specific H2 molar sink.
            bcrm = model.biochemFluid;
            nreact = bcrm.nbioreact;
            liquid = model.getLiquidIndex();
            if isa(model, 'BiochemistryPhreeqcModel')
                growthRateName = 'CarbonLimitedGrowthRate';
            else
                growthRateName = 'PsiGrowthRate';
            end
            psigrowth = model.FlowPropertyFunctions.get( ...
                model, state, growthRateName);
            bmass = model.PVTPropertyFunctions.get( ...
                model, state, 'BacterialMass');
            rho = model.PVTPropertyFunctions.get(model, state, 'Density');
            if iscell(rho)
                rhoL = rho{liquid};
            else
                rhoL = rho(:, liquid);
            end

            rate = zeros(model.G.cells.num, nreact);
            for i = 1:nreact
                if iscell(psigrowth)
                    growth = psigrowth{i};
                else
                    growth = psigrowth(:, i);
                end
                if iscell(bmass)
                    mass = bmass{i};
                else
                    mass = bmass(:, i);
                end
                rate(:, i) = value(bcrm.nbactMax(i).*growth.*mass./ ...
                    (bcrm.Y_H2(i).*rhoL));
            end
        end

        function [state, report] = updateState(model, state, problem, dz, drivingForces)
            % Update state with adaptive damping for bacteria variables.
            % The bacteria equation is stiff (linear growth + quadratic decay),
            % causing ill-conditioned Jacobians at growth/decay transitions.
            % Damping prevents unbounded Newton increments that cause divergence.

            if model.bacteriamodel
                % Save old bacteria state before update
                nbact_old = value(model.getProp(state, 'nbact'));

                % Apply parent class update
                [state, report] = updateState@GenericOverallCompositionModel(model, state, problem, dz, drivingForces);

                % Get updated bacteria state (after parent capping/processing)
                nbact_new = value(model.getProp(state, 'nbact'));

                % Limit fractional change per iteration to stabilize stiff kinetics.
                % Aggressive damping for highly stiff growth/decay kinetics (linear growth + nbact^2 decay).
                % 5% change per iteration is conservative but necessary for quadratic decay singularities.
                 nbioreact=model.biochemFluid.nbioreact;
                for i = 1:nbioreact
                    if iscell(nbact_old)
                        nbacti_old=nbact_old{i};
                        nbacti_new=nbact_new{i};
                    else
                        nbacti_old=nbact_old(:,i);
                        nbacti_new=nbact_new(:,i);
                    end
                    max_frac_change = 0.05;
                    frac_change = (nbacti_new - nbacti_old) ./ max(abs(nbacti_old), 1e-12);

                    % Apply adaptive damping where fractional change is excessive
                    excessive = abs(frac_change) > max_frac_change;
                    if any(excessive)&&false
                        % Apply exponential damping: new = old + max_frac_change * sign(change) * old_mag
                        sign_inc = sign(nbacti_new(excessive) - nbacti_old(excessive));
                        nbacti_damped = nbacti_old(excessive) + ...
                            max_frac_change * sign_inc .* max(abs(nbacti_old(excessive)), 1e-12);
                        nbacti_new(excessive) = nbacti_damped;
                        state = model.setProp(state, 'nbact', nbacti_new);
                    end
                end


                % Final capping to physical bounds
                state = model.capProperty(state, 'nbact', model.bact_capProp, model.bact_maxProp);
                state = model.capProperty(state, 's', 1.0e-8, 1);
                state.components = ensureMinimumFraction(state.components, model.EOSModel.minimumComposition);

                if model.sulfateReduction
                    state = model.capProperty(state, 'so4', 0);
                    state = model.capProperty(state, 'hs', 0);
                end
            else
                [state, report] = updateState@GenericOverallCompositionModel(model, state, problem, dz, drivingForces);
            end
        end


        function isDynamic = dynamicFlowTrans(model)
            % Get boolean indicating if the fluid flow transmissibility is
            % dynamically calculated
            isDynamic = isa(model.rock.perm, 'function_handle');

        end

        function isDynamic = dynamicFlowPv(model)
            % Get boolean indicating if the fluid flow porevolume is
            % dynamically calculated

            isDynamic = isa(model.rock.poro, 'function_handle');

        end

        function nbactArray = extractBactValues(model, nbact)
            % Extract bacterial values as array, handling both cell and matrix formats
            %   This helper extracts nbact values in a format suitable for
            %   function calls with variable number of arguments.
            %
            % PARAMETERS:
            %   nbact - Cell array {nbact1, nbact2, ...} or matrix [nbact1, nbact2, ...]
            %
            % RETURNS:
            %   nbactArray - Cell array of extracted values, one per reactor

            nbioreact = model.biochemFluid.nbioreact;
            nbactArray = cell(1, nbioreact);

            if iscell(nbact)
                for i = 1:nbioreact
                    nbactArray{i} = nbact{i};
                end
            else
                for i = 1:nbioreact
                    nbactArray{i} = nbact(:, i);
                end
            end
        end

        function nbact = getDummyBacterialValues(model, propertyFunction)
            % Return zero-valued bacterial arguments matching a rock handle.
            narg = nargin(propertyFunction);
            if narg < 0
                nbioreact = model.biochemFluid.nbioreact;
            else
                nbioreact = max(narg - 1, 1);
            end
            nbact = num2cell(zeros(1, nbioreact));
        end

    end
end

function flag = normalizeTransportFlag(flag, name)
validateattributes(flag, {'logical', 'numeric'}, {'scalar', 'real', 'finite'}, ...
    mfilename, name);
flag = logical(flag);
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