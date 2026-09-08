classdef BiochemistryPhreeqcModel < BiochemistryModel
    % Extends BiochemistryModel with optional, timestep-level coupling to
    % an external PHREEQC geochemistry engine.
    %
    % SYNOPSIS:
    %   model = BiochemistryPhreeqcModel(G, rock, fluid)
    %   model = BiochemistryPhreeqcModel(G, rock, fluid, compFluid)
    %   model = BiochemistryPhreeqcModel(..., 'pn1', vn1, ...)
    %
    % DESCRIPTION:
    %   BiochemistryPhreeqcModel adds an optional PHREEQC coupling on top
    %   of BiochemistryModel's compositional-flow + Monod-kinetics base.
    %   With phreeqcTimestepCoupling = false (the default) it behaves
    %   exactly like BiochemistryModel. With it enabled, one of two
    %   mutually exclusive backends (phreeqcBackend) is used after every
    %   accepted timestep:
    %
    %     - 'sequential-compositional-phreeqc': PHREEQC owns the MET/ACE/
    %       SRB kinetics and biomass entirely (via runSequentialCompositionalPhreeqcCoupling).
    %       This model's own bacterial reaction source terms
    %       (BactConvertionRate) become a no-op for this backend to avoid
    %       double counting -- see isSequentialCompositionalPhreeqcBackend.
    %
    %     - 'sequential-h2biochem-phreeqc': MRST retains ownership of the
    %       Monod kinetics and biomass transport (as in BiochemistryModel),
    %       and PHREEQC is instead used to periodically re-equilibrate
    %       aqueous/mineral chemistry -- pH, carbonate speciation, and
    %       dolomite/calcite/anhydrite buffering -- via
    %       runSequentialH2BiochemPhreeqcEquilibrium. The result is cached
    %       in sequentialH2BiochemPhreeqcChemistryFeedback and fed back
    %       into the next kinetics evaluation. This backend is driven
    %       either by the outer Picard loop in
    %       simulateSequentialH2BiochemPhreeqc, or, without any outer
    %       iteration, as the two internal building blocks of
    %       SequentialBiochemistryPhreeqcModel's coarse-flow/local-reaction
    %       split (see that class).
    %
    %   It also adds an optional fixed-pH carbonate buffer
    %   (carbonateBuffer/carbonateBufferPH) for cases that need a
    %   simplified HCO3-/CO2 chemistry without a full PHREEQC coupling,
    %   and Ca/Mg aqueous tracers used by the PHREEQC backends.
    %
    % REQUIRED PARAMETERS:
    %   G         - Simulation grid
    %   rock      - Rock properties for the model
    %   fluid     - Fluid model for the simulation
    %   compFluid - Compositional fluid mixture (optional)
    %
    % OPTIONAL PARAMETERS:
    %   'phreeqcTimestepCoupling' - Enable the post-timestep PHREEQC
    %                               coupling (default false; requires
    %                               Windows and a registered IPhreeqcCOM
    %                               server).
    %   'phreeqcBackend'          - 'sequential-compositional-phreeqc'
    %                               (default) or 'sequential-h2biochem-phreeqc'.
    %   'phreeqcDatabaseFile'     - Absolute path to the PHREEQC database
    %                               (PHREEQC_Modified.DAT).
    %   'phreeqcComProgId'        - Registered IPhreeqcCOM server ProgID
    %                               (default 'IPhreeqcCOM.Object').
    %   'phreeqcCouplingOptions'  - Struct of backend-specific coupling
    %                               parameters (kinetic rate constants,
    %                               brine composition, tolerances, ...);
    %                               see getCouplingOptions in the
    %                               corresponding utils/run*.m file for the
    %                               full field list and defaults.
    %   'carbonateBuffer'         - Enable the simplified fixed-pH HCO3-/
    %                               CO2 buffer (default false).
    %   'bacterialDecayOrder'     - 1 (first-order) or 2 (legacy
    %                               quadratic, default) biomass decay.
    %
    % RETURNS:
    %   model - BiochemistryPhreeqcModel class instance
    %
    % SEE ALSO:
    %   BiochemistryModel, SequentialBiochemistryPhreeqcModel,
    %   convertToSequentialBiochemistryPhreeqcModel,
    %   simulateSequentialH2BiochemPhreeqc, ReservoirModel,
    %   ThreePhaseCompositionalModel

    properties
        carbonateBuffer = false;          % Fixed-pH HCO3-/CO2 buffer
        carbonateBufferPH = 6.24;
        carbonateBufferPka1 = 6.35;
        phreeqcTimestepCoupling = false; % Post-timestep PHREEQC equilibration
        phreeqcBackend = 'sequential-compositional-phreeqc'; % or 'sequential-h2biochem-phreeqc'
        phreeqcDatabaseFile = '';          % Absolute path to PHREEQC_Modified.DAT
        phreeqcComProgId = 'IPhreeqcCOM.Object'; % Registered IPhreeqcCOM server ProgID
        phreeqcCouplingOptions = struct(); % Backend-specific coupling parameters (kinetics, brine, tolerances)
        % Set only by simulateSequentialH2BiochemPhreeqc. Keeping the
        % automatic post-step hook disabled prevents an accidental,
        % lagged one-pass chemistry update for this backend.
        sequentialH2BiochemPhreeqcPicardActive = false;
        % Set only internally by SequentialBiochemistryPhreeqcModel on the
        % flow-stage/reaction-stage sub-models it builds. Like
        % sequentialH2BiochemPhreeqcPicardActive, this bypasses the
        % "direct simulateScheduleAD is unsafe" guard for the
        % sequential-h2biochem-phreeqc backend, but for the general
        % coarse-flow/local-reaction split rather than the outer Picard
        % scheme.
        sequentialSplitActive = false;
        bacterialDecayOrder = 2;          % 2: legacy quadratic, 1: first-order
    end

    properties (SetAccess = private)
        % Numerical PHREEQC snapshot used only by the sequential-h2biochem-phreeqc
        % outer Picard iteration. It is deliberately model metadata, not
        % a state variable, so it has no accumulation or AD derivatives.
        sequentialH2BiochemPhreeqcChemistryFeedback = [];
    end

    methods
        %-----------------------------------------------------------------%
        function model = BiochemistryPhreeqcModel(G, rock, fluid, compFluid, biochemFluid, includeWater, backend, varargin)
            % Constructor
            baseOptions = getBiochemistryModelOptions(varargin);
            model = model@BiochemistryModel(G, rock, fluid, compFluid, ...
                biochemFluid, includeWater, backend, baseOptions{:});
            model = merge_options(model, varargin{:});
            model.molecularDiffusion = normalizeTransportFlag( ...
                model.molecularDiffusion, 'molecularDiffusion');
            model.molecularDispersion = normalizeTransportFlag( ...
                model.molecularDispersion, 'molecularDispersion');
            model.bactDiffusion = normalizeTransportFlag( ...
                model.bactDiffusion, 'bactDiffusion');
            model.chemotaxisEffect = normalizeTransportFlag( ...
                model.chemotaxisEffect, 'chemotaxisEffect');
            model.carbonateBuffer = normalizeTransportFlag( ...
                model.carbonateBuffer, 'carbonateBuffer');
            model.phreeqcTimestepCoupling = normalizeTransportFlag( ...
                model.phreeqcTimestepCoupling, 'phreeqcTimestepCoupling');
            model.phreeqcBackend = normalizePhreeqcBackend(model.phreeqcBackend);
            model.sequentialH2BiochemPhreeqcPicardActive = normalizeTransportFlag( ...
                model.sequentialH2BiochemPhreeqcPicardActive, ...
                'sequentialH2BiochemPhreeqcPicardActive');
            assert(~model.sequentialH2BiochemPhreeqcPicardActive, ...
                ['sequentialH2BiochemPhreeqcPicardActive is reserved for ', ...
                 'simulateSequentialH2BiochemPhreeqc and cannot be ', ...
                 'configured on construction.']);
            model.sequentialSplitActive = normalizeTransportFlag( ...
                model.sequentialSplitActive, 'sequentialSplitActive');
            assert(~model.sequentialSplitActive, ...
                ['sequentialSplitActive is reserved for the flow-stage ', ...
                 'and reaction-stage sub-models built internally by ', ...
                 'SequentialBiochemistryPhreeqcModel and cannot be ', ...
                 'configured on construction.']);
            assert(ischar(model.phreeqcComProgId) || ...
                (isstring(model.phreeqcComProgId) && isscalar(model.phreeqcComProgId)), ...
                'phreeqcComProgId must be a character vector or scalar string.');
            model.phreeqcComProgId = char(model.phreeqcComProgId);
            assert(~isempty(strtrim(model.phreeqcComProgId)), ...
                'phreeqcComProgId must identify a registered IPhreeqcCOM server.');
            validateattributes(model.bacterialDecayOrder, {'numeric'}, ...
                {'scalar', 'integer', '>=', 1, '<=', 2}, ...
                mfilename, 'bacterialDecayOrder');
            if model.phreeqcTimestepCoupling
                assert(ispc, ['COM PHREEQC backends require Windows and ', ...
                    'a registered IPhreeqcCOM server.']);
                validateModifiedComConfiguration(model);
            end

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
            if model.phreeqcTimestepCoupling
                assert(model.bacteriamodel && model.carbonateBuffer && model.sulfateReduction, ...
                    ['phreeqcTimestepCoupling requires bacterial HCO3 and ', ...
                     'SO4/HS tracer transport.']);
            end

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
            % As BiochemistryModel.validateModel, plus a guard specific to
            % the 'sequential-compositional-phreeqc' backend: spatial
            % bacterial diffusion/chemotaxis is rejected (PHREEQC's
            % kinetics are purely local, with no notion of biomass
            % transport), and a one-time notice is printed explaining that
            % MRST's own nbact equation is a no-op under this backend.
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
            if model.bacteriamodel && model.isSequentialCompositionalPhreeqcBackend()
                assert(~model.bactDiffusion && ~model.chemotaxisEffect, ...
                    ['phreeqcBackend=''sequential-compositional-phreeqc'' cannot be combined with ', ...
                     'bactDiffusion or chemotaxisEffect: PHREEQC integrates only ', ...
                     'local per-cell kinetics and has no notion of spatial ', ...
                     'bacterial transport. Diffusing/chemotaxis-moving nbact would ', ...
                     'move a quantity that is disconnected from the biomass PHREEQC ', ...
                     'is actually growing in sequentialCompositionalPhreeqcBiomassMET/ACE/SRB.']);
                fprintf(['BiochemistryPhreeqcModel: phreeqcBackend=''sequential-compositional-phreeqc'' is active -- ', ...
                    'MET/ACE/SRB kinetics (growth/decay and component sources) are ', ...
                    'carried out by PHREEQC. MRST''s bacterial mass-balance equation ', ...
                    '(nbact) is not assembled: with zero MRST reaction source and no ', ...
                    'diffusion/chemotaxis transport it would be a no-op every step. ', ...
                    'PsiGrowthRate/CarbonLimitedGrowthRate/BacterialMass remain ', ...
                    'available as diagnostic outputs only.\n']);
            end
            model = validateModel@GenericOverallCompositionModel(model, varargin{:});
        end

        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@GenericOverallCompositionModel(model, varargin{:});

            fluxprops = model.FlowDiscretization;
            pvtprops  = model.PVTPropertyFunctions;
            flowprops = model.FlowPropertyFunctions;

            if model.bacteriamodel
                flowprops = flowprops.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                flowprops = flowprops.setStateFunction('CarbonLimitedGrowthRate', ...
                    CarbonLimitedGrowthRate(model));
                flowprops = flowprops.setStateFunction('PsiDecayRate',  DecayBactRateSRC(model));
                flowprops = flowprops.setStateFunction('BactConvRate',  BactConvertionRate(model));

                % Register bacterial mass as cell property (not a source term)
                pvtprops = pvtprops.setStateFunction('BacterialMass', BacterialMass(model));

                if model.hasMobileAqueousTracers()
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
            state = validateState@ThreePhaseCompositionalModel(model, state);
            if model.bacteriamodel && ~isfield(state, 'nbact')
                nbact0 = 1e6;
                state.nbact = repmat(nbact0, model.G.cells.num, 1);
            end
            if model.carbonateBuffer && ~isfield(state, 'tracerHCO3')
                state.tracerHCO3 = zeros(model.G.cells.num, 1);
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
            if model.phreeqcTimestepCoupling
                % These are state quantities, rather than mutable model
                % properties, because updateAfterConvergence returns only
                % state. They are used by the next timestep's kinetics.
                if ~isfield(state, 'phreeqcPH')
                    state.phreeqcPH = repmat(model.carbonateBufferPH, model.G.cells.num, 1);
                end
                if ~isfield(state, 'phreeqcCarbonatePka1')
                    state.phreeqcCarbonatePka1 = repmat( ...
                        model.carbonateBufferPka1, model.G.cells.num, 1);
                end
                if ~isfield(state, 'tracerCa')
                    state.tracerCa = zeros(model.G.cells.num, 1);
                end
                if ~isfield(state, 'tracerMg')
                    state.tracerMg = zeros(model.G.cells.num, 1);
                end
            end
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            % As BiochemistryModel.getPrimaryVariables. The PHREEQC-owned
            % biomass under the 'sequential-compositional-phreeqc' backend
            % (sequentialCompositionalPhreeqcBiomass*) is plain state
            % data updated after the coupling call, not a primary
            % variable/AD unknown, so it does not appear here.
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
                if model.isSequentialCompositionalPhreeqcBackend()
                    % nbact has no equation to solve for this backend (see
                    % validateModel/isSequentialCompositionalPhreeqcBackend): it stays a
                    % plain state field, not a Newton primary variable.
                    names = [{'pressure'}, cnames(2:end), enames];
                    vars  = [p, z(2:end), evars];
                else
                    nbact = model.getProp(state, 'nbact');
                    nbact = expandMatrixToCell(nbact);
                    bactnames = model.biochemFluid.bactnames;
                    names = [{'pressure'}, cnames(2:end), bactnames, enames];
                    vars  = [p, z(2:end), nbact, evars];
                end
                aqueousNames = model.getAqueousTracerNames();
                for i = 1:numel(aqueousNames)
                    names = [names, aqueousNames(i)]; %#ok<AGROW>
                    vars  = [vars, {model.getProp(state, lower(aqueousNames{i}))}]; %#ok<AGROW>
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
            % As BiochemistryModel.getModelEquations, with two additions:
            % when carbonateBuffer is active, the pre-step aqueous carbon
            % inventory is cached (carbonSubstrateMoles/Dt) so the
            % methanogenic/acetogenic growth rate can be capped by the
            % carbon actually available (see CarbonLimitedGrowthRate); and
            % an hco3Sink accumulator tracks the carbonate drawn down by
            % the buffer, used later in this method to close its balance.
            % None of this touches PHREEQC directly -- PHREEQC coupling
            % happens after the timestep converges (see the run*Coupling
            % utilities and BiochemistryPhreeqcModel's class-level docs).
            %
            % Discretize
            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            if model.bacteriamodel && model.carbonateBuffer && model.reactionsEnabled
                % Use the pre-step aqueous inventory to bound this step's
                % growth. Using the evolving Newton state here would make
                % the rate cap consume its own final-state availability.
                state.carbonSubstrateDt = dt;
                state.carbonSubstrateMoles = getAqueousCarbonMoles(model, state0);
            end
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

            hco3Sink = 0;
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

                if model.carbonateBuffer
                    idxCO2 = find(strcmpi(cnames, 'CO2'), 1);
                    assert(~isempty(idxCO2), ...
                        'carbonateBuffer requires an EOS CO2 component.');

                    hco3 = model.getProp(state, 'hco3');
                    pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
                    s = model.getProp(state, 's');
                    if iscell(s)
                        sL = max(s{L_ix}, 0);
                    else
                        sL = max(s(:, L_ix), 0);
                    end

                    molarMassCO2 = model.EOSModel.CompositionalMixture.molarMass(idxCO2);
                    co2Demand = max(-src_rate{idxCO2}./rhoL./molarMassCO2, 0);
                    hco3Available = pv.*sL.*hco3./dt;
                    hco3Sink = min(co2Demand, hco3Available);

                    % Bicarbonate supplies the CO2 consumed by MET/ACE.
                    % src_rate is later divided by rhoL in the component
                    % equation, so apply the inverse scaling here.
                    src_rate{idxCO2} = src_rate{idxCO2} + ...
                        hco3Sink.*molarMassCO2.*rhoL;
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
                fd = model.FlowDiscretization;
                if model.isSequentialCompositionalPhreeqcBackend()
                    % nbact's mass-balance equation is a provable no-op
                    % under this backend: MET/ACE/SRB kinetics (the
                    % reaction source) are carried out entirely by
                    % PHREEQC, and validateModel forbids combining this
                    % backend with bactDiffusion/chemotaxisEffect (the
                    % only transport terms nbact has). With zero source
                    % and zero flux the equation only re-derives
                    % nbact == nbact0 every step, so skip assembling it
                    % rather than spending a Newton unknown on it.
                    beqs = {}; bnames = {}; btypes = {};
                else
                    % Bacterial mass balance: d(M)/dt + div(flux) = source
                    % where M = pv * S_l * nbact [kg]
                    [beqs, bflux, bnames, btypes] = model.FlowDiscretization.bacteriaConservationEquation(model, state, state0, dt);
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
                             %beqs{1} = model.operators.AccDiv(beqs{1},0);
                        end

                        if model.reactionsEnabled
                            beqs{i} = beqs{i} - src_growthdecay{i};
                        end
                    end
                end

                if model.hasMobileAqueousTracers()
                    % Aqueous tracers are advected with the liquid phase.
                    % SO4/HS retain their SRB sources; HCO3 carries the
                    % existing biological carbon sink, while Ca/Mg have no
                    % MRST reaction source and are updated by PHREEQC.
                    [tacc, tflux, tnames, ttypes] = model.FlowDiscretization.tracerConservationEquation(model, state, state0, dt);
                    if model.reactionsEnabled
                        src_tracer = model.FacilityModel.getAqueousTracerSources(fd, state, state0, dt);
                        hco3Index = find(strcmp(tnames, 'HCO3'), 1);
                        if ~isempty(hco3Index)
                            src_tracer{hco3Index} = src_tracer{hco3Index} - hco3Sink;
                        end
                    end
                    for i = 1:numel(tacc)
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
            forces = validateDrivingForces@GenericOverallCompositionModel(model, forces, varargin{:});
            if isa(model.EOSModel, 'SoreideWhitsonEos')
                forces = validateCompositionalForcesSW(model, forces, varargin{:});
            end
        end

        function state = initStateAD(model, state, vars, names, origin)
            % As BiochemistryModel.initStateAD, except nbact is only
            % pulled from the primary-variable vector when it actually is
            % one; under 'sequential-compositional-phreeqc' it is carried
            % through unchanged from state (see getPrimaryVariables).
            if model.bacteriamodel

                isP = strcmp(names, 'pressure');
                isAD = any(cellfun(@(x) isa(x, 'ADI'), vars));
                state = model.setProp(state, 'pressure', vars{isP});

                removed = isP;

                if model.isSequentialCompositionalPhreeqcBackend()
                    % nbact is not a primary variable for this backend
                    % (see getPrimaryVariables); carry its existing value
                    % through unchanged.
                    nbact = expandMatrixToCell(model.getProp(state, 'nbact'));
                else
                    bactnames=model.biochemFluid.bactnames;
                    nbioreact=model.biochemFluid.nbioreact;
                    nbact=cell(1, nbioreact);
                    for i = 1:nbioreact
                        name = bactnames{i};
                        sub = strcmp(names, name);
                        nbact{i} = vars{sub};
                        removed(sub) = true;
                    end
                end
                state = model.setProp(state, 'nbact', nbact);

                aqueousNames = model.getAqueousTracerNames();
                for i = 1:numel(aqueousNames)
                    isTracer = strcmp(names, aqueousNames{i});
                    if any(isTracer)
                        state = model.setProp(state, lower(aqueousNames{i}), vars{isTracer});
                        removed(isTracer) = true;
                    end
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
                        so4 = model.getProp(state, 'so4');
                        hs = model.getProp(state, 'hs');
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

            if model.bacteriamodel && ~model.isSequentialCompositionalPhreeqcBackend()
                % No nbact equation is assembled for this backend (see
                % getModelEquations/getPrimaryVariables), so there is no
                % '<name> (cell)' residual to look up here.
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
            switch(lower(name))
                case {'nbact', 'bacteriamodel'} %Bacteria model
                    index = ':';
                    fn = 'nbact';
                case  {'so4', 'tracerso4'} % SO4 aqueous tracer (non-volatile, outside the EOS)
                    index = ':';
                    fn = 'tracerSO4';
                case {'hs', 'tracerhs'} % HS aqueous tracer (non-volatile, outside the EOS)
                    index = ':';
                    fn = 'tracerHS';
                case {'hco3', 'tracerhco3'}
                    index = ':';
                    fn = 'tracerHCO3';
                case {'ca', 'tracerca'}
                    index = ':';
                    fn = 'tracerCa';
                case {'mg', 'tracermg'}
                    index = ':';
                    fn = 'tracerMg';
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
            % Like BiochemistryModel.updateAfterConvergence, but skips the
            % H2-consumption diagnostic entirely when PHREEQC owns the
            % kinetics ('sequential-compositional-phreeqc'), and, for the
            % 'sequential-h2biochem-phreeqc' backend, additionally
            % accumulates this step's consumption for the outer Picard
            % loop / coarse-flow-local-reaction substep bookkeeping.
            [state, report] = updateAfterConvergence@GenericOverallCompositionModel(model, state0, state, dt, drivingForces);
            if model.bacteriamodel && ~model.isSequentialCompositionalPhreeqcBackend()
                if model.reactionsEnabled
                    % Preserve the converged reaction source before an optional
                    % PHREEQC split step changes the chemical state.
                    state.h2ConsumptionRate = ...
                        model.computeConvergedH2ConsumptionRate(state);
                    if model.isSequentialH2BiochemPhreeqcBackend()
                        state = accumulateSequentialH2BiochemPhreeqcH2Consumption( ...
                            state0, state, dt, model.G.cells.num, ...
                            model.biochemFluid.nbioreact);
                    end
                else
                    % Reactions are disabled for this stage (e.g. the flow
                    % stage of a coarse-flow/local-reaction split, see
                    % SequentialBiochemistryPhreeqcModel): no reaction
                    % extent occurred here, so the H2 accounting is
                    % carried forward unchanged rather than recomputed or
                    % accumulated from an unused hypothetical rate.
                    state.h2ConsumptionRate = zeros(model.G.cells.num, ...
                        model.biochemFluid.nbioreact);
                    if isfield(state0, 'sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles')
                        state.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles = ...
                            value(state0.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles);
                    end
                end
            end
            if model.phreeqcTimestepCoupling
                % This deliberately runs only after the nonlinear timestep
                % has converged. The helper updates aqueous tracer and
                % mineral state for the following timestep; it never enters
                % a Newton iteration.
                switch model.phreeqcBackend
                    case 'sequential-compositional-phreeqc'
                        state = runSequentialCompositionalPhreeqcCoupling( ...
                            model, state, dt);
                    case 'sequential-h2biochem-phreeqc'
                        if ~model.sequentialH2BiochemPhreeqcPicardActive && ...
                                ~model.sequentialSplitActive
                            error('BiochemistryPhreeqcModel:SequentialH2BiochemPhreeqcRequiresPicardDriver', ...
                                ['phreeqcBackend=''sequential-h2biochem-phreeqc'' must be run with ', ...
                                 'simulateSequentialH2BiochemPhreeqc, or with a ', ...
                                 'SequentialBiochemistryPhreeqcModel-driven local ', ...
                                 'reaction stage. Direct simulateScheduleAD would ', ...
                                 'apply a lagged chemistry split and is intentionally ', ...
                                 'rejected.']);
                        end
                    otherwise
                        error('BiochemistryPhreeqcModel:InvalidPhreeqcBackend', ...
                            'Unsupported PHREEQC backend: %s.', model.phreeqcBackend);
                end
            end
            if model.sulfateReduction
                % Cache the converged dissolved-H2S concentration
                % [mol/m3 liquid] for use as the (lagged) total-sulfide
                % input to the EOS salinity coupling next timestep -- see
                % initStateAD. This avoids a circular dependency between
                % the flash (needs total sulfide -> msalt) and its own
                % output (dissolved H2S).
                % It deliberately follows the PHREEQC update and reflash.
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

                % Backtrack only updates that would cross the artificial
                % lower bound. A generic fractional cap would make the
                % converged biomass depend on the Newton iteration count.
                nbioreact = model.biochemFluid.nbioreact;
                nbact_safe = nbact_new;
                didBacktrack = false;
                lower = model.bact_capProp;
                for i = 1:nbioreact
                    if iscell(nbact_old)
                        nbacti_old = nbact_old{i};
                        nbacti_new = nbact_new{i};
                    else
                        nbacti_old = nbact_old(:, i);
                        nbacti_new = nbact_new(:, i);
                    end

                    crossesLower = nbacti_old > lower & nbacti_new <= lower;
                    if any(crossesLower)
                        alpha = 0.5*(nbacti_old(crossesLower) - lower)./ ...
                            (nbacti_old(crossesLower) - nbacti_new(crossesLower));
                        nbacti_safe_values = nbacti_old(crossesLower) + alpha.* ...
                            (nbacti_new(crossesLower) - nbacti_old(crossesLower));
                        if iscell(nbact_safe)
                            nbact_safe{i}(crossesLower) = nbacti_safe_values;
                        else
                            nbact_safe(crossesLower, i) = nbacti_safe_values;
                        end
                        didBacktrack = true;
                    end
                end
                if didBacktrack
                    state = model.setProp(state, 'nbact', nbact_safe);
                end


                % Final capping to physical bounds
                state = model.capProperty(state, 'nbact', model.bact_capProp, model.bact_maxProp);
                if model.carbonateBuffer
                    state = model.capProperty(state, 'hco3', 0, inf);
                end
                state = model.capProperty(state, 's', 1.0e-8, 1);
                state.components = ensureMinimumFraction(state.components, model.EOSModel.minimumComposition);

                if model.sulfateReduction
                    state = model.capProperty(state, 'so4', 0);
                    state = model.capProperty(state, 'hs', 0);
                end
                if model.phreeqcTimestepCoupling
                    state = model.capProperty(state, 'ca', 0);
                    state = model.capProperty(state, 'mg', 0);
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

        function names = getAqueousTracerNames(model)
            % Mobile non-EOS tracers, all stored as mol/m^3 liquid.
            names = {};
            if model.sulfateReduction
                names = {'SO4', 'HS'};
            end
            if model.carbonateBuffer
                names = [names, {'HCO3'}];
            end
            if model.phreeqcTimestepCoupling
                names = [names, {'Ca', 'Mg'}];
            end
        end

        function active = hasMobileAqueousTracers(model)
            active = model.bacteriamodel && ~isempty(model.getAqueousTracerNames());
        end

        function active = isSequentialCompositionalPhreeqcBackend(model)
            % True only while the optional split COM kinetics is active.
            %
            % The implicit h2-biochem reaction sources are disabled in
            % this mode. The nbact mass-balance equation is not assembled
            % either (validateModel forbids bactDiffusion/chemotaxisEffect
            % here, so it would otherwise be solved as a pure no-op every
            % step); aqueous tracer transport equations remain active.
            % PHREEQC owns separate biomass state in the post-convergence
            % split.
            active = model.phreeqcTimestepCoupling && ...
                strcmp(model.phreeqcBackend, 'sequential-compositional-phreeqc');
        end

        function active = isSequentialH2BiochemPhreeqcBackend(model)
            % True when equilibrium COM chemistry is coupled to MRST's
            % existing Monod reaction model through the Picard driver.
            active = model.phreeqcTimestepCoupling && ...
                strcmp(model.phreeqcBackend, 'sequential-h2biochem-phreeqc');
        end

        function model = setSequentialH2BiochemPhreeqcChemistryFeedback(model, feedback)
            % Store one immutable numerical chemistry snapshot for the
            % next Picard candidate or split reaction stage.
            assert(model.isSequentialH2BiochemPhreeqcBackend() && ...
                (model.sequentialH2BiochemPhreeqcPicardActive || ...
                 model.sequentialSplitActive), ...
                ['h2-biochem chemistry feedback is reserved for the active ', ...
                 'sequential-h2biochem-phreeqc coupling drivers.']);
            model.sequentialH2BiochemPhreeqcChemistryFeedback = ...
                normalizeSequentialH2BiochemPhreeqcChemistryFeedback(feedback, model.G.cells.num);
        end

        function model = clearSequentialH2BiochemPhreeqcChemistryFeedback(model)
            model.sequentialH2BiochemPhreeqcChemistryFeedback = [];
        end

        function feedback = getSequentialH2BiochemPhreeqcChemistryFeedback(model)
            % Return the active coupling stage's fixed numerical chemistry.
            feedback = [];
            if model.isSequentialH2BiochemPhreeqcBackend() && ...
                    (model.sequentialH2BiochemPhreeqcPicardActive || ...
                     model.sequentialSplitActive)
                feedback = model.sequentialH2BiochemPhreeqcChemistryFeedback;
                assert(~model.reactionsEnabled || ~isempty(feedback), ...
                    ['The active sequential-h2biochem-phreeqc solve is missing its ', ...
                     'immutable PHREEQC chemistry feedback snapshot.']);
            end
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

function baseOptions = getBiochemistryModelOptions(options)
% Keep subclass-only options away from the original model constructor.
subclassOptions = { ...
    'carbonateBuffer', 'carbonateBufferPH', 'carbonateBufferPka1', ...
    'phreeqcTimestepCoupling', 'phreeqcBackend', ...
    'phreeqcDatabaseFile', ...
    'phreeqcComProgId', 'phreeqcCouplingOptions', ...
    'sequentialH2BiochemPhreeqcPicardActive', 'sequentialSplitActive', ...
    'bacterialDecayOrder'};
assert(mod(numel(options), 2) == 0, ...
    'BiochemistryPhreeqcModel options must be property/value pairs.');
baseOptions = {};
for i = 1:2:numel(options)
    if ~any(strcmpi(options{i}, subclassOptions))
        baseOptions(end + 1:end + 2) = options(i:i + 1); %#ok<AGROW>
    end
end
end

function backend = normalizePhreeqcBackend(backend)
assert(ischar(backend) || (isstring(backend) && isscalar(backend)), ...
    'phreeqcBackend must be a character vector or scalar string.');
backend = lower(strtrim(char(backend)));
assert(ismember(backend, {'sequential-compositional-phreeqc', ...
    'sequential-h2biochem-phreeqc'}), ...
    ['phreeqcBackend must be ''sequential-compositional-phreeqc'' or ', ...
     '''sequential-h2biochem-phreeqc''; ', ...
     'got ''%s''.'], backend);
end

function feedback = normalizeSequentialH2BiochemPhreeqcChemistryFeedback(feedback, nc)
required = {'pH', 'carbonatePka1', 'hco3Molality', 'co2Molality', ...
    'sulfateMolality'};
assert(isstruct(feedback) && isscalar(feedback) && ...
    all(isfield(feedback, required)), ...
    ['sequential-h2biochem-phreeqc chemistry feedback must be a scalar struct with ', ...
     'pH, carbonatePka1, hco3Molality, co2Molality, and sulfateMolality.']);
for i = 1:numel(required)
    name = required{i};
    values = feedback.(name);
    assert(isnumeric(values) && isreal(values) && ~isa(values, 'ADI'), ...
        'sequential-h2biochem-phreeqc feedback %s must be a real numeric (non-AD) vector.', name);
    values = values(:);
    assert(numel(values) == nc && all(isfinite(values)), ...
        'sequential-h2biochem-phreeqc feedback %s must contain one finite value per cell.', name);
    if ~ismember(name, {'pH', 'carbonatePka1'})
        assert(all(values >= 0), ...
            'sequential-h2biochem-phreeqc feedback %s must be non-negative.', name);
    end
    feedback.(name) = double(values);
end
end

function validateModifiedComConfiguration(model)
databaseFile = model.phreeqcDatabaseFile;
if isempty(databaseFile) && isfield(model.phreeqcCouplingOptions, 'databaseFile')
    databaseFile = model.phreeqcCouplingOptions.databaseFile;
end
assert(ischar(databaseFile) || (isstring(databaseFile) && isscalar(databaseFile)), ...
    'COM PHREEQC backends require an explicitly configured PHREEQC_Modified.DAT path.');
databaseFile = char(databaseFile);
assert(~isempty(strtrim(databaseFile)) && isAbsolutePhreeqcPath(databaseFile), ...
    ['COM PHREEQC backends require an explicit absolute databaseFile path to ', ...
     'PHREEQC_Modified.DAT.']);
assert(isfile(databaseFile), ...
    'COM PHREEQC database not found: %s', databaseFile);
assert(contains(lower(databaseFile), 'phreeqc_modified.dat'), ...
    ['COM PHREEQC backends require PHREEQC_Modified.DAT, not a standard ', ...
     'PHREEQC database: %s'], databaseFile);

comProgId = model.phreeqcComProgId;
if isfield(model.phreeqcCouplingOptions, 'comProgId')
    comProgId = model.phreeqcCouplingOptions.comProgId;
end
assert(ischar(comProgId) || (isstring(comProgId) && isscalar(comProgId)), ...
    'COM PHREEQC backends require an explicitly configured IPhreeqcCOM ProgID.');
assert(~isempty(strtrim(char(comProgId))), ...
    'COM PHREEQC backends require an explicitly configured IPhreeqcCOM ProgID.');
end

function isAbsolute = isAbsolutePhreeqcPath(path)
path = char(path);
isAbsolute = ~isempty(regexp(path, '^[A-Za-z]:[\\/]|^\\\\', 'once')) || ...
    startsWith(path, filesep);
end

function carbonMoles = getAqueousCarbonMoles(model, state)
% Return the finite aqueous CO2 plus HCO3 inventory in each cell.
cnames = model.EOSModel.getComponentNames();
idxCO2 = find(strcmpi(cnames, 'CO2'), 1);
assert(~isempty(idxCO2), 'carbonateBuffer requires an EOS CO2 component.');
assert(isfield(state, 'Z_L'), ...
    'carbonateBuffer requires the liquid EOS Z-factor.');

x = model.getProp(state, 'x');
if iscell(x)
    xCO2 = x{idxCO2};
else
    xCO2 = x(:, idxCO2);
end

liquid = model.getLiquidIndex();
s = model.getProp(state, 's');
if iscell(s)
    sL = max(s{liquid}, 0);
else
    sL = max(s(:, liquid), 0);
end

rhoMolarL = model.EOSModel.PropertyModel.computeMolarDensity( ...
    model.EOSModel, state.pressure, x, state.Z_L, state.T, true);
poreVolume = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
hco3 = max(model.getProp(state, 'hco3'), 0);
carbonMoles = poreVolume.*sL.*(rhoMolarL.*xCO2 + hco3);
end

function state = accumulateSequentialH2BiochemPhreeqcH2Consumption(state0, state, dt, nc, nreact)
% Keep an exact candidate-step reaction extent across adaptive substeps.
%
% NonLinearSolver calls updateAfterConvergence once for each accepted
% ministep. The instantaneous final rate alone cannot represent the
% reaction extent of the enclosing nominal Picard timestep.
if isfield(state0, 'sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles')
    cumulative = value(state0.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles);
else
    cumulative = zeros(nc, nreact);
end
assert(isequal(size(cumulative), [nc, nreact]) && ...
    all(isfinite(cumulative(:)) & cumulative(:) >= 0), ...
    'sequential-h2biochem-phreeqc cumulative H2 consumption has invalid dimensions or values.');
state.sequentialH2BiochemPhreeqcCumulativeH2ConsumptionMoles = ...
    cumulative + max(value(state.h2ConsumptionRate), 0).*dt;
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