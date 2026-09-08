classdef GrowthBactRateSRC < StateFunction
    % Bacterial growth rate computation for compositional simulations
    %
    % SYNOPSIS:
    %   gr = GrowthBactRateSRC(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the specific growth rate coefficient (kinetic term only)
    %   for each microbial population using Monod kinetics.
    %
    %   For methanogens and acetogens, H2 and CO2 are taken from the
    %   liquid-phase EOS composition (x). For sulfate reducers (SRB),
    %   the substrate SO4 is taken from the aqueous tracer `tracerSO4`
    %   and converted to mole fraction using the liquid molar density.
    %   During an active sequential-h2biochem-phreeqc Picard solve, the PHREEQC pH,
    %   DIC/CO2, and sulfate snapshot replaces only those kinetic
    %   substrate terms; all conserved-state balances remain unchanged.

    properties
        % No additional properties
    end

    methods
        function gp = GrowthBactRateSRC(model, varargin)
            % Constructor
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('x', 'state');
            % If SRB is active, we need the tracer and liquid Z-factor
            if isprop(model, 'sulfateReduction') && model.sulfateReduction
                gp = gp.dependsOn('tracerSO4', 'state');
                gp = gp.dependsOn('Z_L', 'state');
            end
            if isprop(model, 'carbonateBuffer') && model.carbonateBuffer
                gp = gp.dependsOn('tracerHCO3', 'state');
                gp = gp.dependsOn('Z_L', 'state');
            end
            rm = model;
            if isprop(model, 'ReservoirModel') && ~isempty(model.ReservoirModel)
                rm = model.ReservoirModel;
            end
            if isprop(rm, 'phreeqcTimestepCoupling') && rm.phreeqcTimestepCoupling
                gp = gp.dependsOn('phreeqcPH', 'state');
                gp = gp.dependsOn('phreeqcCarbonatePka1', 'state');
                gp = gp.dependsOn('phreeqcHCO3Molality', 'state');
            end
            gp.label = '\Psi_{growth}';
        end

        function Psigrowth = evaluateOnDomain(prop, model, state)
            % Compute specific growth rate coefficient [1/s]
            %
            % Returns: Psigrowthmax * axH2 * axsub [1/s]
            % Used with BacterialMass: source = Psigrowth * BacterialMass

            if isprop(model, 'ReservoirModel') && ~isempty(model.ReservoirModel)
                rm = model.ReservoirModel;
            else
                rm = model;
            end
            bcrm = rm.biochemFluid;
            namecp = rm.getComponentNames();
%            namecp = model.EOSModel.getComponentNames();
            nbioreact = bcrm.nbioreact;
            feedback = [];
            if ismethod(rm, 'getSequentialH2BiochemPhreeqcChemistryFeedback')
                feedback = rm.getSequentialH2BiochemPhreeqcChemistryFeedback();
            end

            % Initialize output
            Psigrowth = cell(1, nbioreact);
            [Psigrowth{:}] = deal(0);

            % Get liquid mole fractions (EOS components)
            x = rm.getProp(state, 'x');
            h2Active = [];
            if isa(rm, 'BiochemistryPhreeqcModel') && ...
                    rm.isSequentialH2BiochemPhreeqcBackend()
                threshold = getH2ActivationThreshold(rm);
                idxH2Overall = find(strcmpi(namecp, 'H2') | ...
                    strcmpi(namecp, 'Hydrogen'), 1);
                assert(~isempty(idxH2Overall), ...
                    'Hybrid PHREEQC kinetics require an H2 EOS component.');
                overallH2 = state.components(:, idxH2Overall);
                h2Active = value(overallH2) > threshold;
            end

            % Loop over reactions
            for i = 1:nbioreact
                % Find H2 index (required for all reactions)
                idx_H2 = find(strcmpi(namecp, bcrm.rH2(i)), 1);
                if isempty(idx_H2)
                    continue;   % H2 not found – skip this reaction
                end
                % H2 mole fraction
                if iscell(x)
                    xH2 = x{idx_H2};
                else
                    xH2 = x(:, idx_H2);
                end

                % Handle substrate depending on reaction type
                if strcmp(bcrm.metabolicReaction(i), 'SulfateReducingBacteria')
                    % SRB: substrate is sulfate tracer (not in EOS)
                    if ~isempty(feedback)
                        so4_conc = feedback.sulfateMolality .* rm.EOSModel.rho_water;
                    elseif isfield(state, 'tracerSO4')
                        so4_conc = state.tracerSO4;   % mol/m3 liquid
                    else
                        so4_conc = zeros(size(xH2));
                    end
                    % Compute liquid molar density from Z_L
                    if isfield(state, 'Z_L')
                        Z_L = state.Z_L;
                    else
                        Z_L = ones(size(xH2));
                    end
                    R = 8.314;   % Pa·m3/(mol·K)
                    % pressure, T may be ADI objects; use value() if needed
                    P = state.pressure;
                    T = state.T;
                    rho_molar = P ./ (Z_L .* R .* T);   % mol/m3
                    % Convert tracer to mole fraction
                    xsub = so4_conc ./ rho_molar;
                    % Use the half‑saturation constant for sulfate (already in mol/mol)
                    alphasub = bcrm.alphasub(i);
                else
                    % Methanogens and acetogens: substrate is an EOS component
                    idx_sub = find(strcmpi(namecp, bcrm.rsub(i)), 1);
                    if isempty(idx_sub)
                        continue;
                    end
                    if iscell(x)
                        xsub = x{idx_sub};
                    else
                        xsub = x(:, idx_sub);
                    end
                    alphasub = bcrm.alphasub(i);

                    if isprop(rm, 'carbonateBuffer') && rm.carbonateBuffer && ...
                            strcmpi(bcrm.rsub(i), 'CO2') && ...
                            (~isempty(feedback) || isfield(state, 'tracerHCO3'))
                        if isfield(state, 'Z_L')
                            Z_L = state.Z_L;
                        else
                            Z_L = ones(size(xsub));
                        end
                        R = 8.314;
                        rho_molar = state.pressure ./ (Z_L .* R .* state.T);
                        if ~isempty(feedback)
                            hco3 = feedback.hco3Molality .* rm.EOSModel.rho_water;
                            pH = feedback.pH;
                            pKa1 = feedback.carbonatePka1;
                            if isfield(feedback, 'totalCarbonMolality')
                                totalCarbon = feedback.totalCarbonMolality;
                            else
                                totalCarbon = feedback.co2Molality + ...
                                    feedback.hco3Molality;
                            end
                            % UGFACT/PHREEQC MET and ACE kinetics use
                            % tot("Carbonate(4)"), not free aqueous CO2.
                            xsub = totalCarbon .* ...
                                rm.EOSModel.rho_water ./ rho_molar;
                        else
                            hco3 = state.tracerHCO3;
                            % tracerHCO3 transports non-CO2 DIC so that
                            % carbonate species omitted from the EOS are not
                            % discarded. For pH-dependent kinetics use the
                            % actual bicarbonate species returned by PHREEQC.
                            if isfield(state, 'phreeqcHCO3Molality')
                                hco3 = state.phreeqcHCO3Molality .* rm.EOSModel.rho_water;
                            end
                            if isfield(state, 'phreeqcPH')
                                pH = state.phreeqcPH;
                            else
                                pH = rm.carbonateBufferPH;
                            end
                            if isfield(state, 'phreeqcCarbonatePka1')
                                pKa1 = state.phreeqcCarbonatePka1;
                            else
                                pKa1 = rm.carbonateBufferPka1;
                            end
                        end
                        if isempty(feedback)
                            hydrogenIon = 10.^(-pH);
                            Ka1 = 10.^(-pKa1);
                            co2FromBuffer = hco3 ./ rho_molar .* ...
                                hydrogenIon ./ Ka1;
                            xsub = max(xsub, co2FromBuffer);
                        end
                    end
                end

                % Now compute Monod terms
                alphaH2 = bcrm.alphaH2(i);
                Psigrowthmax = bcrm.Psigrowthmax(i);

                axH2 = xH2 ./ (alphaH2 + xH2);
                axsub = xsub ./ (alphasub + xsub);

                environmentalResponse = 1;
                if ~isempty(feedback)
                    environmentalResponse = microbialEnvironmentalResponse( ...
                        bcrm.metabolicReaction{i}, value(state.T) - 273.15, ...
                        feedback.pH, feedback.tds);
                end
                Psigrowth{i} = Psigrowthmax .* axH2 .* axsub .* ...
                    environmentalResponse;
                if ~isempty(h2Active)
                    Psigrowth{i} = Psigrowth{i} .* h2Active;
                end
            end
        end
    end
end

function threshold = getH2ActivationThreshold(model)
threshold = 1e-3;
options = model.phreeqcCouplingOptions;
name = 'sequentialH2BiochemPhreeqcH2ActivationThreshold';
if isfield(options, name)
    threshold = options.(name);
end
validateattributes(threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative', '<=', 1}, ...
    mfilename, name);
end

function response = microbialEnvironmentalResponse(reaction, temperatureC, pH, tds)
switch reaction
    case 'MethanogenicArchae'
        temperatureLimits = [10, 45, 122];
        pHLimits = [4.1, 7.7, 10.2];
    case 'SulfateReducingBacteria'
        temperatureLimits = [10, 48, 113];
        pHLimits = [1, 7.01, 11.5];
    case 'AcetogenicBacteria'
        temperatureLimits = [25, 38, 72];
        pHLimits = [3.6, 7.04, 9.5];
    otherwise
        error('GrowthBactRateSRC:UnsupportedEnvironmentalResponse', ...
            'No PHREEQC environmental response is defined for reaction "%s".', reaction);
end

temperatureResponse = asymmetricParabolicResponse( ...
    temperatureC, temperatureLimits);
pHResponse = asymmetricParabolicResponse(pH, pHLimits);

tdsResponse = ones(size(tds));
salinityLimited = tds >= 50 & tds <= 300;
tdsResponse(salinityLimited) = ...
    (tds(salinityLimited) - 2*50 + 300).* ...
    (300 - tds(salinityLimited))./(300 - 50)^2;
tdsResponse(tds > 300) = 0;

response = max(temperatureResponse, 0).*max(pHResponse, 0).* ...
    max(tdsResponse, 0);
end

function response = asymmetricParabolicResponse(input, limits)
lower = limits(1);
optimum = limits(2);
upper = limits(3);
response = zeros(size(input));
belowOptimum = input > lower & input < optimum;
response(belowOptimum) = -(input(belowOptimum) - lower).* ...
    (input(belowOptimum) + lower - 2*optimum)./(optimum - lower)^2;
aboveOptimum = input >= optimum & input < upper;
response(aboveOptimum) = -(input(aboveOptimum) - upper).* ...
    (input(aboveOptimum) + upper - 2*optimum)./(optimum - upper)^2;
end