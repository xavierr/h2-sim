function [biochemFluid, model, schedule, state0] = setupH2StorageExample(varargin)
% Set up a 1D UHS case reproducing the simple storage cycle from
% Shojae et al. (2025) – "New flow simulation framework ..."
%
% Grid: 50 x 1 x 1, 50 m length, 1 m width, 1 m height.
% Initial: 150 bar, 60 °C, Sw=0.5508, gas: 90% H2, 10% CH4.
% Schedule: inject H2 (0.565 kg/d) for 50 d, shut-in 150 d, produce (rate with BHP limit 150 bar) for 50 d.
%
% OPTIONAL PARAMETERS (property/value pairs):
%   bacteriamodel       - Enable bacterial growth/decay (default: false)
%   bactDiffusion       - Enable microbial diffusion (default: false)
%   chemotaxisEffect    - Enable bacterial chemotaxis (default: false)
%   molecularDiffusion  - Enable molecular diffusion (default: false)
%   molecularDispersion - Enable mechanical dispersion (default: false)
%   bioClogging         - Enable bio‑clogging (default: false)
%   nbact0              - Initial bacterial concentration (default: 1e-6)
%   ncycles             - Not used (only one cycle)
%
% RETURNS:
%   biochemFluid - Biochemical reaction definition (with stoichiometry & kinetics)
%   model        - BiochemistryModel
%   schedule     - Cyclic schedule (injection, shut-in, production)
%   state0       - Initial state

    require ad-props compositional deckformat h2-biochem

    % --- Parse optional parameters ---
    opt = struct('bacteriamodel', true, ...
                 'bactDiffusion', false, ...
                 'chemotaxisEffect', false, ...
                 'molecularDiffusion', false, ...
                 'molecularDispersion', false, ...
                 'bioClogging', false, ...
                 'nbact0', 10e0, ...
                 'ncycles', 1);
    opt = merge_options(opt, varargin{:});

    %% 1. Grid and rock – 1D horizontal, 50 cells
    G = cartGrid([50, 1, 1], [50, 1, 1]);
    G = computeGeometry(G);
    rock = makeRock(G, 100 * milli * darcy, 0.20);

    %% 2. Fluid components and EOS (5 components)
    compFluid = TableCompositionalMixture( ...
        {'Water', 'Hydrogen', 'CarbonDioxide', 'Methane', 'AceticAcid'}, ...
        {'H2O', 'H2', 'CO2', 'C1', 'CH3COOH'});
    biochemFluid = TableBioChemMixture({'MethanogenicArchae', 'AcetogenicBacteria'}, ...
                                       {'bactM', 'bactA'});

    % ----- OVERRIDE KINETIC PARAMETERS TO MATCH PAPER (MODERATE RATE) -----
    % Keep original database intact; override only the fields we need.
    % Values from Shojae et al. (2025) Table 5, converted to MRST units.

    % 1. MethanogenicArchae
    idxM = strcmp(biochemFluid.metabolicReaction, 'MethanogenicArchae');
    if any(idxM)
        biochemFluid.Psigrowthmax(idxM) = 1.109 / 86400;      % 1.284e-5 s^-1
        biochemFluid.alphaH2(idxM)      = 10e-6 / 55.5;       % 1.8e-7 (mol/mol)
        biochemFluid.alphasub(idxM)     = 230e-6 / 55.5;      % 4.14e-6 (mol/mol)
        biochemFluid.bbact(idxM)        = 0.01 * biochemFluid.Psigrowthmax(idxM); % 1.284e-7
        biochemFluid.Y_H2(idxM)         = 3.0e10;             % cells/mol H2 (scaled)
        %biochemFluid.nbactMax(idxM)     = 1.06e8;             % cells/m^3
    end

    % 2. AcetogenicBacteria
    idxA = strcmp(biochemFluid.metabolicReaction, 'AcetogenicBacteria');
    if any(idxA)
        biochemFluid.Psigrowthmax(idxA) = 0.872 / 86400;      % 1.009e-5 s^-1
        biochemFluid.alphaH2(idxA)      = 2.5e-6 / 55.5;      % 4.5e-8 (mol/mol)
        biochemFluid.alphasub(idxA)     = 115e-6 / 55.5;      % 2.07e-6 (mol/mol)
        biochemFluid.bbact(idxA)        = 0.01 * biochemFluid.Psigrowthmax(idxA); % 1.009e-7
        biochemFluid.Y_H2(idxA)         = 7.0e10;             % cells/mol H2 (scaled)
        %biochemFluid.nbactMax(idxA)     = 1.06e8;             % cells/m^3
    end

    % Optionally, override stoichiometry if you want to be explicit
    % (already correct in database, but we can set it anyway)
    biochemFluid.gamrH2  = [-4, -4];
    biochemFluid.gamrsub = [-1, -2];
    biochemFluid.gampH2O = [ 2,  2];
    biochemFluid.gamp2   = [ 1,  1];
    eos = SoreideWhitsonEos(G, compFluid, 'msalt', 0);  % fresh water

    % Simple fluid (used as a base for transport properties)
    fluid = initSimpleADIFluid('phases', 'OG', ...
        'mu',  [1.3059*centi*poise, 0.01763*centi*poise], ...
        'rho', [999.7, 1.2243] .* kilogram/meter^3, ...
        'pRef', 100*barsa, ...
        'c',   [5.0e-5/barsa, 1.0/barsa], ...
        'n',   [2, 2], ...
        'smin',[0.2, 0.05]);
    Pe = 0.1*barsa;
    fluid.pcOG = @(sg) Pe * max((1 - sg - 0.2) ./ (1 - 0.2), 1e-5).^(-1/2);

    %% 3. Assemble model
    backend = DiagonalAutoDiffBackend('modifyOperators', true);
    model = BiochemistryModel(G, rock, fluid, compFluid, biochemFluid, ...
        true, backend, ...   % explicit = true
        'water', false, 'oil', true, 'gas', true, ...
        'bacteriamodel', opt.bacteriamodel, ...
        'bactDiffusion', opt.bactDiffusion, ...
        'chemotaxisEffect', opt.chemotaxisEffect, ...
        'molecularDiffusion', opt.molecularDiffusion, ...
        'molecularDispersion', opt.molecularDispersion, ...
        'liquidPhase', 'O', 'vaporPhase', 'G');
    model.EOSModel = eos;

    % Bio‑clogging (if requested)
    if opt.bioClogging && opt.bacteriamodel
        nc = [180,180];               % critical bacteria concentration
        cp = [0.5, 0.5];               % scaling coefficient
        nbact0 = opt.nbact0 * ones(1, biochemFluid.nbioreact);
        model = setupBioCloggingModel(model, nbact0, nc, cp, true);
    else
        nc = [180,180];               % critical bacteria concentration
        cp = [0.0, 0.0];               % scaling coefficient
        nbact0 = opt.nbact0 * ones(1, biochemFluid.nbioreact);
        model = setupBioCloggingModel(model, nbact0, nc, cp, false);
    end

    %% 4. Initial state – based on paper's 1D case
    p0 = 150 * barsa;                     % initial pressure
    T0 = 273.15 + 60;                     % 60 °C
    % ---- MODIFIED: your requested initial saturations ----
    Sw = 0.90;                            % Water saturation (90% water, as you said)
    Sg = 1 - Sw;                          % Gas saturation = 10%

    % Gas phase: 100% CH4 (NO H2!)
    z_gas = [0, 0.0, 0.0, 1.0, 0.0];      % H2O=0, H2=0, CO2=0, CH4=1, CH3COOH=0

    % Water phase: pure water + dissolved CO2 (electron acceptor for bacteria)
    CO2_in_water = 7.929e-5;                % ~4e-5 mole fraction (adjust as needed)
    z_water = [1.0 - CO2_in_water, 0.0, CO2_in_water, 0.0, 0.0];

    % Overall mole fractions = Sw * z_water + Sg * z_gas
    z0_overall = Sw * z_water + Sg * z_gas;

    ncell = G.cells.num;
    z0 = repmat(z0_overall, ncell, 1);

    if opt.bacteriamodel
        nbact0 = opt.nbact0 * ones(1, biochemFluid.nbioreact);
        state0 = initCompositionalStateBacteria(model, p0, T0, [Sw, Sg], z0, nbact0, eos);
    else
        state0 = initCompositionalState(G, p0, T0, [Sw, Sg], z0, eos);
    end

    %% 5. Wells – injector at cell 1, producer at cell 50
    % Injection rate: 6.2811 sm3/day of H2 ≈ 0.565 kg/day (using rho_H2=0.0899 kg/m3)
    q_inj = 0.565 * kilogram/day;   % mass rate of H2
    q_prod = -0.565 * kilogram/day;

    % Well components: injection of pure H2
    comp_inj = [0.0, 0.95, 0.05, 0.0, 0.0];   % H2 only

    injCell = 1;
    prodCell = 1;

    % Injector
    W_inj = addWell([], G, rock, injCell, ...
        'Type', 'rate', ...
        'Val', q_inj, ...
        'Sign', 1, ...               % injection
        'components', comp_inj, ...
        'Comp_i', [0, 1], ...
        'Name', 'Injector');
    % Producer – rate control with BHP limit
    W_prod = addWell([], G, rock, prodCell, ...
        'Type', 'rate', ...
        'Val', q_prod, ...
        'Sign', -1, ...              % production
        'components', [0, 1, 0, 0, 0], ... % will be overwritten by reservoir composition
        'Comp_i', [0, 1], ...
        'Name', 'Producer');
    % Add BHP limit: minimum BHP for producer (if pressure drops below, switch to BHP)
    W_prod.lims.bhp = 150 * barsa;

    %% 6. Schedule: injection 50 d, shut-in 150 d, production 50 d
    injDays  = 50 * day;
    shutDays = 150 * day;
    prodDays = 50 * day;

    % Timesteps: 1 day during active periods, 5 days during shut-in
    dt_inj  = 1 * day;
    dt_shut = 5 * day;
    dt_prod = 1 * day;

    % Number of time steps
    nstep_inj  = ceil(injDays / dt_inj);    % 50
    nstep_shut = ceil(shutDays / dt_shut);  % 30  ← fix here
    nstep_prod = ceil(prodDays / dt_prod);  % 50

    scheduleInj  = simpleSchedule(repmat(dt_inj,  nstep_inj,  1), 'W', W_inj);
    W_shut = W_inj;
    W_shut.val = 0;
    W_shut.sign = 0;
    W_shut.status = false;
    scheduleShut = simpleSchedule(repmat(dt_shut, nstep_shut, 1), 'W', W_shut);
    scheduleProd = simpleSchedule(repmat(dt_prod, nstep_prod, 1), 'W', W_prod);

    schedule = combineSchedules(scheduleInj, scheduleShut, scheduleProd, ...
        'makeConsistent', false);

    % Clean up well limits (if any)
    for i = 1:numel(schedule.control)
        for j = 1:numel(schedule.control(i).W)
            schedule.control(i).W(j).lims = [];
        end
    end
    % Re‑apply BHP limit to producer controls (only during production stage)
    for i = 3:numel(schedule.control)
        for j = 1:numel(schedule.control(i).W)
            if strcmp(schedule.control(i).W(j).name, 'Producer')
                schedule.control(i).W(j).lims.bhp = 150 * barsa;
            end
        end
    end
end