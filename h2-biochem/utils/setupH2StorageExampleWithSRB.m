function [biochemFluid, model, schedule, state0] = setupH2StorageExampleWithSRB(varargin)
% Set up a 1D UHS case with methanogenesis, acetogenesis, and sulfate reduction.
%
% OPTIONAL PARAMETERS (property/value pairs):
%   rate                - 'medrate' (default) or 'highrate'
%   bacteriamodel       - Enable bacterial growth/decay (default: true)
%   bactDiffusion       - Enable microbial diffusion (default: false)
%   chemotaxisEffect    - Enable bacterial chemotaxis (default: false)
%   molecularDiffusion  - Enable molecular diffusion (default: false)
%   molecularDispersion - Enable mechanical dispersion (default: false)
%   bioClogging         - Enable bio‑clogging (default: false)
%   nbact0              - Initial bacterial concentration (default: based on rate)
%   initial_SO4         - Initial sulfate molality (mol/kg, default: 0.1)
%   pH                  - pH of the system (default: 7.2)

require ad-props compositional deckformat h2-biochem

% --- Parse optional parameters ---
opt = struct('rate', 'highrate', ...           % 'medrate' or 'highrate'
    'bacteriamodel', true, ...
    'bactDiffusion', false, ...
    'chemotaxisEffect', false, ...
    'molecularDiffusion', false, ...
    'molecularDispersion', false, ...
    'bioClogging', false, ...
    'nbact0', [], ...                         % Will be set based on rate
    'initial_SO4', 0.25, ...
    'pH', 7.2);
opt = merge_options(opt, varargin{:});

% --- Set case-specific parameters ---
switch lower(opt.rate)
    case 'highrate'
        % HIGH RATE CASE (from paper Table 5)
        mu_MET = 4.1 / 86400;          % 4.745e-5 s^-1
        mu_ACE = 1.9 / 86400;          % 2.199e-5 s^-1
        mu_SRB = 5.5 / 86400;          % 6.366e-5 s^-1
        KD_MET = 9e-6 / 55.5;          % 1.62e-7 mol/mol
        KD_ACE = 2.5e-6 / 55.5;        % 4.50e-8 mol/mol
        KD_SRB = 2.9e-6 / 55.5;        % 5.22e-8 mol/mol
        nbact0_default =  60.6;          % High initial biomass
        fprintf('=== HIGH RATE CASE ===\n');
    otherwise
        % MEDRATE CASE (default)
        mu_MET = 1.109 / 86400;        % 1.284e-5 s^-1
        mu_ACE = 0.872 / 86400;        % 1.009e-5 s^-1
        mu_SRB = 1.048 / 86400;        % 1.213e-5 s^-1
        KD_MET = 9e-6 / 55.5;          % 1.62e-7 mol/mol
        KD_ACE = 2.5e-6 / 55.5;        % 4.50e-8 mol/mol
        KD_SRB = 2.9e-6 / 55.5;        % 5.22e-8 mol/mol
        nbact0_default =  60.6;          % Moderate initial biomass
        fprintf('=== MODERATE RATE CASE ===\n');
end

% Set nbact0 if not provided
if isempty(opt.nbact0)
    opt.nbact0 = nbact0_default;
end

%% 1. Grid and rock – 1D horizontal, 50 cells
G = cartGrid([50, 1, 1], [50, 1, 1]);
G = computeGeometry(G);
rock = makeRock(G, 100 * milli * darcy, 0.20);

%% 2. Fluid components
volatileNames   = {'Water', 'Hydrogen', 'CarbonDioxide', 'Methane', 'HydrogenSulfide', 'AceticAcid'};
volatileSymbols = {'H2O',   'H2',       'CO2',           'C1',      'H2S',             'CH3COOH'};
compFluid = TableCompositionalMixture(volatileNames, volatileSymbols);

%% 3. Biochemical reactions
reactNames = {'MethanogenicArchae', 'AcetogenicBacteria', 'SulfateReducingBacteria'};
biomassNames = {'bactM', 'bactA', 'bactS'};
biochemFluid = TableBioChemMixture(reactNames, biomassNames);

% ----- OVERRIDE KINETIC PARAMETERS BASED ON CASE -----
% Methanogens
idxM = strcmp(biochemFluid.metabolicReaction, 'MethanogenicArchae');
if any(idxM)
    biochemFluid.Psigrowthmax(idxM) = mu_MET;
    biochemFluid.alphaH2(idxM)      = KD_MET;
    biochemFluid.alphasub(idxM)     = 230e-6 / 55.5;
    biochemFluid.bbact(idxM)        = 0.01 * biochemFluid.Psigrowthmax(idxM);
    biochemFluid.Y_H2(idxM)         = 0.03.*3.333e12; %5.63789E+11;
    biochemFluid.nbactMax(idxM)     = 1.0e10;
end

% Acetogens
idxA = strcmp(biochemFluid.metabolicReaction, 'AcetogenicBacteria');
if any(idxA)
    biochemFluid.Psigrowthmax(idxA) = mu_ACE;
    biochemFluid.alphaH2(idxA)      = KD_ACE;
    biochemFluid.alphasub(idxA)     = 115.5e-6 / 55.5;
    biochemFluid.bbact(idxA)        = 0.01 * biochemFluid.Psigrowthmax(idxA);
    biochemFluid.Y_H2(idxA)         = 0.07.*3.333e12;%1.31551E+12;
    biochemFluid.nbactMax(idxA)     = 1.0e10;

end

% SRB
idxS = strcmp(biochemFluid.metabolicReaction, 'SulfateReducingBacteria');
if any(idxS)
    biochemFluid.gamrH2(idxS)  = -4;                    % FIXED
    biochemFluid.gamrsub(idxS) = -1;
    biochemFluid.gampH2O(idxS) =  4;
    biochemFluid.gamp2(idxS)   =  1;
    biochemFluid.Psigrowthmax(idxS) = mu_SRB;
    biochemFluid.alphaH2(idxS)      = KD_SRB;
    biochemFluid.alphasub(idxS)     = 2751.5e-6 / 55.5;
    biochemFluid.bbact(idxS)        = 0.01 * biochemFluid.Psigrowthmax(idxS);
    biochemFluid.Y_H2(idxS)         = 0.08*3.333e12;%1.50344E+12;
    biochemFluid.nbactMax(idxS)     = 1.0e9;
end

%% 4. EOS
eos = SoreideWhitsonEos(G, compFluid, ...
    'msalt', 0, ...
    'pH', opt.pH, ...
    'initial_NaCl', 0, ...
    'initial_SO4', opt.initial_SO4, ...
    'rho_water', 1000);

%% 5. Fluid properties (base transport)
%% 5. Fluid properties – TABLE-BASED (matching repository)
% Repository uses table-based relative permeability, not Brooks-Corey.
Sw_table = [0, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, ...
    0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 0.999];
krW_table = [0, 0, 0.002, 0.01, 0.02, 0.033, 0.049, 0.066, 0.09, 0.119, 0.15, ...
    0.186, 0.227, 0.277, 0.33, 0.39, 0.462, 0.54, 0.62, 0.71, 0.8, 0.9, 1];
Sg_table = [0.05, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, ...
    0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 1];
krG_table = [0, 0.013, 0.026, 0.04, 0.058, 0.078, 0.1, 0.126, 0.156, 0.187, ...
    0.222, 0.26, 0.3, 0.348, 0.4, 0.45, 0.505, 0.562, 0.62, 0.68, 0.74, 1];

fluid = initSimpleADIFluid('phases', 'OG', ...
    'mu',  [1.3059*centi*poise, 0.01763*centi*poise], ...
    'rho', [999.7, 1.2243] .* kilogram/meter^3, ...
    'pRef', 150*barsa, ...
    'c',   [5e-5/barsa, 0.0067/barsa], ...
    'n',   [2, 2], ...
    'smin',[0.2, 0.05]);

% Override with table-based relative permeability (repository style)
fluid.krO = @(sw) interpTable(Sw_table, krW_table, sw);
fluid.krG = @(sg) interpTable(Sg_table, krG_table, sg);

% No capillary pressure (repository uses none)
fluid.pcOG = @(sg) 0;

%% 6. Assemble model
backend = DiagonalAutoDiffBackend('modifyOperators', true);
model = BiochemistryModel(G, rock, fluid, compFluid, biochemFluid, ...
    true, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', opt.bacteriamodel, ...
    'bactDiffusion', opt.bactDiffusion, ...
    'chemotaxisEffect', opt.chemotaxisEffect, ...
    'molecularDiffusion', opt.molecularDiffusion, ...
    'molecularDispersion', opt.molecularDispersion, ...
    'liquidPhase', 'O', 'vaporPhase', 'G');
model.EOSModel = eos;
model.OutputStateFunctions{end+1} = 'ComponentPhaseDensity';

% Bio-clogging
nBioReactions = biochemFluid.nbioreact;
if opt.bioClogging && opt.bacteriamodel
    nc = [180, 180, 180];
    cp = [0.5, 0.5, 0.5];
    nbact0 = opt.nbact0 * ones(1, nBioReactions);
    model = setupBioCloggingModel(model, nbact0, nc, cp, true);
else
    nc = [180, 180, 180];
    cp = [0.0, 0.0, 0.0];
    nbact0 = opt.nbact0 * ones(1, nBioReactions);
    model = setupBioCloggingModel(model, nbact0, nc, cp, false);
end

%% 7. Initial state
p0 = 150 * barsa;
T0 = 273.15 + 60;
Sw = 0.5508;
Sg = 1 - Sw;

z_gas = zeros(1, compFluid.getNumberOfComponents());
z_gas(strcmp(compFluid.names, 'C1')) = 1.0;

x_H2O = 1.0 - 4e-5 - 0.015;
x_CH4 = 0.015;
x_CO2 = 4e-5;

z_water = zeros(1, compFluid.getNumberOfComponents());
z_water(strcmp(compFluid.names, 'H2O')) = x_H2O;
z_water(strcmp(compFluid.names, 'C1')) = x_CH4;
z_water(strcmp(compFluid.names, 'CO2')) = x_CO2;

z0_overall = Sw * z_water + Sg * z_gas;
% Paper Table 6: Overall composition is 90% H2O and 10% CH4
% z0_overall = zeros(1, compFluid.getNumberOfComponents());
% z0_overall(strcmp(compFluid.names, 'H2O')) = 0.90;
% z0_overall(strcmp(compFluid.names, 'C1'))  = 0.10;
% All other components (H2, CO2, H2S, CH3COOH) are zero initially
ncell = G.cells.num;
z0 = repmat(z0_overall, ncell, 1);

if opt.bacteriamodel
    nbact0 = opt.nbact0 * ones(1, nBioReactions);
    state0 = initCompositionalStateBacteria(model, p0, T0, [Sw, Sg], z0, nbact0, eos);
else
    state0 = initCompositionalState(G, p0, T0, [Sw, Sg], z0, eos);
end

if isa(model.EOSModel, 'SoreideWhitsonEos') && model.sulfateReduction
    state0.tracerSO4 = repmat(opt.initial_SO4 * model.EOSModel.rho_water, ncell, 1);
    state0.tracerHS  = zeros(ncell, 1);
    state0.h2sDissolvedLag = zeros(ncell, 1);
end

%% 8. Wells
%% Storage Scenario
totTime = 250*day;
steps   = 125;
dt = totTime/steps;
dt = repmat(dt, steps, 1);
schedule = simpleSchedule(dt);
schedule.step.control(26:100)=2; % Control Value to switch to a new schedule
schedule.step.control(101:end)=3; % Control Value to switch to a new schedule
[L, x, y, Z_L, Z_V] = standaloneFlash(p0, T0, [0,1,0,0,0,0], eos);

Bg = 101325/298.15*Z_V*T0/p0;
pv = sum(G.cells.volumes.*rock.poro); % Calculating the pore volume of each cell
rate = 0.005*pv*meter^3/day/Bg; % Surface Rate
W1 = [];
W2 = [];
W3 = [];
tmp = cell(4,1);
schedule.control=struct('W',tmp,'bc',tmp,'src',tmp);
%Injection
W1 = verticalWell(W1, G, rock, 1, 1, 1,...
                'Type', 'rate', 'Val', rate, ...
                'Name', 'Injector','comp_i',[0 1],'sign',+1,'radius',0.1);
W1(1).components = [0, 0.95, 0.05, 0, 0, 0];
%Storage
W2 = verticalWell(W2, G, rock, 1, 1, 1,...
                'Type', 'rate', 'Val', 0, ...
                'Name', 'Shut-in','comp_i',[0 1],'sign',+1,'radius',0.1);
W2(1).components = [0, 1, 0, 0, 0, 0];
%Production
W3 = verticalWell(W3, G, rock, 1, 1, 1,...
                'Type', 'grat', 'Val', -rate, ...
                'Name', 'Producer','comp_i',[0.5 0.5],'sign',-1,'radius',0.1);
W3(1).components = [0, 1, 0, 0, 0, 0];
W3.lims.bhp = p0;

schedule.control(1).W=W1;
schedule.control(2).W=W2;
schedule.control(3).W=W3;

%% 10. Inform user
fprintf('Case: %s\n', opt.rate);
fprintf('  MET μmax = %.2e s^-1 (%.2f day^-1)\n', mu_MET, mu_MET*86400);
fprintf('  ACE μmax = %.2e s^-1 (%.2f day^-1)\n', mu_ACE, mu_ACE*86400);
fprintf('  SRB μmax = %.2e s^-1 (%.2f day^-1)\n', mu_SRB, mu_SRB*86400);
fprintf('  nbact0   = %.1e cells/m3\n', opt.nbact0);
fprintf('  Injection: 95%% H2 + 5%% CO2 (unchanged)\n');
end