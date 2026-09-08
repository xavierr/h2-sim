function [inventory, elements] = computePhreeqcElementInventory(reservoirs)
% Sum audited PHREEQC-boundary elemental inventories in moles per cell.
%
% Required aqueous fields are total inorganic carbon (c4), acetate, sulfate
% (s6), sulfide (s2), calcium, magnesium, and Fe(II)/Fe(III). Gas and mineral
% fields are component/phase moles. The hydrogen field is reactive aqueous H;
% callers may replace the final H column with PHREEQC's solvent-normalized
% system total when aqueous speciation is known only inside PHREEQC.

elements = {'H', 'C', 'S', 'Ca', 'Mg', 'Fe'};
required = {'hydrogen', 'c4', 'acetate', 's6', 's2', 'ca', 'mg', ...
    'fe2', 'fe3', 'gasH2', 'gasCO2', 'gasCH4', 'gasH2S', ...
    'calcite', 'dolomite', 'anhydrite', 'gypsum', 'goethite', ...
    'pyrite', 'brucite', 'portlandite'};
for i = 1:numel(required)
    assert(isfield(reservoirs, required{i}), ...
        'PHREEQC element inventory is missing reservoir field "%s".', required{i});
end

nc = numel(reservoirs.hydrogen);
values = zeros(nc, numel(required));
for i = 1:numel(required)
    v = reservoirs.(required{i});
    v = v(:);
    assert(numel(v) == nc && all(isfinite(v)) && all(v >= -1e-12), ...
        'PHREEQC reservoir "%s" must contain %d finite non-negative moles.', ...
        required{i}, nc);
    values(:, i) = max(v, 0);
end
r = cell2struct(num2cell(values, 1), required, 2);

% Columns are H, C, S, Ca, Mg, Fe. Mineral rows explicitly encode:
% calcite CaCO3; dolomite CaMg(CO3)2; anhydrite CaSO4;
% gypsum CaSO4.2H2O; goethite FeOOH; pyrite FeS2;
% brucite Mg(OH)2; portlandite Ca(OH)2. Quartz is SiO2 and therefore has
% no atom of any audited element.
gas = [2, 0, 0, 0, 0, 0; ... % H2
       0, 1, 0, 0, 0, 0; ... % CO2
       4, 1, 0, 0, 0, 0; ... % CH4
       2, 0, 1, 0, 0, 0];    % H2S
mineral = [0, 1, 0, 1, 0, 0; ... % calcite
           0, 2, 0, 1, 1, 0; ... % dolomite
           0, 0, 1, 1, 0, 0; ... % anhydrite
           4, 0, 1, 1, 0, 0; ... % gypsum
           1, 0, 0, 0, 0, 1; ... % goethite
           0, 0, 2, 0, 0, 1; ... % pyrite
           2, 0, 0, 0, 1, 0; ... % brucite
           2, 0, 0, 1, 0, 0];    % portlandite

inventory = zeros(nc, numel(elements));
inventory(:, 1) = r.hydrogen;
inventory(:, 2) = r.c4 + 2*r.acetate;
inventory(:, 3) = r.s6 + r.s2;
inventory(:, 4) = r.ca;
inventory(:, 5) = r.mg;
inventory(:, 6) = r.fe2 + r.fe3;
inventory = inventory + ...
    [r.gasH2, r.gasCO2, r.gasCH4, r.gasH2S]*gas + ...
    [r.calcite, r.dolomite, r.anhydrite, r.gypsum, r.goethite, ...
     r.pyrite, r.brucite, r.portlandite]*mineral;
end
