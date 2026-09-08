# MRST Bio-Chemistry Module for Hydrogen Storage Simulation

A comprehensive MATLAB Reservoir Simulation Toolbox (MRST) module for simulating hydrogen storage in depleted reservoirs with bio-chemical reactions and compositional modeling.

## Overview

This module extends MRST's capabilities by integrating a bio-chemistry model with the compositional simulator, specifically designed for hydrogen storage applications. It implements the Soreide-Whitson (SW) equation of state fitted to experimental data and enables simulation of microbial activity affecting hydrogen storage operations.

### Key Features

- **Compositional Modeling**: Full compositional simulation for H₂O-H₂-CO₂-CH₄-N₂ mixtures
- **Bio-Chemistry Integration**: Microbial growth and methanogenesis reactions
- **Equation of State**: Soreide-Whitson EoS fitted to experimental data
- **Bacterial Effects**: Simulation of hydrogen loss due to microbial activity

### Optional PHREEQC backends

`setupH2StorageExampleWithSRB_benchmark` supports exactly two PHREEQC
backends. Both require Windows, a registered `IPhreeqcCOM.Object` (or
configured `phreeqcComProgId`), and an explicit absolute
`phreeqcDatabaseFile` path to `PHREEQC_Modified.DAT`.

Set `phreeqcBackend='sequential-compositional-phreeqc'` with
`phreeqcTimestepCoupling=true` to run the post-convergence compositional
kinetics/chemistry split.

The COM split carries separate PHREEQC MET/ACE/SRB biomass
(`N0=1e9`, `Nmax=1e13` cells/kg water) and disables MRST's implicit microbial
reaction sources to avoid double counting; aqueous tracer transport (SO4, HS,
HCO3, Ca, Mg) remains active. `bactDiffusion` and `chemotaxisEffect` are
rejected for this backend: PHREEQC integrates only local per-cell kinetics
with no notion of spatial bacterial transport, so `nbact` cannot be diffused
or chemotaxis-moved independently of the biomass PHREEQC is actually growing.
Since its reaction source is also zero here, MRST does not assemble `nbact`'s
mass-balance equation at all for this backend (it would be a pure no-op every
step); `nbact` is still carried as a state field, and
`PsiGrowthRate`/`CarbonLimitedGrowthRate`/`BacterialMass` remain available as
diagnostic-only outputs. It maps the prescribed selected-output schema
back to tracers, minerals, and EOS inventories before reflashing. This
sequential coupling is not claimed to exactly reproduce any paper or external
benchmark.

Both backends reject a PHREEQC result before updating the state unless H, C, S,
Ca, Mg, and Fe are conserved in every cell. The returned state records the
`nc`-by-6 diagnostics `phreeqcElementBalanceInput`,
`phreeqcElementBalanceOutput`, `phreeqcElementBalanceAbsoluteResidual`,
`phreeqcElementBalanceNormalizedResidual`, and
`phreeqcElementBalancePass`; column names are in
`phreeqcElementBalanceElements`. Configure scalar or six-element tolerances
with `phreeqcElementBalanceAbsoluteTolerance` (default `1e-7` mol) and
`phreeqcElementBalanceRelativeTolerance` (default `1e-8`) in
`phreeqcCouplingOptions`.

The audit covers aqueous analytical totals (including the separate acetate
element), EOS H2/CO2/CH4/H2S, and all configured equilibrium minerals. Hydrogen
uses PHREEQC's system inventory with the fixed initial 1 kg solvent-water
baseline removed, avoiding subtraction of cell-scale solvent hydrogen while
retaining reaction- and hydrate-water changes. The compositional backend's
kinetic biomass is outside the reactive-element inventory because its PHREEQC
definition has `-formula H 0` and defines no C, S, Ca, Mg, or Fe storage;
MET/ACE/SRB kinetic amounts are reaction extents rather than stored products.
MRST `nbact` is likewise outside the equilibrium-only boundary.

`phreeqcBackend='sequential-h2biochem-phreeqc'` retains MRST's biochemical
sources and adds sequential PHREEQC equilibrium feedback.
It also requires a registered IPhreeqcCOM server and an absolute
`PHREEQC_Modified.DAT` path, but contains **no** PHREEQC `RATES` or
`KINETICS`. MRST's existing `state.nbact` Monod model remains the sole
reaction owner: bacterial growth, `BactConvertionRate`, and tracer reaction
sources remain active.

```matlab
[~, model, schedule, state0] = setupH2StorageExampleWithSRB_benchmark( ...
    'phreeqcBackend', 'sequential-h2biochem-phreeqc', ...
    'phreeqcTimestepCoupling', true, ...
    'phreeqcDatabaseFile', 'C:\PHREEQC\database\PHREEQC_Modified.DAT');
[wellSols, states, report] = simulateSequentialH2BiochemPhreeqc( ...
    state0, model, schedule);
```

Run this backend with `simulateSequentialH2BiochemPhreeqc`, rather than
directly with `simulateScheduleAD`. For every nominal schedule timestep, the
wrapper repeats the MRST solve from the fixed timestep-start state, then
equilibrates the already-reacted full component, tracer, gas, and mineral
inventories in PHREEQC. The relaxed preceding PHREEQC pH, DIC/CO2, and sulfate
snapshot is supplied only to MRST's Monod kinetic substrate evaluation; it is
not an MRST nonlinear initial guess or an accumulation state. Iteration stops
only when that chemistry feedback and the MRST-reported reaction extent
(`h2ConsumptionRate * dt`) meet configured tolerances, otherwise it raises a
nonconvergence error. It is therefore a same-timestep outer Picard scheme, not
a lagged post-step chemistry split.

### Biochemical Reaction Model

The module simulates the methanogenesis reaction:
\[
4\text{H}_2 + \text{CO}_2 \longrightarrow \text{CH}_4 + 2\text{H}_2\text{O} + \text{energy}
\]

For detailed methodology and validation, see our publication:
[**Numerical Modeling of Bio-Reactive Transport During Underground Hydrogen Storage**](https://www.sciencedirect.com/science/article/pii/S0360319925039473)

## Installation

### Prerequisites

- **MATLAB**: Version R2021a or newer
- **MRST**: MATLAB Reservoir Simulation Toolbox (2023b or newer)
- **Required MRST Modules**:
  - `compositional`
  - `ad-blackoil` 
  - `ad-core`
  - `ad-props`
  - `h2store`
  -`biochemistry`
