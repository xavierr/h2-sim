=============================================================
Hydrogen Storage Simulation in Aquifers using Black-Oil Model
=============================================================

This module provides simulation tools for modeling a Hydrogen and brine mixture within the
**ad-blackoil** module in the **MATLAB Reservoir Simulation Toolbox (MRST)**.

Our module includes tools for implementing the **Redlich-Kwong** (**RK**) equation of state (EoS)
and generating tabulated PVT data for H₂-brine systems. Additionally, it provides solubility tables
derived from **ePC-Saft** and **Henry's law** EoS for precise phase behavior calculations in
hydrogen storage simulations.

Overview
========

This module is developed to simulate the phase behavior and thermodynamic properties of hydrogen
stored in saline aquifers, with specific consideration of temperature, pressure, and salinity
effects on hydrogen solubility and fluid properties.

Features
========

- **Implementation of RK Equation of State**: The module includes an implementation of the
  Redlich-Kwong (RK) EoS for H₂-brine mixtures, allowing for quick thermodynamic calculations that
  capture real gas behavior and temperature effects.

- **Scripts for tabulating PVT Data for blackoil simulators**: Precomputed PVT tables facilitate
  efficient simulation of H₂-brine mixtures in MRST's ad-blackoil module. The resulting tables are
  also suitable for other blackoil simulators.

- **Solubility tables**: Estimate solubility from both ePC-Saft (data) and Henry-Setschnow
  correlation.

- **Correlations** for estimating the properties of water and H₂ mixtures for different pressure,
  temperature and salinity.

