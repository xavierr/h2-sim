Code reference
==============

Browse the `source repository <https://github.com/xavierr/h2-sim>`_, or use
MATLAB's ``help`` and ``edit`` commands for local function documentation.

Entry points
------------

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - Function
     - Purpose
   * - ``startupH2sim``
     - Initialize code and dependencies
   * - ``exampleBlackOilAquifer2D``
     - 2D structural trapping and five storage cycles
   * - ``simple1DBacterialExampleMetacet``
     - MRST compositional microbial example
   * - ``setupH2StorageExampleWithSRB_benchmark``
     - Configure a multi-reaction benchmark
   * - ``exampleCompositionalBacterial1DPhreeqcMinerals``
     - Original MRST and three PHREEQC mineral configurations
   * - ``exampleSequentialBiochemistryPhreeqc1D``
     - Coarse-flow / local-reaction multirate example
   * - ``exampleCompositionalBacterial2DThreeReactions``
     - Existing 2D reservoir case with MET, ACE, SRB, packed cycles, and selectable scenarios
   * - ``runPhreeqcValidation``
     - Compositional PHREEQC injection benchmark against an external reference
   * - ``exampleH2BrineSolubilityStudy``
     - Phase-partition and salting-out comparisons
   * - ``setupOptimizedLinearSolver``
     - AMGCL CPR including bacterial and tracer unknowns; configurable native threads
   * - ``loadH2simResults``
     - Load result MAT files; stage Windows UNC files locally for HDF5
   * - ``activateH2simAMGCL``
     - Activate an optional locally compiled MRST gateway
   * - ``profilePhreeqcWorkflows``
     - Paired profiling and state comparisons

Output and replotting
---------------------

The aquifer, compositional 2D, validation, and short chemistry examples save their setup and accepted states with MRST
packed simulation output. ``packed_setup.mat`` sits beside the saved states,
well solutions, and timestep reports. Changed inputs select a separate result
directory. Returned results contain the setup and state histories;
``packedProblem`` identifies the black-oil and multirate output locations.
The mineral comparison returns one result per model configuration.

.. code-block:: matlab

   plotBlackOilAquifer2D(result);
   plotH2BrineSolubilityStudy(study);
   plotCompositionalBacterial2DThreeReactions(results);
   plotCompositionalPhreeqcValidation(validation);
   plotPhreeqcMineralComparison(comparison);
   plotSequentialBiochemistryPhreeqc(summary);
   plotPhreeqcBackendComparison(scenarios);

The compositional 2D MAT results retain physical fields while omitting
transient derivative/state-function caches. Complete restart states remain
in the packed directory. A completed MAT result with matching packed setup
can be reused without rereading every restart cache. Use ``loadH2simResults``
when loading large result files through a Windows WSL/network path; it
stages a temporary local copy for HDF5 and removes it after loading.

Plotting functions preserve values and export vector PDF and PNG figures.
``outputDirectory=''`` disables export. Dimensionless distance uses actual
domain edges, rather than stretching cell centers to zero and one.

Hydrogen-consumption percentages refer to microbial reaction consumption,
with the injection denominator stated by the figure. They are distinct from
unrecovered gas or a production recovery factor.

Tests
-----

.. code-block:: matlab

   checks = runtests('h2-biochem/tests');
   assertSuccess(checks);

COM integration tests require Windows and a registered IPhreeqcCOM server.
A skipped test does not validate a backend on that platform. The serial
fallback tests cover both available and missing optional toolbox checks.

See :doc:`phreeqc` for the current conservation-audit limitations.
