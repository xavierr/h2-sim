Quick start
===========

After :doc:`installation`, choose one of the three model families.
All examples below use serial execution.

Hydrogen storage in saline aquifers
-----------------------------------

.. code-block:: matlab

   startupH2sim
   result = exampleBlackOilAquifer2D;
   lastState = result.states{end};
   gasSaturation = lastState.s(:, 2);
   plotBlackOilAquifer2D(result);

The example uses a 2D structural trap with hydrogen build-up and five storage
cycles. Gravity, capillarity, and dissolved hydrogen use the existing aquifer
model and input tables. See :doc:`store` and :doc:`thermodynamics`.

Hydrogen storage in depleted gas reservoirs
-------------------------------------------

.. code-block:: matlab

   results = exampleCompositionalBacterial2DThreeReactions( ...
       'numCycles', 5, 'scenarios', {'bacterial', 'abiotic'}, ...
       'molecularDiffusion', false, 'molecularDispersion', false);
   plotCompositionalBacterial2DThreeReactions(results);

Open the compositional notebook in :doc:`examples` to inspect the setup,
well controls, microbial populations, and transport plots.

Short PHREEQC mineral comparison
--------------------------------

On Windows with IPhreeqcCOM configured:

.. code-block:: matlab

   comparison = exampleCompositionalBacterial1DPhreeqcMinerals( ...
       'gridCells', 4, 'maximumFlowSteps', 2);
   pH = comparison(2).states{end}.phreeqcPH;

Redraw the completed states without solving again:

.. code-block:: matlab

   plotPhreeqcMineralComparison(comparison, ...
       'outputDirectory', fullfile(pwd, 'build', 'my-figures'));

These cases use different kinetic owners. Read :doc:`phreeqc` before
interpreting a backend comparison as a sensitivity to mineral chemistry.
