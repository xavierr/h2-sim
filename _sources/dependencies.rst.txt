Dependencies
============

The standard workflows run serially. Install optional components only when
your selected model needs them.

.. list-table:: What you need
   :header-rows: 1
   :widths: 28 36 36

   * - Component
     - Required for
     - Installation
   * - MATLAB
     - Core examples; R2021a or newer
     - Your MATLAB installation and license
   * - MRST
     - All simulation workflows
     - Pinned Git submodule; initialized by ``startupH2sim``
   * - PhreeqcMatlab
     - Included chemistry integration tools
     - Pinned Git submodule
   * - IPhreeqcCOM
     - The two sequential PHREEQC backends
     - Windows COM installation and registration
   * - ``h2_biogeochemistry.dat``
     - Selected PHREEQC examples
     - Bundled under ``h2-biochem/database``
   * - Parallel Computing Toolbox
     - Optional parallel pressure flashes
     - Optional; serial is the default
   * - AMGCL and a supported C++ compiler
     - Optional acceleration of large compositional cases
     - MRST linear-solver MEX gateway; no Parallel Computing Toolbox needed
   * - Internet access to NIST
     - Pure-component table-generation scripts
     - Needed when fetching data, not for the compact notebooks

Reference validation
--------------------

The PHREEQC validation case additionally uses a separate UGFACT checkout.
Pass its root directory as ``ugfactRoot`` to ``runPhreeqcValidation``.
All selected examples and the reference run use serial execution.

Without Parallel Computing Toolbox
----------------------------------

Leave ``phreeqcReflashParallel=false`` (the default). The selected models and
notebooks do not call ``parpool`` or require worker processes.

.. code-block:: matlab

   summary = exampleSequentialBiochemistryPhreeqc1D( ...
       'maximumFlowSteps', 2, 'phreeqcReflashParallel', false);

If parallel pressure flashes are explicitly requested but the toolbox or its
license is unavailable, H2sim reports the fallback and uses the serial loop.
Parallel execution is an optional optimization, not a different model.

Dependency revisions
--------------------

Git records the exact MRST and PhreeqcMatlab revisions used by the checkout:

.. code-block:: shell

   git submodule status
   git submodule update --init --recursive

H2sim's acetate component parameters live in its own helper; no local edit
to MRST's property table is required. Avoid recursively adding an entire MRST
checkout to MATLAB's path. Its startup and module manager select the needed paths.

PHREEQC platform requirements
-----------------------------

The COM workflows require MATLAB on Windows and a registered IPhreeqcCOM
server. MATLAB on Linux or macOS cannot activate this Windows COM interface. See
:doc:`installation` for the COM setup check.

Database source
---------------

``h2_biogeochemistry.dat`` is an unmodified copy of the database distributed
with `UGFACT <https://github.com/ahmadrezashojaee/UGFACT>`_, developed by
Ahmadreza Shojaee at Heriot-Watt University. H2sim changes only its filename.
The `exact source revision <https://github.com/ahmadrezashojaee/UGFACT/blob/1f8847ea1895cef52a6f004f46a81614a1d57989/examples/database/PHREEQC_Modified.DAT>`_ is recorded alongside the checksum in
``h2-biochem/database/README.md``. The original USGS PHREEQC database credits
are preserved. Cite `UGFACT's framework publication
<https://doi.org/10.1016/j.ijhydene.2025.150453>`_ when using this database.
