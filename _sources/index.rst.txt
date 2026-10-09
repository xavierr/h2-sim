.. raw:: html

   <div class="home-hero">
    <div class="hero-copy">
     <p class="eyebrow">OPEN TOOLS FOR UNDERGROUND HYDROGEN STORAGE</p>
     <h1>Understand recovery.<br>Explore what changes it.</h1>
     <p class="hero-description">Quantify hydrogen losses and compare storage strategies with open MATLAB–MRST workflows for porous reservoirs.</p>
     <div class="hero-actions"><a class="primary-button" href="installation.html">Get started ↗</a><a class="secondary-button" href="examples.html">Explore notebooks →</a></div>
    </div>
    <figure class="hero-figure"><img src="_static/reservoir.svg" alt="Conceptual porous reservoir: hydrogen injection, a gas plume, and local microbial reactions." width="520" height="420"><figcaption>Understand the physics. Inspect the equations.</figcaption></figure>
   </div>

Explore what H2sim can do
=========================

From reservoir cycling to reactive losses, connect a storage question to
simulation results. Select a panel to discover the capabilities and examples.

.. raw:: html

   <div class="assessment-grid">
    <details class="assessment-card">
     <summary>
      <span class="assessment-visual"><img src="_static/home-storage.png" alt="Conceptual illustration of hydrogen storage and cycling in a porous reservoir." loading="lazy"></span>
      <span class="assessment-number">01 · FLOW &amp; CYCLING</span>
      <span class="assessment-title">Storage performance</span>
      <span class="assessment-question">How does hydrogen move through a storage cycle?</span>
      <span class="assessment-toggle"><span class="when-closed">Explore capabilities</span><span class="when-open">Close overview</span><span aria-hidden="true"> +</span></span>
     </summary>
     <div class="assessment-description">
      <p>Simulate hydrogen injection, storage, and withdrawal in saline aquifers and depleted gas reservoirs. Examine plume migration, dissolution, component production rates, and well pressure across repeated cycles. Adapt reservoir properties and operating schedules to compare storage scenarios.</p>
      <a href="store.html">Explore aquifer storage →</a>
      <a href="biochem.html">Explore depleted reservoirs →</a>
     </div>
    </details>
    <details class="assessment-card">
     <summary>
      <span class="assessment-visual"><img src="_static/home-microbial.png" alt="Conceptual illustration of microbes consuming hydrogen in brine-filled sandstone pores." loading="lazy"></span>
      <span class="assessment-number">02 · MICROBIAL REACTIONS</span>
      <span class="assessment-title">Hydrogen losses &amp; gas quality</span>
      <span class="assessment-question">Which reactions consume the stored hydrogen?</span>
      <span class="assessment-toggle"><span class="when-closed">Explore capabilities</span><span class="when-open">Close overview</span><span aria-hidden="true"> +</span></span>
     </summary>
     <div class="assessment-description">
      <p>Compare bacterial and abiotic storage scenarios, track microbial growth, and quantify hydrogen consumption by methanogenesis, acetogenesis, and sulfate reduction. Inspect how these reactions change hydrogen, carbon dioxide, and methane well streams, then test sensitivity to reaction parameters and storage duration.</p>
      <a href="biochem.html">Explore biochemical storage →</a>
      <a href="notebooks/bacterial_storage_2d.html">Open the cycling notebook →</a>
     </div>
    </details>
    <details class="assessment-card">
     <summary>
      <span class="assessment-visual"><img src="_static/home-chemistry.png" alt="Conceptual illustration of aqueous chemistry and mineral reactions in porous sandstone." loading="lazy"></span>
      <span class="assessment-number">03 · REACTIVE CHEMISTRY</span>
      <span class="assessment-title">Chemistry &amp; operating choices</span>
      <span class="assessment-question">When does chemistry change the storage response?</span>
      <span class="assessment-toggle"><span class="when-closed">Explore capabilities</span><span class="when-open">Close overview</span><span aria-hidden="true"> +</span></span>
     </summary>
     <div class="assessment-description">
      <p>Couple reservoir flow and microbial kinetics with PHREEQC to investigate aqueous speciation, pH, and mineral reactions. Use compositional or multirate coupling to study chemistry feedback through injection, idle, and production periods. Compare scenarios to inform operating choices, with thermodynamic models supplying the gas–liquid partitioning.</p>
      <a href="phreeqc.html">Explore PHREEQC workflows →</a>
      <a href="thermodynamics.html">Explore phase behavior →</a>
     </div>
    </details>
   </div>

Build the assessment you need
------------------------------

Start with hydrogen–brine flow in **h2-store**, add compositional transport
and microbial reactions in **h2-biochem**, and include **PHREEQC** when
aqueous and mineral chemistry matter. Thermodynamic models and PVT tables
support the flow calculations. Open equations and saved simulation outputs
let reservoir engineers and researchers inspect and adapt each workflow.

Start small. Build understanding.
=================================

.. grid:: 1 1 2 2
   :gutter: 3

   .. grid-item-card:: Installation that fits your setup
      :link: dependencies
      :link-type: doc

      Required dependencies, optional PHREEQC, and serial execution without
      Parallel Computing Toolbox.

      +++
      Check dependencies →

   .. grid-item-card:: Guided MATLAB notebooks
      :link: examples
      :link-type: doc

      A 2D aquifer cycle, microbial transport, thermodynamics, and PHREEQC validation.

      +++
      Browse examples →

.. admonition:: A serial path for every selected example
   :class: platform-note

   Parallel Computing Toolbox is optional. PHREEQC COM examples require
   Windows MATLAB and a registered IPhreeqcCOM server; MRST-only examples do not.

Research & community
====================

.. raw:: html

   <aside class="release-announcement" aria-label="Upcoming H2sim release">
    <div class="release-logos"><img class="release-h2sim-logo" src="_static/h2sim-logo.svg" alt="H2sim" width="190" height="54"><span>Built on</span><a href="https://www.sintef.no/projectweb/mrst/" aria-label="MRST — MATLAB Reservoir Simulation Toolbox"><img class="release-mrst-logo" src="_static/mrst-logo.png" alt="MRST — MATLAB Reservoir Simulation Toolbox" width="992" height="116"></a></div>
    <p class="eyebrow">COMING SOON · FIRST RELEASE</p>
    <h2>H2sim is getting ready to launch.</h2>
    <p>Built on the <strong>MATLAB Reservoir Simulation Toolbox (MRST)</strong>, H2sim brings together open tools for hydrogen storage in porous reservoirs: hydrogen–brine flow, microbial reactions, PHREEQC coupling, and thermodynamic models, with selected examples and guided MATLAB notebooks.</p>
    <p>Our first release is in preparation. The release date and download details will be announced here.</p>
    <a class="primary-button" href="https://github.com/xavierr/h2-sim">Follow H2sim on GitHub ↗</a>
   </aside>


.. grid:: 1 1 2 2
   :gutter: 3

   .. grid-item-card:: Publications & presentations
      :link: publications
      :link-type: doc

      Papers, conference talks, and posters behind the simulation workflows.

      +++
      Explore the research →

   .. grid-item-card:: Coming up · HydroGEMM 2026
      :link: events
      :link-type: doc

      **16–18 November · La Rochelle, France**

      Elyes Ahmed’s invited talk on MRST, PHREEQC, and reactive hydrogen storage.

      +++
      Talks, releases & events →

.. toctree::
   :hidden:
   :caption: Getting started

   installation
   dependencies
   quickstart

.. toctree::
   :hidden:
   :caption: Models

   store
   biochem
   phreeqc

.. toctree::
   :hidden:
   :caption: Learn and explore

   thermodynamics
   examples

.. toctree::
   :hidden:
   :caption: Research & community

   publications
   events

.. toctree::
   :hidden:
   :caption: Reference

   reference
   bibliography
   funding
