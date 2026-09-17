============
Installation
============

The source code of H2sym can be installed using git. In this way, you can easily keep track of the
last developments.

First, open a terminal and navigate to the directory where you would like to install H2sym. Then, clone the repository
using the following command, which will include all the dependencies as submodules

.. code-block:: shell

   git clone --recurse-submodules https://github.com/xavierr/h2-sim.git

Then, run :code:`startupH2sim.m` 

.. _MRST: https://www.sintef.no/Projectweb/MRST/


Update existing installation
============================

In the case where we alread have installed H2sym and you want to update to the latest version. As usual in git, you
will do that by running

.. code-block:: shell

   git pull

In addition to that, the dependencies that are given through git submodules. The are not updated often but, if it is the
case, you will need to run in addition to the previous command,

.. code-block:: shell

   git submodule update --recursive

Dependencies
============

The depencies are included as submodules so that no special installation for those is needed.

The package is built upon

- **MATLAB**: Version R2021a or newer
- **MRST**: MATLAB Reservoir Simulation Toolbox (2023b or newer), see `MRST`_
- **MRST Modules**:

  - ``compositional``
  - ``ad-blackoil``
  - ``ad-core``
  - ``ad-props``

