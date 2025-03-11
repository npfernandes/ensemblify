.. _Usage:

💻 Usage
========

Ensemblify offers four main modules, all of which can be accessed either through the command line or from inside a Python script/Jupyter Notebook.

The `generation` module
-----------------------

With the `generation` module, you can generate conformational ensembles for your protein of interest.

Before generating an ensemble, you must create a parameters file either through the provided `parameters form <https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html>`_ or directly by editing the provided `parameters file template <docs/assets/parameters_template.yaml>`_. Check the :ref:`parameters file setup <Parameters File Setup>` section for more details.

To generate an ensemble, provide Ensemblify with the path to your parameters file.

Using the `ensemblify` command in a terminal:

.. code-block:: bash

    ensemblify generation -p parameters_file.yaml

Inside a Python script or Jupyter Notebook:

.. code-block:: python

    from ensemblify.generation import generate_ensemble
    generate_ensemble('parameters_file.yaml')

Check the interactive `Generation Module <examples/02_generation_module.ipynb>`_ notebook for detailed usage examples.

.. _Parameters File Setup:

Setting up your parameters file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An `.html form <https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html>`_ is provided to aid you in building your parameters file.

.. raw:: html

   <details>
   <summary>Parameters Form Preview</summary>

.. image:: ../assets/parameters_form_preview.svg

.. raw:: html

   </details>

If you prefer to create your own parameters file from scratch, a `template file <docs/assets/parameters_template.yaml>`_ is also provided.

The `conversion` module
-----------------------

With the `conversion` module, you can convert your generated .pdb structures into a .xtc trajectory file, enabling you to easily store and analyze your conformational ensemble.

To do this, provide the name for your created trajectory, the directory where the ensemble is stored and the directory where the trajectory file should be created.

Using the `ensemblify` command in a terminal:

.. code-block:: bash

    ensemblify conversion -j trajectory_name -e ensemble_dir -t trajectory_dir

Inside a Python script or Jupyter Notebook:

.. code-block:: python

    from ensemblify.conversion import ensemble2traj
    ensemble2traj('trajectory_name','ensemble_dir','trajectory_dir')

Check the interactive `Conversion Module <examples/03_conversion_module.ipynb>`_ notebook for detailed usage examples.

The `analysis` module
---------------------

With the `analysis` module, you can create an interactive graphical dashboard displaying structural information calculated from the conformational ensemble of your protein of interest.

To do this, provide your ensemble in trajectory format, your trajectory's topology file and the name you want to use for your protein in the graphical dashboard.

Using the `ensemblify` command in a terminal:

.. code-block:: bash

    ensemblify analysis -trj trajectory.xtc -top topology.pdb -tid trajectory_name

Inside a Python script or Jupyter Notebook:

.. code-block:: python

    from ensemblify.analysis import analyze_trajectory
    analyze_trajectory('trajectory.xtc','topology.pdb','trajectory_name')

Check the interactive `Analysis Module <examples/04_analysis_module.ipynb>`_ notebook for detailed usage examples.

The `reweighting` module
------------------------

With the `reweighting` module, you can use experimental SAXS data to reweigh your conformational ensemble following the Bayesian Maximum Entropy method [12]_.

To do this, provide your ensemble in trajectory format, your trajectory's topology file, the name you want to use for your protein in the graphical dashboard and your experimental SAXS data.

Using the `ensemblify` command in a terminal:

.. code-block:: bash

    ensemblify reweighting -trj trajectory.xtc -top topology.pdb -tid trajectory_name -exp exp_SAXS_data.dat

Inside a Python script or Jupyter Notebook:

.. code-block:: python

    from ensemblify.reweighting import reweight_ensemble
    reweight_ensemble('trajectory.xtc','topology.pdb','trajectory_name','exp_SAXS_data.dat')

Check the interactive `Reweighting Module <examples/05_reweighting_module.ipynb>`_ notebook for detailed usage examples. 

References
----------

.. [12] S. Bottaro , T. Bengsten and K. Lindorff-Larsen, "Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach," pp. 219-240, Feb. 2020. In: Z. Gáspári, (eds) *Structural Bioinformatics*, *Methods in Molecular Biology*, vol. 2112, Humana, New York, NY. `Link <https://doi.org/10.1007/978-1-0716-0270-6_15>`_ 
