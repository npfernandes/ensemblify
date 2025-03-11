🌐 Overview
===========

Ensemblify: A Python package for generating ensembles of intrinsically disordered regions of AlphaFold or user defined models
-----------------------------------------------------------------------------------------------------------------------------

.. image:: ../assets/ensemblify_presentation.svg
   :width: 100%
   :align: center

💡 What is Ensemblify?
----------------------

**Ensemblify** is a Python package that can generate protein conformational ensembles by sampling dihedral angle values from a three-residue fragment database and inserting them into flexible regions of a protein of interest (e.g. intrinsically disordered regions (IDRs)).

It supports both user-defined models and AlphaFold [1]_ predictions, using predicted Local Distance Difference Test (pLDDT) and Predicted Aligned Error (PAE) confidence metrics to guide conformational sampling. Designed to enhance the study of IDRs, it allows flexible customization of sampling parameters and works with single or multi-chain proteins, offering a powerful tool for protein structure research. Ensemble analysis and reweighting with experimental data is also available through interactive graphical dashboards.

🧰 How do I install Ensemblify?
-------------------------------

Step-by-step instructions for installing Ensemblify are available in the :ref:`Installation <Installation>` section.

After installing Ensemblify, make sure to visit the :ref:`Tripeptide Database <Tripeptide Database>` section to learn where you can get the database files required for ensemble generation.

💻 How can I use Ensemblify?
----------------------------

Ensemblify can be used either as a Command Line Interface (CLI) like so:

..  code-block:: bash
   
    conda activate ensemblify_env
    ensemblify [options]

or as a library inside a Python script or Jupyter notebook:

..  code-block:: python

    import ensemblify as ey
    ey.do_cool_stuff()

Check the :ref:`Usage <Usage>` section for more details.

You can also check out the interactive `Quick Reference Guide <examples/01_quick_reference_guide.ipynb>`_ notebook.

🔎 How does Ensemblify work?
----------------------------

A general overview of Ensemblify, descriptions of employed methods and applications can be found in the Ensemblify paper:

    PAPER 

📖 References
-------------

.. [1] J. Jumper, R. Evans, A. Pritzel et al., "Highly accurate protein structure prediction with AlphaFold," *Nature*, vol. 596, pp. 583-589, 2021. `Link <https://doi.org/10.1038/s41586-021-03819-2>`_