# Ensemblify: A Python package for generating ensembles of intrinsically disordered regions of AlphaFold or user defined models

<img src="docs/assets/ensemblify_presentation.svg" width="100%"/>

<div align="justify">

### 💡 What is Ensemblify?

**Ensemblify** is a Python package that can generate protein conformational ensembles by sampling dihedral angle values from a three-residue fragment database and inserting them into flexible regions of a protein of interest (*e.g.* intrinsically disordered regions (IDRs)).

It supports both user-defined models and AlphaFold<sup>[[1]](https://doi.org/10.1038/s41586-021-03819-2)</sup> predictions, using the predicted Local Distance Difference Test (pLDDT) and Predicted Aligned Error (PAE) confidence metrics to guide conformational sampling. Designed to enhance the study of IDRs, it allows flexible customization of sampling parameters and works with single or multi-chain proteins, offering a powerful tool for protein structure research. Ensemble analysis and reweighting with experimental data is also available through interactive graphical dashboards.

### 🧰 How do I install Ensemblify?
Step-by-step instructions for installing Ensemblify are available in the [Documentation](#https://ensemblify.readthedocs.io/latest/installation.html).

### 💻 How can I use Ensemblify?
Ensemblify can be used either as a Command Line Interface (CLI):

    conda activate ensemblify_env
    ensemblify [options]

or as a Python library inside a script or Jupyter notebook:

    import ensemblify as ey
    ey.do_cool_stuff()

Check the [Documentation](https://ensemblify.readthedocs.io/latest/usage.html) for more details.

You can also check out the [Quick Reference Guide](examples/quick_reference_guide.ipynb) notebook for a basic rundown of Ensemblify's features. 

### 🔎 How does Ensemblify work?
A general overview of Ensemblify, descriptions of employed methods and applications can be found in the Ensemblify paper:

    PAPER

## 🗃 Tripeptide Database

Ensemblify provides a three-residue fragment (tripeptide) database from which to sample dihedral angle values.

This database is provided separately from the Ensemblify source-code.

You can learn where to get it and more about its creation in the [Documentation](https://ensemblify.readthedocs.io/latest/database.html).

## 📚 Accessing Documentation

Ensemblify's documentation is available together with an API reference at https://ensemblify.readthedocs.io.
Alternatively, the source-code contains docstrings with relevant information. 

## 🗨️ Citation and Publications

If you use Ensemblify, please cite its original publication:

    PUB

## 🤝 Acknowledgements

We would like to thank the DeepMind team for developing AlphaFold.

We would also like to thank the team at the Juan Cortés lab in the LAAS-CNRS institute for creating the tripeptide database used in the development of this tool. Check out their work at https://moma.laas.fr/.

## ✍️ Authors

**Nuno P. Fernandes** (Main Developer) [[GitHub]](https://github.com/npfernandes?tab=repositories)

**Tiago Lopes Gomes** (Initial prototyping, Supervisor) [[GitHub]](https://github.com/TiagoLopesGomes?tab=repositories)

**Tiago N. Cordeiro** (Supervisor) [[GitHub]](https://github.com/CordeiroLab?tab=repositories)

</div>
