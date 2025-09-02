# 🌐 Overview

<img src="../assets/ensemblify_presentation.svg" width="100%"/>

## 💡 What is Ensemblify?

**Ensemblify** is a Python package that can generate protein conformational ensembles by sampling dihedral angle values from a three-residue fragment database and inserting them into flexible regions of a protein of interest (e.g. intrinsically disordered regions (IDRs)).

It supports both user-defined models and AlphaFold<sup>[[1]](#ref1)</sup> predictions, using predicted Local Distance Difference Test (pLDDT) and Predicted Aligned Error (PAE) confidence metrics to guide conformational sampling. Designed to enhance the study of IDRs, it allows flexible customization of sampling parameters and works with single or multi-chain proteins, offering a powerful tool for protein structure research. Ensemble analysis and reweighting with experimental data is also available through interactive graphical dashboards.

----

## 🧰 How do I install Ensemblify?
Step-by-step instructions for installing Ensemblify are available in the [Installation](installation.md#-installation) section.

After installing Ensemblify, if you wish to generate conformational ensembles, make sure to visit the [Tripeptide Database](database.md#-tripeptide-database) section to learn where you can get the required database.

----

## 💻 How can I use Ensemblify?
Ensemblify can be used either as a Command Line Interface (CLI):

    conda activate ensemblify_env
    ensemblify [options]

or as a Python library inside a script or Jupyter notebook:

    import ensemblify as ey
    ey.show_config()

Check the [Usage](usage.md#-usage) section for more details.

----

## 🔎 How does Ensemblify work?
A general overview of Ensemblify, descriptions of employed methods and applications can be found in the Ensemblify [pre-print](https://www.biorxiv.org/content/10.1101/2025.08.26.672300v1) and accompanying support information.

----

## 📖 References

<a id="ref1">[1]</a> J. Jumper, R. Evans, A. Pritzel et al., "Highly accurate protein structure prediction with AlphaFold," *Nature*, vol. 596, pp. 583–589, 2021. [[Link](https://doi.org/10.1038/s41586-021-03819-2)]