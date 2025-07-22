# 💻 Usage

Ensemblify offers four main modules, all of which can be accessed either through the command line or from inside a Python script/Jupyter Notebook.

```{toctree}
:titlesonly:
:maxdepth: 1
:hidden:
usage_generation.md
usage_conversion.md
usage_analysis.md
usage_reweighting.md
usage_clash_checking.md
```
- **[`generation`](usage_generation.md#the-generation-module) module:** generate a conformational ensemble for your protein of interest.
- **[`conversion`](usage_conversion.md#the-conversion-module) module:** convert your generated .pdb structures into a single .xtc trajectory file, facilitating ensemble storage and analysis.
- **[`analysis`](usage_analysis.md#the-analysis-module) module:** create an interactive graphical dashboard displaying structural information calculated from the conformational ensemble of your protein of interest.
- **[`reweighting`](usage_reweighting.md#the-reweighting-module) module:** use experimental SAXS data to reweight your conformational ensemble following the Bayesian/Maximum Entropy method.

----

Ensemblify also offers minor modules that can be accessed through the command line or from inside a Python script/Jupyter Notebook.

- **[`clash_checking`](usage_clash_checking.md#the-clash_checking-module) module:** check previously generated ensembles (even ones generated not using Ensemblify) for steric clashes, outputting detailed reports.
