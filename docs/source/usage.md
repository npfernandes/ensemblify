# 💻 Usage

Ensemblify offers four main modules, all of which can be accessed either through the command line or from inside a Python script/Jupyter Notebook.

## The `generation` module
  
With the `generation` module, you can generate conformational ensembles for your protein of interest.

Before generating an ensemble, you must create a parameters file either:

- Using the provided [parameters form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html);
- Directly, by editing the provided [parameters file template](docs/assets/parameters_template.yaml).
    
Check the [parameters file setup](#-setting-up-your-parameters-file) section for more details.

To generate an ensemble, provide Ensemblify with the path to your parameters file.

Using the `ensemblify` command in a terminal:

    ensemblify generation -p parameters_file.yaml

Inside a Python script or Jupyter Notebook:

    from ensemblify.generation import generate_ensemble
    generate_ensemble('parameters_file.yaml')

Check the interactive [Generation Module](examples/02_generation_module.ipynb) notebook for detailed usage examples. 

###  📝 Setting up your parameters file
  
An [.html form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html) is provided to aid you in building your parameters file.

<details>  
  <summary><b>Parameters Form Preview</b></summary>

  ![alt text](docs/assets/parameters_form_preview.svg)

</details>

If you prefer to create your own parameters file from scratch, a [template file](docs/assets/parameters_template.yaml) is also provided.

## The `conversion` module
  
With the `conversion` module, you can convert your generated .pdb structures into a .xtc trajectory file, enabling you to easily store and analyze your conformational ensemble.

To do this, provide the name for your created trajectory, the directory where the ensemble is stored and the directory where the trajectory file should be created.

Using the `ensemblify` command in a terminal:

    ensemblify conversion -j trajectory_name -e ensemble_dir -t trajectory_dir

Inside a Python script or Jupyter Notebook:

    from ensemblify.conversion import ensemble2traj
    ensemble2traj('trajectory_name','ensemble_dir','trajectory_dir')

Check the interactive [Conversion Module](examples/03_conversion_module.ipynb) notebook for detailed usage examples.

## The `analysis` module

With the `analysis` module, you can create an interactive graphical dashboard displaying structural information calculated from the conformational ensemble of your protein of interest.

To do this, provide your ensemble in trajectory format, your trajectory's topology file and the name you want to use for your protein in the graphical dashboard.

Using the `ensemblify` command in a terminal:

    ensemblify analysis -trj trajectory.xtc -top topology.pdb -tid trajectory_name

Inside a Python script or Jupyter Notebook:

    from ensemblify.analysis import analyze_trajectory
    analyze_trajectory('trajectory.xtc','topology.pdb','trajectory_name')

Check the interactive [Analysis Module](examples/04_analysis_module.ipynb) notebook for detailed usage examples.

## The `reweighting` module

With the `reweighting` module, you can use experimental SAXS data to reweigh your conformational ensemble following the Bayesian Maximum Entropy method [[12]](#ref12).

To do this, provide your ensemble in trajectory format, your trajectory's topology file, the name you want to use for your protein in the graphical dashboard and your experimental SAXS data.

Using the `ensemblify` command in a terminal:

    ensemblify reweighting -trj trajectory.xtc -top topology.pdb -tid trajectory_name -exp exp_SAXS_data.dat

Inside a Python script or Jupyter Notebook:

    from ensemblify.reweighting import reweight_ensemble
    reweight_ensemble('trajectory.xtc','topology.pdb','trajectory_name','exp_SAXS_data.dat')

Check the interactive [Reweighting Module](examples/05_reweighting_module.ipynb) notebook for detailed usage examples.
