# The `reweighting` module

With the `reweighting` module, you can use experimental SAXS data to reweigh your conformational ensemble following the Bayesian/Maximum Entropy method <sup>[[12]](#ref12)</sup>.

To do this, provide your ensemble in trajectory format, your trajectory's topology file, the name you want to use for your protein in the graphical dashboard and your experimental SAXS data.

Using the `ensemblify` command in a terminal:

    ensemblify reweighting -trj trajectory.xtc -top topology.pdb -tid trajectory_name -exp exp_SAXS_data.dat

Inside a Python script or Jupyter Notebook:

    from ensemblify.reweighting import reweight_ensemble
    reweight_ensemble('trajectory.xtc','topology.pdb','trajectory_name','exp_SAXS_data.dat')

Check the interactive [Reweighting Module](examples/05_reweighting_module.ipynb) notebook for detailed usage examples.

## References

<a id="ref12">[12]</a> S. Bottaro , T. Bengsten and K. Lindorff-Larsen, "Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach," pp. 219-240, Feb. 2020. In: Z. Gáspári, (eds) *Structural Bioinformatics*, *Methods in Molecular Biology*, vol. 2112, Humana, New York, NY. [[Link](https://doi.org/10.1007/978-1-0716-0270-6_15)]
