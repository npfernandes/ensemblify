# The `reweighting` module

With the `reweighting` module, you can use experimental data to reweight your conformational ensemble following the Bayesian/Maximum Entropy (BME) method <sup>[[12]](#ref12)</sup>.

Fitting to experimental data and calculated ensemble structural properties are presented in a user-friendly interactive graphical dashboard.

Calculations are done for the ensemble before and after reweighting, facilitating comparisons.

## Reweight your conformational ensemble using experimental SAXS data

To use experimental SAXS data to reweight your conformational ensemble following the BME method, provide Ensemblify with:

- your ensemble in trajectory format;
- your trajectory's corresponding topology file;
- the name you want to use for your protein in the resulting graphical dashboard;
- the experimental SAXS data of your protein.

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify reweighting -trj trajectory.xtc -top topology.pdb -tid protein_name -exp exp_SAXS_data.dat
   ```

   ```{code-tab} python Python
   from ensemblify.reweighting import reweight_ensemble
   reweight_ensemble('trajectory.xtc','topology.pdb','trajectory_name','exp_SAXS_data.dat')
   ```
````

----

## References

<a id="ref12">[12]</a> S. Bottaro , T. Bengsten and K. Lindorff-Larsen, "Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach," pp. 219-240, Feb. 2020. In: Z. Gáspári, (eds) *Structural Bioinformatics*, *Methods in Molecular Biology*, vol. 2112, Humana, New York, NY. [[Link](https://doi.org/10.1007/978-1-0716-0270-6_15)]
