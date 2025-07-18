# The `analysis` module

With the `analysis` module, you can calculate structural properties averaged across your generated ensemble.

The resulting information is then presented in a user-friendly interactive graphical dashboard.

## Analyze your ensemble with an interactive graphical dashboard

To create an interactive graphical dashboard with structural information calculated from your conformational ensemble, provide Ensemblify with:

- your ensemble in trajectory format;
- your trajectory's corresponding topology file;
- the name you want to use for your protein in the created graphical dashboard.

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify analysis -trj trajectory.xtc -top topology.pdb -tid protein_name
   ```

   ```{code-tab} python Python
   from ensemblify.analysis import analyze_trajectory
   analyze_trajectory('trajectory.xtc','topology.pdb','protein_name')
   ```
````
