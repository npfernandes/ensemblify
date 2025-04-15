# The `conversion` module
  
With the `conversion` module, you can convert your generated .pdb structures into a .xtc trajectory file, enabling you to easily store and analyze your conformational ensemble.

## Convert a conformational ensemble to trajectory format

To convert your generated .pdb structures into a .xtc trajectory file, provide Ensemblify with:

- the name for your created trajectory;
- the directory where the ensemble is stored;
- the directory where the trajectory file should be created.

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify conversion -e ensemble_dir -t trajectory_dir -i trajectory_name
   ```

   ```{code-tab} python Python
   from ensemblify.conversion import ensemble2traj
   ensemble2traj('ensemble_dir','trajectory_dir','trajectory_name')
   ```
````
