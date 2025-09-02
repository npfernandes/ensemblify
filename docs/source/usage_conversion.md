# The `conversion` module
  
With the `conversion` module, you can convert your generated .pdb structures into a .xtc trajectory file.

This enables much easier storage and analysis of generated ensembles.

## Convert a conformational ensemble to trajectory format

To convert your generated .pdb structures into a single .xtc trajectory file, provide Ensemblify with:

- the directory where the generated ensemble is stored;
- the directory where the trajectory file should be created.
- the name for the trajectory file that will be created;

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify conversion -e ensemble_dir -t trajectory_dir -i trajectory_name
   ```

   ```{code-tab} python Python
   from ensemblify.conversion import ensemble2traj
   ensemble2traj('ensemble_dir','trajectory_dir','trajectory_name')
   ```
````
