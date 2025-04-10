# The `conversion` module
  
With the `conversion` module, you can convert your generated .pdb structures into a .xtc trajectory file, enabling you to easily store and analyze your conformational ensemble.

To do this, provide:
- the name for your created trajectory;
- the directory where the ensemble is stored;
- the directory where the trajectory file should be created.

Using the `ensemblify` command in a terminal:

   ```{code-block} console
   $ ensemblify conversion -j trajectory_name -e ensemble_dir -t trajectory_dir
   ```

Inside a Python script or Jupyter Notebook:

   ```{code-block} python
   from ensemblify.conversion import ensemble2traj
   ensemble2traj('trajectory_name','ensemble_dir','trajectory_dir')
   ```