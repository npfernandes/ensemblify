# The `clash_checking` module

The `clash_checking` module offers a streamlined method to check if previously generated ensembles contain steric clashes using PULCHRA.

## Check if there are steric clashes present in a conformational ensemble

To check for steric clashes in your conformational ensemble, provide Ensemblify with the path to the directory where your .pdb files are stored.

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify clash_checking -e ensemble_dir
   ```

   ```{code-tab} python Python
   from ensemblify.clash_checking import check_steric_clashes
   check_steric_clashes('ensemble_dir')
   ```
````
