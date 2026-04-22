# The `clash_checking` module

The `clash_checking` module offers a streamlined method to check if previously generated ensembles contain steric clashes using PULCHRA.

## Check a conformational ensemble for steric clashes

To check for steric clashes in your conformational ensemble, provide Ensemblify with:

<table class="tg">
<thead>
  <tr>
    <th>CLI Parameter</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-parameter"> Ensemble directory (--ensembledir, -e) </td>
    <td class="tg-description"> Directory where your .pdb files are stored. </td>
  </tr>
</tbody>
</table>

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify clash_checking -e ensemble_dir
   ```

   ```{code-tab} python Python
   from ensemblify.clash_checking import check_steric_clashes

   check_steric_clashes('ensemble_dir')
   ```
````
