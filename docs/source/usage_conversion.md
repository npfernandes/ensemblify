# The `conversion` module
  
With the `conversion` module, you can convert your generated .pdb structures into a .xtc trajectory file.

This enables much easier storage and analysis of generated ensembles.

## Convert a conformational ensemble to trajectory format

To convert your generated .pdb structures into a single .xtc trajectory file, provide Ensemblify with:

<table class="tg">
<thead>
  <tr>
    <th>CLI Parameter</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-parameter"> Ensemble Directory (--ensembledir, -e) </td>
    <td class="tg-description"> Directory where the generated ensemble is stored. </td>
  </tr>
  <tr>
    <td class="tg-parameter"> Trajectory Directory (--trajectorydir, -t) </td>
    <td class="tg-description"> Directory where the trajectory file should be created. </td>
  </tr>
  <tr>
    <td class="tg-parameter"> Trajectory ID (--trajectoryid, -i) </td>
    <td class="tg-description"> Name for the trajectory file that will be created. </td>
  </tr>
</tbody>
</table>

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify conversion \
       -e ensemble_dir \
       -t trajectory_dir \
       -i trajectory_name
   ```

   ```{code-tab} python Python
   from ensemblify.conversion import ensemble2traj

   ensemble2traj(
       'ensemble_dir',
       'trajectory_dir',
       'trajectory_name'
   )
   ```
````
