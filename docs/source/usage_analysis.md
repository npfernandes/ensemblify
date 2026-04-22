# The `analysis` module

With the `analysis` module, you can calculate structural properties averaged across your generated ensemble.

The resulting information is then presented in a user-friendly interactive graphical dashboard.

## Analyze your ensemble with an interactive graphical dashboard

To create an interactive graphical dashboard with structural information calculated from your conformational ensemble, provide Ensemblify with:

<table class="tg">
<thead>
  <tr>
    <th>CLI Parameter</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-parameter">Trajectory (--trajectory, -trj)</td>
    <td class="tg-description">Your generated ensemble in trajectory format.</td>
  </tr>
  <tr>
    <td class="tg-parameter">Topology (--topology, -top)</td>
    <td class="tg-description">Your trajectory's corresponding topology file.</td>
  </tr>
  <tr>
    <td class="tg-parameter">Trajectory ID (--trajectoryid, -tid)</td>
    <td class="tg-description">Name used to identify your protein in the created graphical dashboard.</td>
  </tr>
</tbody>
</table>

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify analysis \
       -trj trajectory.xtc \
       -top topology.pdb \
       -tid protein_name
   ```

   ```{code-tab} python Python
   from ensemblify.analysis import analyze_trajectory

   analyze_trajectory(
       'trajectory.xtc',
       'topology.pdb',
       'protein_name'
   )
   ```
````
