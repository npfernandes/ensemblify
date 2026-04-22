# The `reweighting` module

With the `reweighting` module, you can use experimental data to reweight your conformational ensemble following the Bayesian/Maximum Entropy (BME) method <sup>[[12]](#ref12)</sup>.

Fitting to experimental data and calculated ensemble structural properties are presented in a user-friendly interactive graphical dashboard.

Calculations are done for the ensemble before and after reweighting, facilitating comparisons.

## Reweight your conformational ensemble using experimental SAXS data

To use experimental SAXS data to reweight your conformational ensemble following the BME method, provide Ensemblify with:

<table class="tg">
<thead>
  <tr>
    <th>CLI Parameter</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-parameter"> Trajectory (--trajectory, -trj) </td>
    <td class="tg-description"> Your generated ensemble in trajectory format. </td>
  </tr>
  <tr>
    <td class="tg-parameter"> Topology (--topology, -top) </td>
    <td class="tg-description"> Your trajectory's corresponding topology file. </td>
  </tr>
  <tr>
    <td class="tg-parameter"> Trajectory ID (--trajectoryid, -tid) </td>
    <td class="tg-description"> Name used to identify your protein in the created graphical dashboard. </td>
  </tr>
  <tr>
    <td class="tg-parameter"> Experimental SAXS data (--expdata, -exp) </td>
    <td class="tg-description"> Experimental SAXS data of your protein. </td>
  </tr>
</tbody>
</table>

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify reweighting \
       -trj trajectory.xtc \
       -top topology.pdb \
       -tid protein_name \
       -exp exp_SAXS_data.dat
   ```

   ```{code-tab} python Python
   from ensemblify.reweighting import reweight_ensemble

   reweight_ensemble(
       'trajectory.xtc',
       'topology.pdb',
       'trajectory_name',
       'exp_SAXS_data.dat'
   )
   ```
````

----

## References

<a id="ref12">[12]</a> S. Bottaro , T. Bengsten and K. Lindorff-Larsen, "Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach," pp. 219-240, Feb. 2020. In: Z. Gáspári, (eds) *Structural Bioinformatics*, *Methods in Molecular Biology*, vol. 2112, Humana, New York, NY. [[DOI](https://doi.org/10.1007/978-1-0716-0270-6_15)]
