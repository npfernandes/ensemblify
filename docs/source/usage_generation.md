# The `generation` module
  
With the `generation` module, you can generate a conformational ensemble for your protein of interest.

## 📝 Setting up your parameters file

Before generating an ensemble, you must create a parameters file either:

- Using the provided [.html form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html);

- Directly, by editing the provided [parameters file template](../../examples/input_parameters/parameters_template.yml).

Below you can find a description of the minimum required parameters.

<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;}
.tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  overflow:hidden;padding:10px 5px;word-break:normal;}
.tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}
.tg .tg-43zj{border-color:inherit;color:#656565;text-align:center;vertical-align:middle}
.tg .tg-18eh{border-color:#000000;font-weight:bold;text-align:center;vertical-align:middle}
.tg .tg-d05w{border-color:#000000;color:#333333;text-align:center;vertical-align:middle}
.tg .tg-xnth{background-color:#2f4171;border-color:#000000;color:#ffffff;font-size:medium;font-weight:bold;text-align:center;
  vertical-align:middle}
.tg .tg-uzvj{border-color:inherit;font-weight:bold;text-align:center;vertical-align:middle}
.tg .tg-j844{border-color:inherit;color:#333333;text-align:center;vertical-align:middle}
</style>
<table class="tg"><thead>
  <tr>
    <th class="tg-xnth" colspan="2">Parameter</th>
    <th class="tg-xnth">Description</th>
  </tr></thead>
<tbody>
  <tr>
    <td class="tg-18eh" colspan="2">Job Name</td>
    <td class="tg-d05w">Name for generated files and folders.</td>
  </tr>
  <tr>
    <td class="tg-uzvj" rowspan="3">Input Structure/Sequence</td>
    <td class="tg-43zj">Full-length structure<br>available or Full IDP</td>
    <td class="tg-j844">Path to structure/sequence in .pdb/.txt format.</td>
  </tr>
  <tr>
    <td class="tg-43zj" rowspan="2">Folded domains<br>available, missing IDRs</td>
    <td class="tg-j844">Path to FASTA file(s) containing the sequences of all<br>(folded + disordered) protein domains, from N- to C-terminal.<br>Can be provided either in FASTA or Multi-FASTA format.</td>
  </tr>
  <tr>
    <td class="tg-j844">Path to PDB file(s) containing the structures of all folded<br>protein domains, from N- to C-terminal. Can be provided<br>as a single PDB file with multiple MODEL entries.</td>
  </tr>
  <tr>
    <td class="tg-uzvj" colspan="2">Size of Ensemble</td>
    <td class="tg-j844">Desired number of conformers in the generated ensemble.</td>
  </tr>
  <tr>
    <td class="tg-uzvj" colspan="2">Database(s)</td>
    <td class="tg-j844">Mapping of database IDs to the path of their respective<br>database files. Currently supported file formats<br>include .pkl, .csv and .parquet.</td>
  </tr>
  <tr>
    <td class="tg-uzvj" colspan="2">Sampling Target(s)</td>
    <td class="tg-j844">Protein regions to be targeted for conformational sampling.<br>You must assign for each desired sampling target a protein<br>chain, a range of residue numbers, a database ID to sample<br>from (matching one defined in Databases(s)) and a sampling<br>mode ('Tripeptide', if neighbouring residue information<br>is to be considered, or 'Single Residue' otherwise).</td>
  </tr>
  <tr>
    <td class="tg-uzvj" colspan="2">Output Path</td>
    <td class="tg-j844">Path to desired output directory. A directory named Job Name<br>will be created here, with all generated files and folders.</td>
  </tr>
</tbody></table>

## Generate a conformational ensemble

To generate an ensemble, simply provide Ensemblify with the path to your parameters file.

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify generation -p parameters_file.yml
   ```

   ```{code-tab} python Python
   from ensemblify.generation import generate_ensemble
   generate_ensemble('parameters_file.yml')
   ```
````
