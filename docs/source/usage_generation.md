# The `generation` module
  
With the `generation` module, you can generate a conformational ensemble for your protein of interest.

## 📝 Setting up your parameters file

Before generating an ensemble, you must create a parameters file either:

- Using the provided [.html form](../../examples/parameters_form.html);

- Directly, by editing the provided [parameters file template](../../examples/parameters_template.yml).

Below you can find a description of the minimum required parameters.

<table class="tg"><thead>
  <tr>
    <th colspan="2">Parameter</th>
    <th>Description</th>
  </tr></thead>
<tbody>
  <tr>
    <td class="tg-parameter" colspan="2">Job Name</td>
    <td class="tg-description">Name for generated files and folders.</td>
  </tr>
  <tr>
    <td class="tg-parameter" rowspan="3">Input Structure/Sequence</td>
    <td class="tg-paramsec">Full-length structure<br>available or Full IDP</td>
    <td class="tg-description">Path to structure/sequence in .pdb/.txt format. A UniProt accession<br>code can also be provided, in which case the corresponding structure<br>from the AlphaFold Protein Structure Database is retrieved.</td>
  </tr>
  <tr>
    <td class="tg-paramsec" rowspan="2">Folded domains<br>available, missing IDRs<br>(Uses modelling module)</td>
    <td class="tg-description">Path to FASTA file(s) containing the sequences of all<br>(folded + disordered) protein domains, from N- to C-terminal.<br>Can be provided either in FASTA or Multi-FASTA format.</td>
  </tr>
  <tr>
    <td class="tg-description">Path to PDB file(s) containing the structures of all folded<br>protein domains, from N- to C-terminal. Can be provided<br>as a single PDB file with multiple MODEL entries.</td>
  </tr>
  <tr>
    <td class="tg-parameter" colspan="2">Size of Ensemble</td>
    <td class="tg-description">Desired number of conformers in the generated ensemble.</td>
  </tr>
  <tr>
    <td class="tg-parameter" colspan="2">Database(s)</td>
    <td class="tg-description">Mapping of database IDs to the path of their respective<br>database files. Currently supported file formats<br>include .pkl, .csv and .parquet.</td>
  </tr>
  <tr>
    <td class="tg-parameter" colspan="2">Sampling Target(s)</td>
    <td class="tg-description">Protein regions to be targeted for conformational sampling.<br>You must assign for each desired sampling target a protein<br>chain, a range of residue numbers, a database ID to sample<br>from (matching one defined in Databases(s)) and a sampling<br>mode ('Tripeptide', if neighbouring residue information<br>is to be considered, or 'Single Residue' otherwise).</td>
  </tr>
  <tr>
    <td class="tg-parameter" colspan="2">Output Path</td>
    <td class="tg-description">Path to desired output directory. A directory named Job Name<br>will be created here, with all generated files and folders.</td>
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
