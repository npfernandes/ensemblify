# CLI Workflow
In this short example, we will:

- Generate an ensemble of 10 conformations for Histatin 5, an intrinsically disordered protein with 24 residues;
- Convert the generated ensemble to trajectory format;
- Analyze the generated ensemble by calculating structural properties;
- Reweight the generated ensemble using experimental SAXS data, and re-calculate structural properties using the optimal set of conformer weights.

## 1. Generating an ensemble
First you must build you parameters file using the provided [template](docs/assets/parameters_template.yaml), optionally with the aid of the more user-friendly HTML [form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html).
Then, you can generate an ensemble using the Ensemblify CLI.

Assuming you have an Ensemblify parameters file `params.yaml` in the current working directory, you can run:

   ```{code-block} console
   $ ensemblify gen -p params.yaml
   Generating ensemble of 10 valid pdbs using 31 processor cores... 
   Ensemblified!: 100%|██████████| 10/10 [00:07<00:00,  1.43valid_pdb/s]   
   There are 10 valid pdbs, 14 were discarded ( 14 clashed | 0 violated constraints).
   Ensemble Generation Finished!
   ```

This command will create a folder in the current working directory named after the `job_name` in the parameters file.
Assuming `job_name` is 'Hst5', the ensemble will be located in `./Hst5/ensemble/valid_pdbs`.

## 2. Converting the generated ensemble to trajectory format
Assuming the current working directory is where the .pdb files making up the ensemble are stored:

   ```{code-block} console
   $ ensemblify con -e . -t ./trajectory -i Hst5
   Hst5 Trajectory creation complete! : 100%|██████████| 4/4 [00:00<00:00, 475.09step/s]
   ```

## 3. 