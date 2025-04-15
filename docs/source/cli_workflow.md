# CLI Workflow
In this short example, we will:

- [Generate](#1-generate-a-protein-conformational-ensemble) an ensemble of 10 conformations for Histatin 5, an intrinsically disordered protein (IDP) with 24 residues;
- [Convert](#2-convert-the-generated-ensemble-to-trajectory-format) the generated ensemble to trajectory format;
- [Analyze](#3-analyze-your-generated-ensemble) the generated ensemble by calculating structural properties;
- [Reweight](#4-reweight-your-generated-ensemble-using-experimental-saxs-data) the generated ensemble using experimental SAXS data, and re-calculate structural properties using the optimal set of conformer weights.

----

## 1. Generate a protein conformational ensemble

To generate a conformational ensemble for your protein of interest, you simply need to provide Ensemblify with your parameters file.

For this example, we will use a parameters file created to sample the conformational space of Histatin5, a small IDP that can normally be found in saliva, where it acts as a defense against fungal infections.

1. Build you parameters file using the provided [template](../assets/parameters_template.yaml), optionally with the aid of the more user-friendly HTML [form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html).

      In this example, we are using a parameters file `params.yaml` with the following information:

      ```{code-block} console
      {'job_name': 'Hst5',
       'sequence': 'DSHAKRHHGYKRKFHEKHHSHRGY',
       'size': 10,
       'databases': {'coil': '<path_to_database>'},
       'targets': {'A': [['MC', [1, 24], 'coil', 'TRIPEPTIDE']]},
       'output_path': '.'}
      ```

2. Generate an ensemble via the Ensemblify CLI:
    
      Assuming you have an Ensemblify parameters file `params.yaml` in the current working directory, you can run:

      ```{code-block} console
      (ensemblify_env) $ ensemblify gen -p params.yaml
      Generating ensemble of 10 valid pdbs using 31 processor cores... 
      Ensemblified!: 100%|██████████| 10/10 [00:07<00:00,  1.43valid_pdb/s]   
      There are 10 valid pdbs, 14 were discarded ( 14 clashed | 0 violated constraints).
      Ensemble Generation Finished!
      ```

      This command will create a directory in the current working directory named after the `job_name` parameter in `params.yaml`.
      Assuming `job_name` is 'Hst5', a directory named 'Hst5' will be created and the ensemble will be located in `./Hst5/ensemble/valid_pdbs`.

----

## 2. Convert the generated ensemble to trajectory format

After generating an ensemble, you can convert your set of PDB files into XTC format, a compressed **trajectory** format from GROMACS.
This allows for much more efficient storage of created ensembles with minimal loss of structural information.
Additionally, we can take advantage of the many available methods for analyzing files in this format.

To do so, you need only provide the location of your stored ensemble and the directory where you want to store your trajectory, optionally defining a prefix identifier for the created file.
Along with the created trajectory, one of the structures of the generated ensemble is saved as a **topology** PDB file that contains atomic connectivity information, as this is not stored in the trajectory and is required when using the `analysis` module.

1. Assuming the current working directory is the previously created 'Hst5' directory:

   ```{code-block} console
   (ensemblify_env) $ ensemblify con -e ./ensemble/valid_pdbs -t ./trajectory -i Hst5
   Hst5 Trajectory creation complete! : 100%|██████████| 4/4 [00:00<00:00, 475.09step/s]
   ```

   The created **trajectory** and **topology** files will be in `./trajectory/Hst5_trajectory.xtc` and `./trajectory/Hst5_top.xtc`, respectively.

----

## 3. Analyze your generated ensemble

After creating your trajectory, you can use it to create an interactive analysis dashboard with different plots and figures which will aid you in the structural analysis of your protein using your created ensemble.
To do this, specify the location of your **trajectory** and **topology** files and the output directory where you want to store your interactive dashboard and the figures and data used in its creation.

1. Assuming the current working directory is the previously created 'Hst5' directory:

   ```{code-block} console
   (ensemblify_env) $ ensemblify ana -trj ./trajectory/Hst5_trajectory.xtc -top ./trajectory/Hst5_top.xtc -tid Hst5 -out ./analysis
   Analyzing Hst5 trajectory...
   Calculating ramachandran data for Hst5...
   Calculating contact matrix for Hst5...
   Calculating contact matrix...: 100%|██████████| 10/10 [00:00<00:00, 810.87it/s]
   Calculating distance matrix for Hst5...
   Calculating distance matrix... : 100%|██████████| 10/10 [00:00<00:00, 2737.62it/s]
   Calculating secondary structure assignment frequency matrix for Hst5...
   Calculating structural metrics data for Hst5...
   Calculating rg...
   Calculating eed...
   Calculating dmax...
   Creating Hst5 analysis figures...
   Building ['Hst5'] analysis dashboard...
   Ensemble analysis calculation has finished. Please consult the interactive analysis_dashboard.html figure.
   ```

   You could now open your `analysis_dashboard.html` file in your browser and interpret your results.

----

## 4. Reweight your generated ensemble using experimental SAXS data

After generating an ensemble, you can use experimental SAXS data to reweight it.
This will create an interactive reweighting dashboard with comparisons between your uniformly weighted and reweighted ensembles, both in regards to fitting to experimental data and the calculation of structural properties of your protein.

To do this, specify the location of your **trajectory**, **topology** and experimental SAXS data files and the output directory where you want to store your interactive dashboard, along with the figures and data used in its creation.

1. Assuming the current working directory is the previously created 'Hst5' directory, and the `Hst5_SAXS.dat` file is one directory up:

   ```{code-block} console
   (ensemblify_env) $ ensemblify rew -trj ./trajectory/Hst5_trajectory.xtc -top ./trajectory/Hst5_top.xtc -tid Hst5 -exp ../Hst5_SAXS.dat -out ./reweighting
   Processing Hst5 experimental data file...
   Experimental errors on SAXS intensities have been corrected with BIFT using scale factor 1.0.
   Calculating Hst5 SAXS data... : 100%|██████████| 10/10 [00:00<00:00, 260.36it/s]
   Applying BME reweighting to Hst5 ensemble with theta values [1, 10, 20, 50, 75, 100, 200, 400, 750, 1000, 5000, 10000] ...
   Reweighting ensemble... : 100%|██████████| 12/12 [00:00<00:00, 93.66it/s]
   Please analyze the provided interactive figure (effective_frames_fit.html) and input the desired value(s) for the theta parameter.
   If more than one value, please separate them using a comma.
   Chosen theta value(s): 200.
   No contact matrix data was provided.
   Calculating contact matrix...: 100%|██████████| 10/10 [00:00<00:00, 910.24it/s]
   Calculating reweighted contact matrix...: 100%|██████████| 10/10 [00:00<00:00, 883.94it/s]
   No distance matrix data was provided.
   Calculating distance matrix... : 100%|██████████| 10/10 [00:00<00:00, 1757.51it/s]
   Calculating reweighted distance matrix... : 100%|██████████| 10/10 [00:00<00:00, 1744.50it/s]
   No secondary structure assignment frequency matrix data was provided.
   Calculating reweighted secondary structure assignment frequency matrix...
   No structural metrics distributions data was provided.
   Calculating rg...
   Calculating eed...
   Calculating dmax...
   Creating Hst5 reweighted interactive figures...
   Building Hst5 reweighting dashboard...
   Ensemble reweighting has finished. Please refer to the interactive reweighting_dashboard.html figure for analysis.
   ```
   
   You could now open your `reweighting_dashboard.html` file in your browser and interpret your results.