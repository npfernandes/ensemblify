# The `pipeline` module

The `pipeline` module offers a convenient way to use all of Ensemblify's main modules in sequence, with a single command.

## Use Ensemblify's main features sequentially

The main features of the Ensemblify Python library include:

- Ensemble Generation using `ensemblify.generation.generate_ensemble`
- Trajectory Creation using `ensemblify.conversion.ensemble2traj`
- Trajectory Analysis using `ensemblify.analysis.analyze_trajectory`
- Ensemble Reweighting using `ensemblify.reweighting.reweight_ensemble`

To use all the main features of the Ensemblify Python library in sequence, provide:

- the path to your [parameters file](#-setting-up-your-parameters-file);
- a flag stating you want to analyze the generated ensemble;
- experimental SAXS data of your protein.

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify pipeline -p parameters.yaml -a -e exp_SAXS_data.dat
   ```

   ```{code-tab} python Python
   from ensemblify.pipeline import ensemblify_pipeline
   ensemblify_pipeline('parameters.yaml', True, 'exp_SAXS_data.dat')
   ```
````
