# The `generation` module
  
With the `generation` module, you can generate a conformational ensemble for your protein of interest.

## 📝 Setting up your parameters file

Before generating an ensemble, you must create a parameters file either:

- Using the provided [.html form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html);

- Directly, by editing the provided [parameters file template](../../examples/input_parameters/parameters_template.yml).

## Generate a conformational ensemble

To generate an ensemble, simply provide Ensemblify with the path to your parameters file.

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify generation -p parameters_file.yaml
   ```

   ```{code-tab} python Python
   from ensemblify.generation import generate_ensemble
   generate_ensemble('parameters_file.yaml')
   ```
````
