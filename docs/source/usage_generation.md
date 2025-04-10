# The `generation` module
  
With the `generation` module, you can generate conformational ensembles for your protein of interest.

##  📝 Setting up your parameters file

Before generating an ensemble, you must create a parameters file either:

- Using the provided [parameters form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html);
- Directly, by editing the provided [parameters file template](../assets/parameters_template.yaml).

An [.html form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html) is provided to aid you in building your parameters file.

<details>  
  <summary><b>Parameters Form Preview</b></summary>

  ![alt text](../assets/parameters_form_preview.svg)

</details>

If you prefer to create your own parameters file from scratch, a [template file](docs/assets/parameters_template.yaml) is also provided.

## Generate a conformational ensemble

To generate an ensemble, provide Ensemblify with the path to your parameters file.

Using the `ensemblify` command in a terminal:

   ```{code-block} console
   $ ensemblify generation -p parameters_file.yaml
   ```

Inside a Python script or Jupyter Notebook:

   ```{code-block} python
   from ensemblify.generation import generate_ensemble
   generate_ensemble('parameters_file.yaml')
   ```
