# Informational commands
In this page you will find a collection of helpful commands that will provide information about each available Ensemblify module.

## `ensemblify help`
Access information about the available Ensemblify modules.

   ```{code-block} console
   (ensemblify_env) $ ensemblify help
   usage: ensemblify {generation, conversion, analysis, reweighting, pipeline, clash_checking, help} [module options]

   Command-line tool to access the modules of the Ensemblify Python library.

   positional arguments:

       help (h)                   Show this message and exit.
       generation (g, gen)        Access the generation module.
       conversion (c, con)        Access the conversion module.
       analysis (a, ana)          Access the analysis module.
       reweighting (r, rew)       Access the reweighting module.
       pipeline (ppl)             Access the pipeline module.
       clash_checking (cch)       Access the clash checking module.
   ```

## `ensemblify generation --help`
Access information about the Ensemblify `generation` module.

   ```{code-block} console
   (ensemblify_env) $ ensemblify generation --help
   usage: ensemblify {generation, gen, g} [options]

   The generation module of the Ensemblify Python library.
   
   options:
       -h, --help        show this help message and exit
       -p, --parameters  Path to parameters file (.yaml).
   ```

## `ensemblify conversion --help`
Access information about the Ensemblify `conversion` module.

   ```{code-block} console
   (ensemblify_env) $ ensemblify conversion --help
   usage: ensemblify {conversion, con, c} [options]

   The conversion module of the Ensemblify Python library.
   
   options:
       -h, --help           show this help message and exit
       -e, --ensembledir    Path to directory where ensemble files (.pdb) are located. Defaults to current working directory.
       -t, --trajectorydir  Path to directory where trajectory file (.xtc) will be created. Defaults to current working directory.
       -i, --trajectoryid   Prefix for created trajectory file (.xtc). Defaults to None.
       -s, --size           (Optional) Number of .pdb files to use for trajectory creation. Defaults to all .pdb files in the ensemble directory.
   ```

## `ensemblify analysis --help`
Access information about the Ensemblify `analysis` module.

   ```{code-block} console 
   (ensemblify_env) $ ensemblify analysis --help
   usage: ensemblify {analysis, ana, a} [options]

   The analysis module of the Ensemblify Python library.
   
   options:
       -h, --help              show this help message and exit
       -trj, --trajectory      Path(s) to trajectory file(s) (.xtc).
       -top, --topology        Path(s) to topology file(s) (.pdb).
       -tid, --trajectoryid    Prefix identifier(s) for trajectory file(s).
       -out, --outputdir       (Optional) Path to output directory. Defaults to current working directory.
       -ram, --ramachandran    (Optional) Whether to calculate a dihedral angles matrix. Defaults to True.
       -dmx, --distancematrix  (Optional) Whether to calculate a distance matrix. Defaults to True.
       -cmx, --contactmatrix   (Optional) Whether to calculate a contact matrix. Defaults to True.
       -ssf, --ssfrequency     (Optional) Whether to calculate a secondary structure assignment frequency matrix. Defaults to True.
       -rgy, --radiusgyration  (Optional) Whether to calculate a radius of gyration distribution. Defaults to True.
       -mxd, --maxdist         (Optional) Whether to calculate a maximum distance distribution. Defaults to True.
       -eed, --endtoend        (Optional) Whether to calculate an end-to-end distance distribution. Defaults to True.
       -cmd, --centermassdist  (Optional) Pair(s) of MDAnalysis selection strings for which to calculate a center of mass distance distribution. Defaults to None.
       -cls, --colors          (Optional) List of color hexcodes to use, one for each analyzed trajectory. Defaults to colorblind friendly color palette.
   ```

## `ensemblify reweighting --help`
Access information about the Ensemblify `reweighting` module.

   ```{code-block} console
   (ensemblify_env) $ ensemblify reweighting --help
   usage: ensemblify {reweighting, rew, r} [options]

   The reweighting module of the Ensemblify Python library.
   
   options:
       -h, --help              show this help message and exit
       -trj, --trajectory      Path to trajectory file (.xtc).
       -top, --topology        Path to topology file (.pdb).
       -tid, --trajectoryid    Prefix identifier for trajectory file.
       -exp, --expdata         Path to experimental SAXS data file (.dat).
       -out, --outputdir       (Optional) Path to output directory. Defaults to current working directory.
       -tht, --theta           (Optional) List of values to try as the theta parameter in BME. Defaults to [1, 10, 20, 50, 75, 100, 200, 400, 750, 1000, 5000, 10000]
       -cmx, --contactmatrix   (Optional) Path to calculated contact matrix file (.csv). Defaults to None.
       -dmx, --distancematrix  (Optional) Path to calculated distance matrix file (.csv). Defaults to None.
       -ssf, --ssfrequency     (Optional) Path to calculated secondary structure frequency matrix file (.csv). Defaults to None.
       -met, --metrics         (Optional) Path to calculated structural metrics file (.csv). Defaults to None.
       -rgy, --compare_rg      (Optional) Whether to calculate and compare uniform/reweighted radius of gyration distributions. Defaults to True.
       -mxd, --compare_dmax    (Optional) Whether to calculate and compare uniform/reweighted maximum distance distributions. Defaults to True.
       -eed, --compare_eed     (Optional) Whether to calculate and compare uniform/reweighted end-to-end distance distributions. Defaults to True.
       -cmd, --compare_cmdist  (Optional) Pair(s) of MDAnalysis selection strings for which to calculate and compare uniform/reweighted center of mass distance distributions. Defaults to None.
   ```

## `ensemblify pipeline --help`
Access information about the Ensemblify `pipeline` module.

   ```{code-block} console
   (ensemblify_env) $ ensemblify pipeline --help
   usage: ensemblify {pipeline, ppl} [options]

   The pipeline module of the Ensemblify Python library.
   
   options:
       -h, --help        show this help message and exit
       -p, --parameters  Path to parameters file (.yaml).
       -a, --analysis    (Optional) Whether to perform the analysis of the ensemble. Defaults to False.
       -e, --expdata     (Optional) Path to experimental SAXS data file (.dat). Defaults to None.
   ```

## `ensemblify clash_checking --help`
Access information about the Ensemblify `clash_checking` module.

   ```{code-block} console
   (ensemblify_env) $ ensemblify clash_checking --help
   usage: ensemblify {clash_checking, cch} [options]

   The clash checking module of the Ensemblify Python library.
   
   options:
       -h, --help             show this help message and exit
       -e, --ensembledir      Path to directory where ensemble .pdb structures are stored. Defaults to current working directory.
       -s, --samplingtargets  (Optional) Path to file (.yaml) with sampling targets: mapping of chain letters to residue ranges. Defaults to None.
       -i, --inputstructure   (Optional) Path to input structure (.pdb) used to generate the ensemble. Defaults to None.
   ```
