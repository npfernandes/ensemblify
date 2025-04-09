# 👾 Command Line Interface

You can access the Ensemblify Command Line Interface (CLI) using the `ensemblify` command.

    >>> ensemblify [options]

## Informational commands

### `ensemblify help`
Access information about the available Ensemblify modules.

   ```{code-block} console
   $ ensemblify help
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

### `ensemblify generation --help`
Access information about the Ensemblify `generation` module.

   ```{code-block} console
   $ ensemblify generation --help
   usage: ensemblify {generation, gen, g} [options]
   
   The generation module of the Ensemblify Python library.
   
   options:
       -h, --help        show this help message and exit
       -p, --parameters  Path to parameters file (.yaml).
   ```

### `ensemblify conversion --help`
Access information about the Ensemblify `conversion` module.

   ```{code-block} console
   $ ensemblify conversion --help
   usage: ensemblify {conversion, con, c} [options]
   
   The conversion module of the Ensemblify Python library.
   
   options:
       -h, --help           show this help message and exit
       -e, --ensembledir    Path to directory where ensemble files (.pdb) are located.
       -t, --trajectorydir  Path to directory where trajectory file (.xtc) will be created.
       -i, --trajectoryid   Name for created trajectory file (.xtc).
       -s, --size           (Optional) Number of frames of created trajectory file (.xtc).
   ```

### `ensemblify analysis --help`
Access information about the Ensemblify `analysis` module.

    >>> ensemblify analysis --help
    usage: ensemblify {analysis, ana, a} [options]

    The analysis module of the Ensemblify Python library.

    options:
        -h, --help                     show this help message and exit
        -trj, --trajectory             Path(s) to trajectory file(s) (.xtc).
        -top, --topology               Path(s) to topology file(s) (.pdb).
        -tid, --trajectoryid           Prefix identifier(s) for trajectory file(s).
        -out, --outputdir              (Optional) Path to output directory.
        -rma, --ramachandran           (Optional) Whether to calculate a dihedral angles matrix. Defaults to True.
        -dm, --distancematrix          (Optional) Whether to calculate a distance matrix. Defaults to True.
        -cm, --contactmatrix           (Optional) Whether to calculate a contact matrix. Defaults to True.
        -ssf, --ssfrequency            (Optional) Whether to calculate a secondary structure assignment frequency matrix. Defaults to True.
        -rg, --radiusgyration          (Optional) Whether to calculate a radius of gyration distribution. Defaults to True.
        -dmax, --maxdistance           (Optional) Whether to calculate a maximum distance distribution. Defaults to True.
        -eed, --endtoend               (Optional) Whether to calculate an end-to-end distance distribution. Defaults to True.
        -cmdist, --centremassdistance  (Optional) Pair(s) of MDAnalysis selection strings for which to calculate a centre of mass distance distribution. Defaults to None.
        -colors, --colorpalette        (Optional) List of color hexcodes to use, one for each analyzed trajectory.

### `ensemblify reweighting --help`
Access information about the Ensemblify `reweighting` module.

    >>> ensemblify reweighting --help
    usage: ensemblify {reweighting, rew, r} [options]

    The reweighting module of the Ensemblify Python library.

    options:
        -h, --help                  show this help message and exit
        -trj, --trajectory          Path to trajectory file (.xtc).
        -top, --topology            Path to topology file (.pdb).
        -tid, --trajectoryid        Prefix identifier for trajectory file.
        -exp, --expdata             Path to experimental SAXS data file (.dat).
        -out, --outputdir           (Optional) Path to output directory. Defaults to current working directory.
        -tht, --theta               (Optional) List of values to try as the theta parameter in BME.
        -cm, --contactmatrix        (Optional) Path to calculated contact matrix file (.csv).
        -dm, --distancematrix       (Optional) Path to calculated distance matrix file (.csv).
        -ss, --ssfrequency          (Optional) Path to calculated secondary structure frequency matrix file (.csv).
        -m, --metrics               (Optional) Path to calculated structural metrics file (.csv).
        -crg, --compare_rg          (Optional) Whether to calculate and compare uniform/reweighted radius of gyration distributions. Defaults to True.
        -cdmax, --compare_dmax      (Optional) Whether to calculate and compare uniform/reweighted maximum distance distributions. Defaults to True.
        -ceed, --compare_eed        (Optional) Whether to calculate and compare uniform/reweighted end-to-end distance distributions. Defaults to True.
        -ccmdist, --compare_cmdist  (Optional) Pair(s) of MDAnalysis selection strings for which to calculate and compare uniform/reweighted centre of mass distance distributions. Defaults to None.

### `ensemblify pipeline --help`
Access information about the Ensemblify `pipeline` module.

    >>> ensemblify pipeline --help
    usage: ensemblify {pipeline, ppl} [options]

    The pipeline module of the Ensemblify Python library.

    options:
        -h, --help        show this help message and exit
        -p, --parameters  Path to parameters file (.yaml).
        -a, --analysis    (Optional) Whether to perform the analysis of the ensemble. Defaults to False.
        -e, --expdata     (Optional) Path to experimental SAXS data file (.dat).

### `ensemblify clash_checking --help`
Access information about the Ensemblify `clash_checking` module.

    >>> ensemblify clash_checking --help
    usage: ensemblify {clash_checking, cch} [options]

    The clash checking module of the Ensemblify Python library.

    options:
        -h, --help             show this help message and exit
        -e, --ensembledir      Path to directory where ensemble .pdb structures are stored.
        -s, --samplingtargets  (Optional) Path to file (.yaml) with sampling targets: mapping of chain letters to residue ranges.
        -i, --inputstructure   (Optional) Path to input structure (.pdb) used to generate the ensemble.
