![alt text](docs/assets/ensemblify_logo.svg)

# Ensemblify: A Python package for generating ensembles of intrinsically disordered regions of AlphaFold or user defined models

![alt text](docs/assets/ensemblify_presentation.svg)

## What is Ensemblify?

Ensemblify is a python package that can generate protein conformational ensembles by sampling dihedral angle values from a three-residue fragment database and inserting them into flexible regions of a protein of interest (e.g intrinsically disordered regions (IDRs)).

It supports both user-defined models and AlphaFold predictions, using PAE and pLDDT confidence metrics to guide the conformational sampling process. Designed to enhance the study of IDRs, it allows flexible customization of sampling parameters and works with single or multi-chain proteins, offering a powerful tool for protein structure research. Ensemble analysis and reweighting with experimental data is also available through interactive graphical dashboards.

## How do I install Ensemblify?
Step-by-step instructions for installing Ensemblify are available in the [Installation](#installation) section.

## How can I use Ensemblify?
Ensemblify can be used either as a Command Line Interface (CLI) like so:

    conda activate ensemblify_env
    ensemblify [options]

or as a library inside a Python script or Jupyter notebook:

    import ensemblify as ey
    ey.do_cool_stuff()

Check the [Usage](#usage) section for more details and links to example Jupyter notebooks.

## How does Ensemblify work?
A general overview of Ensemblify, descriptions of employed methods and applications can be found in the Ensemblify paper:

    PAPER

<details><summary>

# Installation

</summary>    

## Ensemblify Python Package
It is heavily recommended to install the `ensemblify` Python package in a dedicated virtual environment.

You can create a new virtual environment using your favorite virtual environment manager. Examples shown will use `conda`. If you want to download `conda` you can do so through their [website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). We recommend [miniconda](https://docs.anaconda.com/miniconda/#quick-command-line-install), a free minimal installer for conda.

To install the `ensemblify` package, you can follow these commands:

1. Choose your current working directory where the `ensemblify` package will be installed. If you prefer, in your home directory you can create a new directory dedicated to the ensemblify installation and navigate into it by running:

    ```bash
    mkdir -p ~/ensemblify_installation
    cd ~/ensemblify_installation
    ```

2. Download and extract the `ensemblify` source code from this repository:

    ```bash
    wget https://github.com/npfernandes/ensemblify/archive/refs/heads/main.zip
    unzip main.zip
    ```

3. Create your `ensemblify_env` [Conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) with all of Ensemblify's python dependencies installed by using the provided [environment file](environment.yml) (recommended):

    ```bash
    cd ensemblify_main
    conda env create -f environment.yml
    ```

    or by creating the environment and installing the necessary python packages directly (not recommended):

    ```bash
    conda env create --name ensemblify_env
    conda activate ensemblify_env
    conda install --channel=conda-forge biopython contact_map MDAnalysis mdtraj numpy pandas plotly pyyaml scikit-sklearn scipy tqdm
    pip install ray["default"]
    ```

4. Install the `ensemblify` python package into your newly created environment.

    ```bash
    conda activate ensemblify_env
    pip install -e .
    ```

<!-- If you want to create a new `ensemblify_env` environment with Ensemblify and all its necessary python dependencies already installed you can use the provided conda environment file by running:

    wget https://raw.githubusercontent.com/npfernandes/ensemblify/main/environment.yml
    conda env create -f environment.yml -->

<!-- The `ensemblify_env` environment should be activated before installing Ensemblify in the usual manner:

    conda activate ensemblify_env    
    conda install ensemblify

Alternatively, Ensemblify is available via the Python Package Index:
    
    conda activate ensemblify_env   
    pip install ensemblify --upgrade -->

## Third Party Software
Each of Ensemblify's modules has different dependencies to third party software, so if you only plan on only using a certain module you do not have to install software required for others. The requirements are:

- `generation` module: [PyRosetta](#pyrosetta), [FASPR](#faspr) and [PULCHRA](#pulchra).

- `conversion` module: [GROMACS](#gromacs), [Pepsi-SAXS](#pepsi-saxs) and optionally [BIFT](#bift).

- `analysis` module: no third_party software required.

- `reweighting` module: no third_party software required.

### PyRosetta
PyRosetta is a Python-based interface to the powerful Rosetta molecular modeling suite [[1]](#ref1). Its functionalities are used and extended through Ensemblify in order to generate conformational ensembles. You can install it by following these commands:

1. Activate your `ensemblify_env` conda environment:

    ```bash
    conda activate ensemblify_env
    ```
    If you have not yet created it, check the [Ensemblify Python Package](#ensemblify-python-package) section.

2. Install the [`pyrosetta-installer`](https://pypi.org/project/pyrosetta-installer/) Python package, kindly provided by RosettaCommons, to aid in the `pyrosetta` installation:
    ```bash
    pip install pyrosetta-installer 
    ```

3. Use `pyrosetta-installer` to download (~ 1.6 GB) and install `pyrosetta` (note the distributed and serialization parameters):
    
    ```bash
    python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(distributed=True,serialization=True)'
    ```

4. To test your `pyrosetta` installation, you can type in a terminal:

    ```bash
    python -c 'import pyrosetta.distributed; pyrosetta.distributed.init()'
    ```

If this step does not produce a complaint or error, your installation has been successful.

Remember to re-activate the `ensemblify_env` conda environment each time you wish to run code that uses `pyrosetta`.

### FASPR

FASPR is an ultra-fast and accurate program for deterministic protein sidechain packing [[2]](#ref2). To compile the provided FASPR source-code, you can follow these commands:

For UNIX or Linux users:

1. Navigate to where the FASPR source code is located:

    ```bash
    cd ~/ensemblify_installation/ensemblify-main/src/ensemblify/third_party/FASPR-master/
    ```
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/FASPR-master/ -->

2. Compile the FASPR source code:

    ```bash
    g++ -O3 --fast-math -o FASPR src/*.cpp
    ```

3. If you are using a bash shell, you can register `faspr` as an alias for your FASPR executable by running:

    ```bash
    echo "alias faspr='$(realpath FASPR)'" >> ~/.bashrc
    source ~/.bashrc
    ```

For MacOS users:

1. Navigate to where the FASPR source code is located:

    ```bash
    cd ~/ensemblify_installation/ensemblify-main/src/ensemblify/third_party/FASPR-master/
    ```
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/FASPR-master/ -->

2. Compile the FASPR source code:

    ```bash
    g++ -03 -fast-math -o FASPR src/*.cpp
    ```

    or, if you get an error

    ```bash
    g++ -03 -o FASPR src/*.cpp
    ```

3. If you are using a bash shell, you can register `faspr` as an alias for your FASPR executable by running:

    ```bash
    echo "alias faspr='$(realpath FASPR)'" >> ~/.bashrc
    source ~/.bashrc
    ```

### PULCHRA
PULCHRA (PowerfUL CHain Restoration Algorithm) is a program for reconstructing full-atom protein models from reduced representations [[3]](#ref3). To compile the provided PULCHRA modified source-code, you can follow these commands:

1. Navigate to where the PULCHRA source code is located:

    ```bash
    cd ~/ensemblify_installation/ensemblify-main/src/ensemblify/third_party/pulchra-master/
    ```
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/pulchra-master/ -->

2. Compile the PULCHRA source code:

    ```bash
    cc -O3 -o pulchra pulchra_CHANGED.c pulchra_data.c -lm
    ```
    Some warnings will show up on your screen, this is normal and they can be ignored.

3. If you are using a bash shell, you can register `pulchra` as an alias for your PULCHRA executable by running:

    ```bash
    echo "alias pulchra='$(realpath pulchra)'" >> ~/.bashrc
    source ~/.bashrc
    ```

### GROMACS
GROMACS is a molecular dynamics package mainly designed for simulations of proteins, lipids, and nucleic acids [[4]](#ref4).
It comes with a large selection of flexible tools for trajectory analysis and the output formats are also supported by all major analysis and visualisation packages.

To download and compile the GROMACS source code from their [website](https://ftp.gromacs.org/gromacs/gromacs-2024.2.tar.gz) you can follow these commands:

1. Create and navigate into your desired GROMACS installation directory, for example:

    ```bash
    mkdir ~/GROMACS
    cd ~/GROMACS
    ```
<!-- mkdir $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/GROMACS/ -->
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/GROMACS/ -->

2. Download the GROMACS source code from their website:

    ```bash
    wget -O gromacs-2024.2.tar.gz https://zenodo.org/records/11148655/files/gromacs-2024.2.tar.gz?download=1
    ```

3. Follow the [GROMACS installation instructions](https://manual.gromacs.org/documentation/current/install-guide/index.html) to compile the GROMACS source code (this could take a while):

    ```bash
    tar xfz gromacs-2024.2.tar.gz
    cd gromacs-2024.2
    mkdir build
    cd build
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -j $(nproc)
    make
    make check
    sudo make install
    source /usr/local/gromacs/bin/GMXRC
    ```

The `gmx` command should already be registered as an alias for your GROMACS installation.

### PEPSI-SAXS
Pepsi-SAXS (Polynomial Expansions of Protein Structures and Interactions - SAXS) is an adaptive method for rapid and accurate computation of small-angle X-ray scattering (SAXS) profiles from atomistic protein models [[5]](#ref5).

To download the Pepsi-SAXS executable from their [website](https://team.inria.fr/nano-d/software/pepsi-saxs/) you can follow these commands:

For UNIX or Linux users:

1. Create and navigate into your desired Pepsi-SAXS installation directory, for example:

    ```bash
    mkdir ~/Pepsi-SAXS/
    cd ~/Pepsi-SAXS/
    ```
<!-- mkdir $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/Pepsi-SAXS/Linux_3.0/ -->
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/Pepsi-SAXS/Linux_3.0/ -->

2. Download and extract the Pepsi-SAXS Linux executable:

    ```bash
    wget -O Pepsi-SAXS-Linux.zip https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/Linux/3.0/Pepsi-SAXS-Linux.zip
    unzip Pepsi-SAXS-Linux.zip
    ```

3. If you are using a bash shell, you can register `pepsi_saxs` as an alias for your Pepsi-SAXS executable by running:

    ```bash
    echo "alias pepsi_saxs='$(realpath Pepsi-SAXS)'" >> ~/.bashrc
    source ~/.bashrc
    ```

For MacOS users:

1. Create and navigate into your desired Pepsi-SAXS installation directory, for example:

    ```bash
    mkdir ~/Pepsi-SAXS/
    cd ~/Pepsi-SAXS/
    ```
<!-- mkdir $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/Pepsi-SAXS/MacOS_2.6/ -->
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/Pepsi-SAXS/MacOS_2.6/ -->

2. Download and extract the Pepsi-SAXS MacOS executable:

    ```bash
    curl -O Pepsi-SAXS-MacOS.zip https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/MacOS/2.6/Pepsi-SAXS.zip
    unzip Pepsi-SAXS-MacOS.zip
    ```

3. If you are using a bash shell, you can register `pepsi_saxs` as an alias for your Pepsi-SAXS executable by running:

    ```bash
    echo "alias pepsi_saxs='$(realpath Pepsi-SAXS)'" >> ~/.bashrc
    source ~/.bashrc
    ```

### BIFT
Bayesian indirect Fourier transformation (BIFT) of small-angle experimental data allows for an estimation of parameters that describe the data [[6]](#ref6). Larsen *et al.* show in [[7]](#ref7) that BIFT can identify whether the experimental error in small-angle scattering data is over or underestimated. Here we use their implementation of this method to make this determination and scale the error values accordingly.

To compile the provided BIFT source code, you can follow these commands:

1. Navigate to where the BIFT source code is located:
    
    ```bash
    cd ~/ensemblify_installation/ensemblify-main/src/ensemblify/third_party/BIFT/
    ```
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third_party/BIFT/ -->

2. Compile the BIFT source code:

    ```bash
    gfortran -march=native -O3 bift.f -o bift
    ```
    the `-march=native` flag may be replaced with `-m64` or `-m32`, and it may be necessary to include the `-static` flag depending on which system you are on.

3. If you are using a bash shell, you can register `bift` as an alias for your bift executable by running:

    ```bash
    echo "alias bift='$(realpath bift)'" >> ~/.bashrc
    source ~/.bashrc
    ```

</details>

<details>  
  <summary>
  
  # Usage
  
  </summary>

Ensemblify offers four different modules, all of which can be acessed either through the command line or from inside a Python script or Jupyter Notebook.

## The `generation` module
With the `generation` module, you can generate conformational ensembles for your protein of interest.

Before generating an ensemble, you must create a parameters file either through the provided [parameters form](https://github.com/npfernandes/ensemblify/releases/download/v0.0.1-downloads/parameters_form.html) or directly by editing the provided [parameters file template](https://github.com/npfernandes/ensemblify/blob/main/docs/assets/parameters_template.yaml). Check the [parameters file setup](#setting-up-your-parameters-file) section for more details.

To generate an ensemble, provide Ensemblify with the path to your parameters file.

Using the `ensemblify` command in a terminal:

    ensemblify -g -p parameters_file.yaml

Inside a Python script or Jupyter Notebook:

    import ensemblify as ey
    ey.generate_ensemble('parameters_file.yaml')

Check the [example notebooks](examples/README.md) for detailed instructions.

### Setting up your parameters file
An [.html form](https://github.com/npfernandes/ensemblify/blob/main/docs/assets/parameters_form.html) is provided to aid you in building your parameters file.
<details>  
  <summary><b>Parameters Form Preview</b></summary>

  ![alt text](docs/assets/parameters_form_preview.svg)
</details>
<br>

If you prefer to create your own parameters file from scratch, a [template file](https://github.com/npfernandes/ensemblify/blob/main/docs/assets/parameters_template.yaml) is also provided.

## The `conversion` module
With the `conversion` module, you can convert your generated .pdb structures into a .xtc trajectory file, enabling you to easily store and manipulate your conformational ensemble.

To convert your ensemble of .pdb structures to a trajectory .xtc format, provide the name for your created trajectory, the directory where the ensemble is stored and the directory where the trajectory file should be created.

Using the `ensemblify` command in a terminal:

    ensemblify -c trajectory_name -i ensemble_dir -o trajectory_dir

Inside a Python script or Jupyter Notebook:

    import ensemblify as ey
    ey.ensemble2traj('trajectory_name','ensemble_dir','trajectory_dir')

Check the [example notebooks](examples/README.md) for detailed instructions.

## The `analysis` module
With the `analysis` module, you can create an interactive graphical dashboard displaying structural information calculated from the conformational ensemble of your protein of interest.

To create an interactive graphical dashboard, provide the path to your ensemble in trajectory format, the path to your trajectory's topology file, the name you want to use for your protein and what structural metrics you want to calculate (by default, all available metrics are calculated).

Using the `ensemblify` command in a terminal:

    ensemblify -a trajectory.xtc topology.pdb trajectory_name

Inside a Python script or Jupyter Notebook:

    import ensemblify as ey
    ey.analyze_trajectory('trajectory.xtc','topology.pdb','trajectory_name')

Check the [example notebooks](examples/README.md) for detailed instructions.

## The `reweighting` module
With the `reweighting` module, you can use experimental SAXS data to reweigh your ensemble following the BME method described in https://github.com/KULL-Centre/BME.

@incollection{bottaro2020integrating,
title={Integrating molecular simulation and experimental data: a Bayesian/maximum entropy reweighting approach},
author={Bottaro, Sandro and Bengtsen, Tone and Lindorff-Larsen, Kresten},
booktitle={Structural Bioinformatics},
pages={219--240},
year={2020},
publisher={Springer}
}

Using the `ensemblify` command in a terminal:

    ensemblify -r trajectory.xtc topology.pdb trajectory_name experimental_SAXS_data.dat

Inside a Python script or Jupyter Notebook:

    import ensemblify as ey
    ey.reweigh_ensemble('trajectory.xtc','topology.pdb','trajectory_name','experimental_SAXS_data.dat')

Check the [example notebooks](examples/README.md) for detailed instructions.

</details>

<details>
  <summary>
    
  # API Documentation
    
  </summary>

Ensemblify's API documentation is available at https://ensemblify.readthedocs.io/en/latest/index.html.

</details>

<details>  
  <summary>
  
  # Tripeptide Database
  
  </summary>   

Ensemblify provides a three-residue fragment (tripeptide) database from which to sample dihedral angles. This database was created and made available by González-Delgado *et al.* and, as described in [[8]](#ref8), it was built by extracting dihedral angles from structures taken from the SCOPe [[9]](#ref9) [[10]](#ref10) 2.07 release, a curated database of high-resolution experimentally determined protein structures.
In total, 6,740,433 tripeptide dihedral angle values were extracted, making up the *all* dataset. A structurally filtered dataset, *coil*, was generated by removing tripeptides contained in α-helices or β-strands, reducing the number of tripeptide dihedral angle values to 3,141,877.

## Using your own database
Ensemblify can sample dihedral angles from any file in a supported format (currently .parquet, .pkl or .csv), structured according to [Database Structure](#database-structure). Tripeptide sampling mode will only work if a tripeptide database is provided. However, single residue sampling mode will work even when you provide a tripeptide database.

### Database Structure

#### Tripeptide Database
Your database must contain at least 10 columns: 9 containing the Phi, Psi and Omega angles for each residue of the triplet (**in radians**) and 1 with the string identification of the fragment they make up. Any additional columns will be ignored.

| FRAG | OMG1 | PHI1 | PSI1 | OMG2 | PHI2 | PSI2 | OMG3 | PHI3 | PSI3 |
| :---: | :---: | :---: | :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |
| AAA | 3.136433 | -1.696219 | 1.100253 | -3.140388 | -2.765840 | 2.675006 | 3.140606 | -2.006085 | 2.063136 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| VYV | -3.135116 | -2.503945 | -0.949731 | -3.119968 | 1.407456 | 1.979130 | -3.112883 | -2.592680 | 2.573798 |

#### Single Residue Database
Your database must contain at least 4 columns: 3 containing the Phi, Psi and Omega angles for each residue (**in radians**)  and 1 with the string identification of the residue. Any additional columns will be ignored. Note the '2' suffix in the column names which help with compatibility between single residue and tripeptide sampling modes.

| FRAG | OMG2 | PHI2 | PSI2 |
| :---: | :---: | :---: | :---: |
| A | -3.140388 | -2.765840 | 2.675006 |
| ... | ... | ... | ... |
| Y | -3.119968 | 1.407456 | 1.979130 |

</details>

<details>  
  <summary>
  
  # Citation and Publications
  
  </summary>

If you use Ensemblify, please cite its original publication:

    PUB

</details>

<details>  
  <summary>
  
  # Acknowledgements
  
  </summary> 

We would like to thank the DeepMind team for developing AlphaFold.
We would also like to thank the team at the Juan Cortés lab in the LAAS-CNRS institute for creating the tripeptide database used in the development of this tool. Check out their work at https://moma.laas.fr/.

</details>

<details>  
  <summary>
  
  # References
  
  </summary>    

<a id="ref1">[1]</a> S. Chaudhury, S. Lyskov and J. J. Gray, "PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta," *Bioinformatics*, vol. 26, no. 5, pp. 689-691, Mar. 2010 [[Link](https://doi.org/10.1093/bioinformatics/btq007)]

<a id="ref2">[2]</a> X. Huang, R. Pearce and Y. Zhang, "FASPR: an open-source tool for fast and accurate protein side-chain packing," *Bioinformatics*, vol. 36, no. 12, pp. 3758-3765, Jun. 2020 [[Link](https://doi.org/10.1093/bioinformatics/btaa234)]

<a id="ref3">[3]</a> P. Rotkiewicz and J. Skolnick, "Fast procedure for reconstruction of full-atom protein models from reduced representations," *Journal of Computational Chemistry*, vol. 29, no. 9, pp. 1460-1465, Jul. 2008 [[Link](https://doi.org/10.1002/jcc.20906)] 

<a id="ref4">[4]</a> S. Pronk, S. Páll, R. Schulz, P. Larsson, P. Bjelkmar, R. Apostolov, M.R. Shirts, and J.C. Smith et al., “GROMACS 4.5: A high-throughput and highly parallel open source molecular simulation toolkit,” *Bioinformatics*, vol. 29, no. 7, pp. 845–854 2013 [[Link](https://doi.org/10.1093/bioinformatics/btt055)].

<a id="ref5">[5]</a> S. Grudinin, M. Garkavenko and A. Kazennov, "Pepsi-SAXS: an adaptive method for rapid and accurate computation of small-angle X-ray scattering profiles," *Structural Biology*, vol. 73, no. 5, pp. 449-464, May 2017 [[Link](https://doi.org/10.1107/S2059798317005745)]

<a id="ref6">[6]</a> B. Vestergaard and S. Hansen, "Application of Bayesian analysis to indirect Fourier transformation in small-angle scattering," *Journal of Applied Crystallography*, vol. 39, no. 6, pp. 797-804, Dec. 2006 [[Link](https://doi.org/10.1107/S0021889806035291)] 

<a id="ref7">[7]</a> A. H. Larsen and M. C. Pedersen, "Experimental noise in small-angle scattering can be assessed using the Bayesian indirect Fourier transformation," *Journal of Applied Crystallography*, vol. 54, no. 5, pp. 1281-1289, Oct. 2021 [[Link](https://doi.org/10.1107/S1600576721006877)]

<a id="ref8">[8]</a> J. González-Delgado , P. Bernadó , P. Neuvial and J. Cortés, "Statistical proofs of the interdependence between nearest neighbor effects on polypeptide backbone conformations," *Journal of Structural Biology*, vol. 214, no. 4, p. 107907, Dec. 2022 [[Link](https://doi.org/10.1016/j.jsb.2022.107907)]

<a id="ref9">[9]</a> N. K. Fox, S. E. Brenner and J. M. Chandonia, "SCOPe: Structural Classification of Proteins—extended, integrating SCOP and ASTRAL data and classification of new structures," *Nucleic Acids Research*, vol. 42, no. D1, pp. D304-D309, Jan. 2014 [[Link](https://doi.org/10.1093/nar/gkt1240)] 

<a id="ref10">[10]</a> J. M. Chandonia, N. K. Fox and S. E. Brenner, "SCOPe: classification of large macromolecular structures in the structural classification of proteins—extended database," *Nucleic Acids Research*, vol. 47, no. D1, pp. D475–D481, Jan. 2019 [[Link](https://doi.org/10.1093/nar/gky1134)]

<!-- [9] A. Estaña, N. Sibille, E. Delaforge, M. Vaisset, J. Cortés and P. Bernadó, "Realistic Ensemble Models of Intrinsically Disordered Proteins Using a Structure-Encoding Coil Database," *Structure* vol.27, no.2, pp. 381-391.e2, Feb. 2019 [[Link](https://doi.org/10.1016/j.str.2018.10.016)] -->

<!-- [10] J. M. Chandonia, L. Guan, S. Lin, C. Yu, N. K. Fox and S. E. Brenner, "SCOPe: Improvements to the Structural Classification of Proteins—extended Database to facilitate Variant Interpretation and Machine Learning," *Nucleic Acids Research*, vol. 50, no. D1, pp. D553–D559, Jan. 2022 [[Link](https://doi.org/10.1093/nar/gkab1054)] -->

</details>