# Ensemblify: A Python package for generating ensembles of intrinsically disordered regions of AlphaFold or user defined models

# Overview
Ensemblify works by sampling dihedral angle values from a three-residue fragment database and repeatedly inserting them into flexible regions of a protein of interest, generating conformational ensembles.

# Usage
To use the package you simply need to import it and call the desired function. If you want to generate an ensemble, provide the path to the parameters file that you created using the provided [parameters form](#setting-up-your-parameters-file).

Inside a Python script or Jupyter Notebook:

    import ensemblify as ey
    ey.generate_ensemble(PARAMETERS_FILEPATH)

In a terminal window:

    ensemblify -p PARAMETERS_FILEPATH

## Setting up your parameters file
An [.html form](docs/assets/parameters_form.html) is provided to aid you in building your parameters file.
<details open>  
  <summary><b>Parameters Form Preview</b></summary>

  ![alt text](docs/assets/parameters_form_preview_.svg)
</details>
<br>

If you prefer to create your own parameters file from scratch, a [template file](docs/assets/parameters_template.yaml) is also provided.

<details><summary>

# Installation

</summary>    

## Ensemblify Python Package
It is heavily recommended to install Ensemblify in a dedicated virtual environment. You can create a new virtual environment using your favorite virtual environment manager. Examples shown will use `conda`. If you want to download `conda` you can do so through their [website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

First create a folder to store downloaded files throughout the installation:

    mkdir -p ~/ensemblify_installation
    cd ~/ensemblify_installation

Then create your `ensemblify_env` conda environment with all of Ensemblify's python dependencies installed by running:

    wget https://raw.githubusercontent.com/npfernandes/ensemblify/main/environment.yml
    conda env create -f environment.yml

Finally download the ensemblify source code from this repository and install it into your newly created `ensemblify_env` conda environment.

    wget https://github.com/npfernandes/ensemblify/archive/refs/heads/main.zip
    unzip main.zip
    mv ensemblify-main ensemblify
    conda activate ensemblify_env
    pip install -e .

<!-- You can create a new environment named `ensemblify_env`, download the ensemblify source code from this repository and install it in your `ensemblify_env` environment (along with all its python dependencies) by running:

    conda env create ensemblify_env
    conda activate ensemblify_env
    wget https://github.com/npfernandes/ensemblify/archive/refs/heads/main.zip
    unzip main.zip
    cd ensemblify-main
    pip install -e .

If you want to create a new `ensemblify_env` environment with Ensemblify and all its necessary python dependencies already installed you can use the provided conda environment file by running:

    wget https://raw.githubusercontent.com/npfernandes/ensemblify/main/environment.yml
    conda env create -f environment.yml -->

<!-- The `ensemblify_env` environment should be activated before installing Ensemblify in the usual manner:

    conda activate ensemblify_env    
    conda install ensemblify

Alternatively, Ensemblify is available via the Python Package Index:
    
    conda activate ensemblify_env   
    pip install ensemblify --upgrade -->

## Third-Party Software
Each of Ensemblify's modules has different dependencies to third-party software, so if you only plan on only using a certain module you do not have to install software required for others. The dependencies are:

- `generation` module: [PyRosetta](#pyrosetta), [FASPR](#faspr) and [PULCHRA](#pulchra).

- `conversion` module: [GROMACS](#gromacs), [Pepsi-SAXS](#pepsi-saxs) and optionally [BIFT](#bift).

- `analysis` module: no third-party software required.

- `reweighting` module: no third-party software required.

### PyRosetta
PyRosetta is a Python-based interface to the powerful Rosetta molecular modeling suite [[1]](#ref1). Its functionalities are used and extended through Ensemblify in order to generate conformational ensembles.

Start by activating your newly created `ensemblify_env` conda environment:

    conda activate ensemblify_env

Then, install the [`pyrosetta-installer`](https://pypi.org/project/pyrosetta-installer/) Python package, kindly provided by RosettaCommons, to aid in the `pyrosetta` installation.

    pip install pyrosetta-installer 

Use `pyrosetta-installer` to download (~ 1.6 GB) and install `pyrosetta` (note the distributed and serialization parameters):

    python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(distributed=True,serialization=True)'

To test your `pyrosetta` installation, you can type in a terminal:

    python -c 'import pyrosetta.distributed; pyrosetta.distributed.init()'

If this step does not produce a complaint or error, your installation has been successful.

Remember to re-activate the `ensemblify_env` conda environment each time you wish to run code that uses `pyrosetta`.

### FASPR

FASPR is an ultra-fast and accurate program for deterministic protein sidechain packing [[2]](#ref2).

You will need to compile the provided FASPR source-code.

For UNIX or Linux users run:

    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/FASPR-master/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/FASPR-master/ -->
    g++ -O3 --fast-math -o FASPR src/*.cpp

For Mac users run:
    
    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/FASPR-master/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/FASPR-master/ -->
    g++ -03 -fast-math -o FASPR src/*.cpp

or

    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/FASPR-master/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/FASPR-master/ -->
    g++ -03 -o FASPR src/*.cpp

To register `faspr` as an alias for your FASPR executable in your bash shell run:

    PATH_TO_FASPR_EXECUTABLE=$(realpath FASPR)
    echo "alias faspr='$PATH_TO_FASPR_EXECUTABLE'" >> ~/.bashrc

### PULCHRA
PULCHRA (PowerfUL CHain Restoration Algorithm) is a program for reconstructing full-atom protein models from reduced representations [[3]](#ref3).

You will need to compile the provided PULCHRA modified source-code, by running:

    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/pulchra-master/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/pulchra-master/ -->
    cc -O3 -o pulchra pulchra_CHANGED.c pulchra_data.c -lm

To register `pulchra` as an alias for your PULCHRA executable in your bash shell run:

    PATH_TO_PULCHRA_EXECUTABLE=$(realpath pulchra)
    echo "alias pulchra='$PATH_TO_PULCHRA_EXECUTABLE'" >> ~/.bashrc 

### GROMACS
GROMACS is a molecular dynamics package mainly designed for simulations of proteins, lipids, and nucleic acids [[4]](#ref4).
It comes with a large selection of flexible tools for trajectory analysis and the output formats are also supported by all major analysis and visualisation packages.

You will need to download and compile the GROMACS source code from their [website](https://ftp.gromacs.org/gromacs/gromacs-2024.2.tar.gz).

For UNIX or Linux users run:

    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/GROMACS/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/GROMACS/ -->
    wget -O gromacs-2024.2.tar.gz https://zenodo.org/records/11148655/files/gromacs-2024.2.tar.gz?download=1

After downloading the GROMACS source-code, you can build it following the [GROMACS installation instructions](https://manual.gromacs.org/documentation/current/install-guide/index.html).

For UNIX or Linux users run:

    tar xfz gromacs-2024.2.tar.gz
    cd gromacs-2024.2
    mkdir build
    cd build
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
    make
    make check
    sudo make install
    source /usr/local/gromacs/bin/GMXRC

The `gmx` command should already be registered as an alias for your GROMACS installation.

### PEPSI-SAXS
Pepsi-SAXS (Polynomial Expansions of Protein Structures and Interactions - SAXS) is an adaptive method for rapid and accurate computation of small-angle X-ray scattering (SAXS) profiles from atomistic protein models [[5]](#ref5).

You will need to download and unzip the Pepsi-SAXS executable from their [website](https://team.inria.fr/nano-d/software/pepsi-saxs/).

For UNIX or Linux users run:

    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/Pepsi-SAXS/Linux_3.0/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/Pepsi-SAXS/Linux_3.0/ -->
    wget -O Pepsi-SAXS-Linux.zip https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/Linux/3.0/Pepsi-SAXS-Linux.zip
    unzip Pepsi-SAXS-Linux.zip

For MacOS users run:

    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/Pepsi-SAXS/MacOS_2.6/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/Pepsi-SAXS/MacOS_2.6/ -->
    curl -O Pepsi-SAXS-MacOS.zip https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/MacOS/2.6/Pepsi-SAXS.zip
    unzip Pepsi-SAXS-MacOS.zip

To register `pepsi-saxs` as an alias for your Pepsi-SAXS executable in your bash shell run:

    PATH_TO_Pepsi-SAXS_EXECUTABLE=$(realpath Pepsi-SAXS)
    echo "alias pepsi-saxs='PATH_TO_Pepsi-SAXS_EXECUTABLE'" >> ~/.bashrc

### BIFT
Bayesian indirect Fourier transformation (BIFT) of small-angle experimental data allows for an estimation of parameters that describe the data [[6]](#ref6). Larsen *et al.* show in [[7]](#ref7) that BIFT can identify whether the experimental error in small-angle scattering data is over or underestimated. Here we use their implementation of this method to make this determination and scale the error values accordingly.

You will need to compile the provided BIFT source-code by running:
    
    conda activate ensemblify_env
    cd ~/ensemblify_installation/ensemblify/third-party/BIFT/
<!-- cd $CONDA_PREFIX/lib/python3.10/ensemblify/third-party/BIFT/ -->
    gfortran -march=native -O3 bift.f -o bift

- the `-march=native` flag may be replaced with `-m64` or `-m32`, and it may be necessary to include the `-static` flag depending on which system you are on.

To register `bift` as an alias for your BIFT executable in your bash shell run:

    PATH_TO_BIFT_EXECUTABLE=$(realpath bift)
    echo "alias bift='PATH_TO_BIFT_EXECUTABLE'" >> ~/.bashrc
</details>

<details>  
  <summary>
  
  # Tripeptide Database
  
  </summary>   

Ensemblify provides a three-residue fragment (tripeptide) database from which to sample dihedral angles. As described by González-Delgado *et al.* in [[8]](#ref8), this database was built by extracting dihedral angles from structures taken from the SCOPe [[9]](#ref9) [[10]](#ref10) 2.07 release, a curated database of high-resolution experimentally determined protein structures.
In total, 6,740,433 tripeptide dihedral angle values were extracted, making up the *all* dataset. A structurally filtered dataset, *coil*, was generated by removing tripeptides contained in α-helices or β-strands, reducing the number of tripeptide dihedral angle values to 3,141,877.

## Using your own database
Ensemblify can sample dihedral angles from any file in a supported format (currently .pkl or .csv), structured according to [Database Structure](#database-structure). Tripeptide sampling mode will only work if a tripeptide database is provided. However, single residue sampling mode will work even when you provide a tripeptide database.

### Database Structure

#### Tripeptide Database
Your database must contain at least 10 columns: 9 containing the Phi, Psi and Omega angles for each residue of the triplet (**in radians**) and 1 with the string identification of the fragment they make up. Any additional columns will be ignored.

| fragment | Omega_res_1 | Phi_res_1 | Psi_res_1 | Omega_res_2 | Phi_res_2 | Psi_res_2 | Omega_res_3 | Phi_res_3 | Phi_res_3 |
| :---: | :---: | :---: | :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |
| AAA | 3.136433 | -1.696219 | 1.100253 | -3.140388 | -2.765840 | 2.675006 | 3.140606 | -2.006085 | 2.063136 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| VYV | -3.135116 | -2.503945 | -0.949731 | -3.119968 | 1.407456 | 1.979130 | -3.112883 | -2.592680 | 2.573798 |

#### Single Residue Database
Your database must contain at least 4 columns: 3 containing the Phi, Psi and Omega angles for each residue (**in radians**)  and 1 with the string identification of the residue. Any additional columns will be ignored. Note the '2' suffix in the column names which help with compatibility between single residue and tripeptide sampling modes.

| residue | Omega_res_2 | Phi_res_2 | Psi_res_2 |
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