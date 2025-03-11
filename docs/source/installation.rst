.. _Installation:

🧰 Installation
===============

Ensemblify Python Package
-------------------------

It is **heavily** recommended to install the `ensemblify` Python package in a dedicated virtual environment.

You can create a new virtual environment using your favorite virtual environment manager. Examples shown will use `conda`. If you want to download `conda` you can do so through their `website <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. We recommend `miniconda <https://www.anaconda.com/docs/getting-started/miniconda/install>`_, a free minimal installer for conda.

To install the `ensemblify` Python package, you can follow these commands:

#. Get the `ensemblify` source code. To do this you:

   - Install Git if you haven't already:

     - On LINUX: ``sudo apt-get install git``
     - On macOS: Install Xcode Command Line Tools or use Homebrew: ``brew install git``

   - Clone this repository and ``cd`` into it:

   .. code-block:: bash

      git clone https://github.com/npfernandes/ensemblify.git
      cd ensemblify

#. Create your `ensemblify_env` `Conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ with all of Ensemblify's python dependencies installed by using the provided `environment file <environment.yml>`_ (recommended):

   .. code-block:: bash

      conda env create -f environment_LINUX.yml # or environment_MACOS.yml, for MacOS users
      conda activate ensemblify_env

   or by creating the environment and installing the necessary python packages manually (not recommended):

   .. code-block:: bash

      conda create --channel=conda-forge --name ensemblify_env python=3.10 MDAnalysis=2.6.1 mdtraj=1.9.9 numpy=1.26.4 pandas=2.2.2 pyarrow=13.0.0 scikit-learn=1.4.2 scipy=1.12.0 tqdm=4.66.2
      conda activate ensemblify_env
      pip install biopython==1.81 plotly==5.23.0 pyyaml==6.0.1 "ray[default]"==2.33.0

#. Install the `ensemblify` python package into your newly created environment.

   .. code-block:: bash

      pip install .

.. Alternatively, Ensemblify is available via the Python Package Index:

..    .. code-block:: bash

..       conda activate ensemblify_env   
..       pip install -U ensemblify

Third Party Software
--------------------

Each of Ensemblify's modules has different dependencies to third party software, so if you only plan on using a certain module you do not have to install software required for others. The requirements are:

- `generation` module: :ref:`PyRosetta <PyRosetta>`, :ref:`FASPR <FASPR>` and :ref:`PULCHRA <PULCHRA>`.
- `conversion` module: :ref:`GROMACS <GROMACS>`, :ref:`Pepsi-SAXS <Pepsi-SAXS>` and optionally :ref:`BIFT <BIFT>`.
- `analysis` module: no other software required.
- `reweighting` module: no other software required.

.. _PyRosetta:

PyRosetta
^^^^^^^^^

PyRosetta is a Python-based interface to the powerful Rosetta molecular modeling suite [2]_. Its functionalities are used through Ensemblify in order to generate conformational ensembles. You can install it by following these commands:

#. Activate your `ensemblify_env` conda environment:

   .. code-block:: bash

      conda activate ensemblify_env

   If you have not yet created it, check the `Ensemblify Python Package <#ensemblify-python-package>`_ section.

#. Install the `pyrosetta-installer <https://pypi.org/project/pyrosetta-installer/>`_ Python package, kindly provided by RosettaCommons, to aid in the `pyrosetta` installation:

   .. code-block:: bash

      pip install pyrosetta-installer 

#. Use `pyrosetta-installer` to download (~ 1.6 GB) and install `pyrosetta` (note the distributed and serialization parameters):

   .. code-block:: bash

      python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(distributed=True,serialization=True)'

#. To test your `pyrosetta` installation, you can type in a terminal:

   .. code-block:: bash

      python -c 'import pyrosetta.distributed; pyrosetta.distributed.init()'

   If this step does not produce a complaint or error, your installation has been successful.

   Remember to re-activate the `ensemblify_env` conda environment each time you wish to run code that uses `pyrosetta`.

.. _FASPR:

FASPR
^^^^^

FASPR is an ultra-fast and accurate program for deterministic protein sidechain packing [3]_. To compile the provided FASPR source-code, you can follow these commands:

For UNIX or Linux users:

#. Activate your `ensemblify_env` conda environment:

   .. code-block:: bash

      conda activate ensemblify_env

   If you have not yet created it, check the `Ensemblify Python Package <#ensemblify-python-package>`_ section.

#. Navigate to where the FASPR source code is located:

   .. code-block:: bash

      cd src/ensemblify/third_party/FASPR-master/ # assuming this cloned repository is your current working directory

#. Compile the FASPR source code:

   .. code-block:: bash

      g++ -O3 --fast-math -o FASPR src/*.cpp

#. You can add an environment variable with the path to your FASPR executable to your shell configuration file by running:

   .. code-block:: bash

      echo "export FASPR_PATH='$(realpath FASPR)'" >> ~/.bashrc # Or ~/.zshrc, depending on the shell
      source ~/.bashrc # Or ~/.zshrc, depending on the shell
      echo $FASPR_PATH # to check if the variable has been set correctly

   this will allow Ensemblify to know where your FASPR executable is located.

For MacOS users:

#. Activate your `ensemblify_env` conda environment:

   .. code-block:: bash

      conda activate ensemblify_env

   If you have not yet created it, check the `Ensemblify Python Package <#ensemblify-python-package>`_ section.

#. Navigate to where the FASPR source code is located:

   .. code-block:: bash

      cd src/ensemblify/third_party/FASPR-master/ # assuming this cloned repository is your current working directory

#. Compile the FASPR source code:

   .. code-block:: bash

      g++ -03 -fast-math -o FASPR src/*.cpp

   or, if you get an error

   .. code-block:: bash

      g++ -03 -o FASPR src/*.cpp

#. Add an environment variable with the path to your FASPR executable to your shell configuration file by running:

   .. code-block:: bash

      echo "export FASPR_PATH='$(realpath FASPR)'" >> ~/.bashrc # Or ~/.zshrc, depending on the shell
      source ~/.bashrc # Or ~/.zshrc, depending on the shell
      echo $FASPR_PATH # to check if the variable has been set correctly

   this will allow Ensemblify to know where your FASPR executable is located.

.. _PULCHRA:

PULCHRA
^^^^^^^

PULCHRA (PowerfUL CHain Restoration Algorithm) is a program for reconstructing full-atom protein models from reduced representations [4]_. To compile the provided PULCHRA modified source-code, you can follow these commands:

#. Activate your `ensemblify_env` conda environment:

   .. code-block:: bash

      conda activate ensemblify_env

   If you have not yet created it, check the `Ensemblify Python Package <#ensemblify-python-package>`_ section.

#. Navigate to where the PULCHRA source code is located:

   .. code-block:: bash

      cd src/ensemblify/third_party/pulchra-master/ # assuming this cloned repository is your current working directory

#. Compile the PULCHRA source code:

   .. code-block:: bash

      cc -O3 -o pulchra pulchra_CHANGED.c pulchra_data.c -lm

   Do not be alarmed if some warnings show up on your screen; this is normal and they can be ignored.

#. Add an environment variable with the path to your PULCHRA executable to your shell configuration file by running:

   .. code-block:: bash

      echo "export PULCHRA_PATH='$(realpath pulchra)'" >> ~/.bashrc # Or ~/.zshrc, depending on the shell
      source ~/.bashrc # Or ~/.zshrc, depending on the shell
      echo $PULCHRA_PATH # to check if the variable has been set correctly

   this will allow Ensemblify to know where your PULCHRA executable is located.

.. _GROMACS:

GROMACS
^^^^^^^

GROMACS is a molecular dynamics package mainly designed for simulations of proteins, lipids, and nucleic acids [5]_. It comes with a large selection of flexible tools for trajectory analysis and the output formats are also supported by all major analysis and visualisation packages.

To download and compile the GROMACS source code from their `website <https://ftp.gromacs.org/gromacs/gromacs-2024.2.tar.gz>`_ you can follow these commands:

#. Create and navigate into your desired GROMACS installation directory, for example:

   .. code-block:: bash

      mkdir -p ~/software/GROMACS
      cd ~/software/GROMACS

#. Download the GROMACS source code from their website:

   .. code-block:: bash

      wget -O gromacs-2024.2.tar.gz https://zenodo.org/records/11148655/files/gromacs-2024.2.tar.gz?download=1

#. Follow the `GROMACS installation instructions <https://manual.gromacs.org/documentation/current/install-guide/index.html>`_ to compile the GROMACS source code (this could take a while):

   .. code-block:: bash

      tar xfz gromacs-2024.2.tar.gz
      cd gromacs-2024.2
      mkdir build
      cd build
      cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
      make -j $(nproc)
      make check
      sudo make install
      source /usr/local/gromacs/bin/GMXRC

   Environment variables that will allow Ensemblify to know where GROMACS is located will have already been added to your shell configuration file.

.. _Pepsi-SAXS:

Pepsi-SAXS
^^^^^^^^^^

Pepsi-SAXS (Polynomial Expansions of Protein Structures and Interactions - SAXS) is an adaptive method for rapid and accurate computation of small-angle X-ray scattering (SAXS) profiles from atomistic protein models [6]_.

To download the Pepsi-SAXS executable from their `website <https://team.inria.fr/nano-d/software/pepsi-saxs/>`_ you can follow these commands:

For UNIX or Linux users:

#. Create and navigate into your desired Pepsi-SAXS installation directory, for example:

   .. code-block:: bash

      mkdir -p ~/software/Pepsi-SAXS/
      cd ~/software/Pepsi-SAXS/

#. Download and extract the Pepsi-SAXS Linux executable:

   .. code-block:: bash

      wget -O Pepsi-SAXS-Linux.zip https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/Linux/3.0/Pepsi-SAXS-Linux.zip
      unzip Pepsi-SAXS-Linux.zip

#. Add an environment variable with the path to your Pepsi-SAXS executable to your shell configuration file by running:

   .. code-block:: bash

      echo "export PEPSI_SAXS_PATH='$(realpath Pepsi-SAXS)'" >> ~/.bashrc # Or ~/.zshrc, depending on the shell
      source ~/.bashrc # Or ~/.zshrc, depending on the shell
      echo $PEPSI_SAXS_PATH # to check if the variable has been set correctly

   this will allow Ensemblify to know where your Pepsi-SAXS executable is located.

For MacOS users:

#. Create and navigate into your desired Pepsi-SAXS installation directory, for example:

   .. code-block:: bash

      mkdir -p ~/software/Pepsi-SAXS/
      cd ~/software/Pepsi-SAXS/

#. Download and extract the Pepsi-SAXS MacOS executable:

   .. code-block:: bash

      curl -O Pepsi-SAXS-MacOS.zip https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/MacOS/2.6/Pepsi-SAXS.zip
      unzip Pepsi-SAXS-MacOS.zip

#. Add an environment variable with the path to your Pepsi-SAXS executable to your shell configuration file by running:

   .. code-block:: bash

      echo "export PEPSI_SAXS_PATH='$(realpath Pepsi-SAXS)'" >> ~/.bashrc # Or ~/.zshrc, depending on the shell
      source ~/.bashrc # Or ~/.zshrc, depending on the shell
      echo $PEPSI_SAXS_PATH # to check if the variable has been set correctly

   this will allow Ensemblify to know where your Pepsi-SAXS executable is located.

.. _BIFT:

BIFT
^^^^

Bayesian indirect Fourier transformation (BIFT) of small-angle experimental data allows for an estimation of parameters that describe the data [7]_. Larsen *et al.* show in [8]_ that BIFT can identify whether the experimental error in small-angle scattering data is over or underestimated. Here we use their implementation of this method to make this determination and scale the error values accordingly.

To compile the provided BIFT source code, you can follow these commands:

#. Activate your `ensemblify_env` conda environment:

   .. code-block:: bash

      conda activate ensemblify_env

   If you have not yet created it, check the `Ensemblify Python Package <#ensemblify-python-package>`_ section.

#. Navigate to where the BIFT source code is located:

   .. code-block:: bash

      cd src/ensemblify/third_party/BIFT/ # assuming this cloned repository is your current working directory

#. Compile the BIFT source code:

   .. code-block:: bash

      gfortran -march=native -O3 bift.f -o bift

   the `-march=native` flag may be replaced with `-m64` or `-m32`, and it may be necessary to include the `-static` flag depending on which system you are on.

#. Add an environment variable with the path to your BIFT executable to your shell configuration file by running:

   .. code-block:: bash

      echo "export BIFT_PATH='$(realpath bift)'" >> ~/.bashrc # Or ~/.zshrc, depending on the shell
      source ~/.bashrc # Or ~/.zshrc, depending on the shell
      echo $BIFT_PATH # to check if the variable has been set correctly

   this will allow Ensemblify to know where your BIFT executable is located.

Do not forget to visit the :ref:`Tripeptide Database <Tripeptide Database>` section to learn where you can get the database files that are required for conformational ensemble generation. 

References
----------

.. [2] S. Chaudhury, S. Lyskov and J. J. Gray, "PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta," *Bioinformatics*, vol. 26, no. 5, pp. 689-691, Mar. 2010 `Link <https://doi.org/10.1093/bioinformatics/btq007>`_

.. [3] X. Huang, R. Pearce and Y. Zhang, "FASPR: an open-source tool for fast and accurate protein side-chain packing," *Bioinformatics*, vol. 36, no. 12, pp. 3758-3765, Jun. 2020 `Link <https://doi.org/10.1093/bioinformatics/btaa234>`_

.. [4] P. Rotkiewicz and J. Skolnick, "Fast procedure for reconstruction of full-atom protein models from reduced representations," *Journal of Computational Chemistry*, vol. 29, no. 9, pp. 1460-1465, Jul. 2008 `Link <https://doi.org/10.1002/jcc.20906>`_

.. [5] S. Pronk, S. Páll, R. Schulz, P. Larsson, P. Bjelkmar, R. Apostolov, M.R. Shirts, and J.C. Smith et al., “GROMACS 4.5: A high-throughput and highly parallel open source molecular simulation toolkit,” *Bioinformatics*, vol. 29, no. 7, pp. 845–854, 2013 `Link <https://doi.org/10.1093/bioinformatics/btt055>`_

.. [6] S. Grudinin, M. Garkavenko and A. Kazennov, "Pepsi-SAXS: an adaptive method for rapid and accurate computation of small-angle X-ray scattering profiles," *Structural Biology*, vol. 73, no. 5, pp. 449-464, May 2017 `Link <https://doi.org/10.1107/S2059798317005745>`_

.. [7] B. Vestergaard and S. Hansen, "Application of Bayesian analysis to indirect Fourier transformation in small-angle scattering," *Journal of Applied Crystallography*, vol. 39, no. 6, pp. 797-804, Dec. 2006 `Link <https://doi.org/10.1107/S0021889806035291>`_

.. [8] A. H. Larsen and M. C. Pedersen, "Experimental noise in small-angle scattering can be assessed using the Bayesian indirect Fourier transformation," *Journal of Applied Crystallography*, vol. 54, no. 5, pp. 1281-1289, Oct. 2021 `Link <https://doi.org/10.1107/S1600576721006877>`_
