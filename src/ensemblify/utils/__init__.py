"""
Utils - `ensemblify.utils`
=============================================

:Author(s): Nuno P. Fernandes
:Year: 2024
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains auxilliary functions and classes used in other modules that can be useful
in other applications.

Example applications
--------------------

- Read a .pdb file into a `pandas.DataFrame`
  ------------------------------------------
  The `ensemblify.utils.df_from_pdb` function can be used to read a .pdb file and get a
  `pandas.DataFrame` with its contents. For example, we can get the contents of a .pdb file
  of Histatin5, an intrinsically disordered protein (IDP) with 24 aminoacid residues.
  Assuming the path to the required input .pdb file is assigned to a variable named HST5_PDB,
  you should run:

```
>>> import ensemblify as ey
>>> hst5_df = ey.df_from_pdb(HST5_PDB)
```

- Write a .pdb file from a `pandas.DataFrame`
  -------------------------------------------
  The `ensemblify.utils.df_to_pdb` function can be used to write a .pdb file from the content
  of a `pandas.DataFrame`. For example, we can write the contents of the Hst5 DataFrame we created
  before back into a .pdb file. Assuming the path to the .pdb file to be created is assigned to a
  variable named NEW_HST5_PDB, you should run:

```
>>> import ensemblify as ey
>>> ey.df_to_pdb(hst5_df, NEW_HST5_PDB)
```

- Extract chain information from a .pdb file
  ------------------------------------------
  The `ensemblify.utils.extract_pdb_info` function can be used to extract from a .pdb file
  information regarding the number of protein chains present, which chain letters identify them,
  their starting residue numbers and their size. For example, we can extract information from a
  .pdb file of Histatin5, an intrinsically disordered protein (IDP) with 24 aminoacid residues.
  Assuming the path to the required .pdb file is assigned to a variable named HST5_PDB, you should
  run:

```
>>> import ensemblify as ey
>>> ey.extract_pdb_info(HST_PDB)
```

Available Functions
----------------
- `df_from_pdb`

      Convert the information in a .pdb file into a pandas DataFrame using BioPDB.

- `df_to_pdb`

      Write content of a DataFrame containing PDB file info as a .pdb file using BioPDB.

- `extract_pdb_info`

      Extract from a .pdb file info about number of chains, chain letters, starting residue
      numbers and chain size.

- `cleanup_pdbs`

      Delete all .pdb files in the given list.

- `kde`

      Calculate a Kernel Density Estimate (KDE) distribution for a given dataset.

"""

from ensemblify.utils.misc import kde
from ensemblify.utils.pdb_manipulation import (
    cleanup_pdbs,
    df_from_pdb,
    df_to_pdb,
    extract_pdb_info,
)

__all__ = ['df_from_pdb',
           'df_to_pdb',
           'extract_pdb_info',
           'cleanup_pdbs',
           'kde']
