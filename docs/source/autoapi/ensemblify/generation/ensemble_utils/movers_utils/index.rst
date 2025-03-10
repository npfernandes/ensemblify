ensemblify.generation.ensemble_utils.movers_utils
=================================================

.. py:module:: ensemblify.generation.ensemble_utils.movers_utils

.. autoapi-nested-parse::

   Auxiliary functions for Custom PyRosetta Movers and to read database files into memory.



Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.movers_utils.read_database
   ensemblify.generation.ensemble_utils.movers_utils.trim_database
   ensemblify.generation.ensemble_utils.movers_utils.optimize_database
   ensemblify.generation.ensemble_utils.movers_utils.setup_databases
   ensemblify.generation.ensemble_utils.movers_utils.get_ss_bounds


Module Contents
---------------

.. py:function:: read_database(database_path: str) -> pandas.DataFrame

   Read a database file into a pandas DataFrame.

   If possible, read only the desired set of columns into a pandas.DataFrame
   (depends on database file format).

   :param database_path: Filepath to database file.
   :type database_path: :py:class:`str`

   :returns:     Database as a pandas.DataFrame.
   :rtype: pd.DataFrame


.. py:function:: trim_database(database: pandas.DataFrame, columns_to_keep: list[str])

   Removes columns in a database whose names are not in a given list.

   Modifies the database in place.

   :param database: Target database in DataFrame format.
   :type database: :py:class:`pd.DataFrame`
   :param columns_to_keep: Column names to keep in the database.
   :type columns_to_keep: :py:class:`list[str]`


.. py:function:: optimize_database(database: pandas.DataFrame) -> dict[str, pandas.DataFrame]

   Reduce a database's memory usage to what is strictly necessary.

   Datatypes of database's columns are optimized and database is broken
   into 20 pieces, one for each aminoacid residue.

   :param database: Unoptimized dihedral angle database.
   :type database: :py:class:`pd.DataFrame`

   :returns:     Mapping of aminoacid 1lettercode to their corresponding
                 dihedral angle values in the optimized database.
   :rtype: dict[str,pd.DataFrame]


.. py:function:: setup_databases(databases_paths: dict[str, str]) -> dict[str, dict[str, pandas.DataFrame]]

   Setup the databases the movers can access during sampling of dihedral angles.

   Databases are read into memory, trimmed into only necessary columns and optimized
   to use the least amount of memory possible.

   :param databases_paths: Mapping of database_ids to filepaths where the specified databases are stored.
   :type databases_paths: :py:class:`dict[str,str]`

   :returns:     Mapping of database_ids to mappings of aminoacid 1lettercode to their
                 corresponding dihedral angle values in the optimized database.
   :rtype: dict[str,dict[str,pd.DataFrame]]


.. py:function:: get_ss_bounds(secondary_structure: str) -> tuple[tuple[int, int], tuple[int, int]]

   Return the allowed range for the phi and psi angle values of a given secondary structure.

   :param secondary_structure: Identifier for a protein secondary structure.
   :type secondary_structure: :py:class:`str`

   :returns:

                 phi_bounds (tuple[int,int]):
                     Tuple with the lower and upper bounds for phi dihedral angle values
                     for the secondary structure in question.
                 psi_bounds (tuple[int,int]):
                     Tuple with the lower and upper bounds for psi dihedral angle values
                     for the secondary structure in question.
   :rtype: tuple[tuple[int,int],tuple[int,int]]


