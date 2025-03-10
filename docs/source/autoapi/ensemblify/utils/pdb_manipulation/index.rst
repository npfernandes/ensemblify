ensemblify.utils.pdb_manipulation
=================================

.. py:module:: ensemblify.utils.pdb_manipulation

.. autoapi-nested-parse::

   Auxiliary functions to read, manipulate and write .pdb files.



Functions
---------

.. autoapisummary::

   ensemblify.utils.pdb_manipulation.df_from_pdb
   ensemblify.utils.pdb_manipulation.df_to_pdb
   ensemblify.utils.pdb_manipulation.extract_pdb_info
   ensemblify.utils.pdb_manipulation.cleanup_pdbs


Module Contents
---------------

.. py:function:: df_from_pdb(pdb: str) -> pandas.DataFrame

   Convert the information in a .pdb file into a pandas DataFrame using BioPDB.

   :param pdb: Path to the .pdb file.
   :type pdb: :py:class:`str`

   :returns:     The given pdb's information in DataFrame format.
   :rtype: pd.DataFrame


.. py:function:: df_to_pdb(df: pandas.DataFrame, output_pdb_filename: str)

   Write content of a DataFrame containing PDB file info as a .pdb file using BioPDB.

   :param df: DataFrame containing PDB information.
   :type df: :py:class:`pd.DataFrame`
   :param output_pdb_filename: Path to the output .pdb.
   :type output_pdb_filename: :py:class:`str`


.. py:function:: extract_pdb_info(pdb: str) -> dict[int, tuple[str, int, int]]

   Extract from a .pdb file info about number of chains, chain letters, starting residue
   numbers and chain size.

   :param topology: Path to .pdb topology file.
   :type topology: :py:class:`str`

   :returns:     Mapping of chain numbers to their letter, starting residue number and chain size.
   :rtype: dict[int,tuple[str,int,int]]


.. py:function:: cleanup_pdbs(pdbs: list[str])

   Delete all .pdb files in the given list.

   :param pdbs: Paths to .pdb files to delete.
   :type pdbs: :py:class:`list[str]`


