ensemblify.generation.ensemble_utils.pdb_processing
===================================================

.. py:module:: ensemblify.generation.ensemble_utils.pdb_processing

.. autoapi-nested-parse::

   Auxiliary functions for processing .pdb files from sampling output.



Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.pdb_processing.check_clashes
   ensemblify.generation.ensemble_utils.pdb_processing.setup_logger
   ensemblify.generation.ensemble_utils.pdb_processing.apply_faspr_single
   ensemblify.generation.ensemble_utils.pdb_processing.apply_rewrite_single
   ensemblify.generation.ensemble_utils.pdb_processing.apply_pulchra_single
   ensemblify.generation.ensemble_utils.pdb_processing.apply_restore_single
   ensemblify.generation.ensemble_utils.pdb_processing.process_pdb


Module Contents
---------------

.. py:function:: check_clashes(sampled_pdb: str, pulchra_output_buffer: str, sampling_targets: dict[str, tuple[tuple[str, tuple[int, Ellipsis], str, str]]], input_clashes: list[tuple[str, str]] | None) -> bool

   Check if there are recorded steric clashes in given PULCHRA output.

   Clashes present in input structure are not considered.
   Clashes are only considered when at least one residue belongs to a sampled region.

   :param sampled_pdb: Filepath to .pdb file output from conformational sampling.
   :type sampled_pdb: :py:class:`str`
   :param pulchra_output_buffer: Stdout from applying PULCHRA to the sampled .pdb structure.
   :type pulchra_output_buffer: :py:class:`str`
   :param sampling_targets: Mapping of chain identifiers to sampled residue numbers.
   :type sampling_targets: :py:class:`dict[str,tuple[tuple[str,tuple[int,...],str,str]]]`
   :param input_clashes: Clashes present in the sampling input structure, that will be ignored if
                         present in the given PULCHRA output.
   :type input_clashes: :py:class:`list[tuple[str,str]] | None`

   :returns:     True if the given PULCHRA output mentions clashes not already present in sampling
                 input structure, False otherwise.
   :rtype: bool


.. py:function:: setup_logger(pdb: str, log_file: str) -> logging.Logger

   Setup a Logger object for this pdb, with output to log_file.

   :param pdb: Path to .pdb file, will be the name of the logger.
   :type pdb: :py:class:`str`
   :param log_file: Filepath to log file.
   :type log_file: :py:class:`str`

   :returns:     Logger object.
   :rtype: logger (logging.Logger)


.. py:function:: apply_faspr_single(faspr_path: str, pdb: str) -> str | None

   Apply FASPR to a .pdb file. Log outcome.

   :param faspr_path: Path to FASPR executable or its alias.
   :type faspr_path: :py:class:`str`
   :param pdb: Path to .pdb file.
   :type pdb: :py:class:`str`

   :returns:     Path to .pdb file from FASPR output, with filename equal to
                 input .pdb with the suffix '_faspr' added.
   :rtype: str | None

   :raises subprocess.CalledProcessError if FASPR was not applied succesfully.:


.. py:function:: apply_rewrite_single(pdb: str) -> str

   Convert a .pdb file into single chain with sequential numbering.

   Necessary when a multichain .pdb will be input into Pulchra, as it does
   not support multiple chain structures.
   The output modified version has the _rewrite suffix added to its name.

   :param pdb: Path to input .pdb file for conversion.
   :type pdb: :py:class:`str`

   Returns
       str:
           Path to modified .pdb. Filename is the same as the input, with
           _rewrite suffix added.


.. py:function:: apply_pulchra_single(pulchra_path: str, pdb: str) -> tuple[str, str] | tuple[None, None]

   Apply PULCHRA to a .pdb file. Log outcome.

   :param pulchra_path: Path to PULCHRA executable or its alias.
   :type pulchra_path: :py:class:`str`
   :param pdb: Path to .pdb file.
   :type pdb: :py:class:`str`

   :returns:

                 rebuilt_filename (str | None):
                     Path to PULCHRA output structure. Same filename as input .pdb
                     with added .rebuilt suffix.
                 pulchra_output (str | None):
                     PULCHRA stdout used for later clash checking.
   :rtype: tuple[str,str] | tuple[None,None]

   :raises subprocess.CalledProcessError if PULCHRA was not applied sucessfully.:


.. py:function:: apply_restore_single(pdb: str, reference_pdb: str) -> str

   Restore chain, residue number and B-Factor info to pdb from reference pdb.

   Restore chain, residue numbering and B-Factor information to a post-Pulchra .pdb
   file, following the information in a reference .pdb file (either the first output
   of sampling process or the sampling input .pdb).

   :param pdb: Path to the PULCHRA output .pdb structure (ending in .rebuilt suffix).
   :type pdb: :py:class:`str`
   :param reference_pdb: Path to .pdb file to use as reference for restoring the chains and
                         residue numbering.
   :type reference_pdb: :py:class:`str`

   :returns:     Path to the .pdb structure with restored chain and residue numbering.
                 Filename matches that of input, with _restores suffix added.
   :rtype: str


.. py:function:: process_pdb(sampled_pdb: str | None, faspr_path: str, pulchra_path: str, input_clashes: list[tuple[str, str]], sampling_targets: dict[str, tuple[tuple[str, tuple[int, Ellipsis], str, str]]], valid_pdbs_dir: str, goal_ensemble_size: int) -> bool | None

   Repack the side-chains and check for steric clashes in a sampled .pdb structure.

   Side-chain repacking is done by passing the structure through FASPR.
   The resulting .pdb file is then rewritten into a single chain with sequential residue numbering
   before being passed into PULCHRA, as it does not support multi-chain structures.
   Clash checking is done by passing the structure through PULCHRA and checking its output.
   If no clashes are present, the resulting .pdb file has its chain and residue numbering
   information restored to its original status.

   :param sampled_pdb: Sampled .pdb structure, unprocessed.
   :type sampled_pdb: :py:class:`str, optional`
   :param faspr_path: Path to FASPR executable or its alias.
   :type faspr_path: :py:class:`str`
   :param pulchra_path: Path to PULCHRA executable or its alias.
   :type pulchra_path: :py:class:`str`
   :param input_clashes: List of clashes present in the input structure that, if present, will be ignored.
   :type input_clashes: :py:class:`list[tuple[str,str]]`
   :param sampling_targets: Mapping of chain letters to target regions for sampling.
   :type sampling_targets: :py:class:`dict[str,tuple[tuple[str,tuple[int,...],str,str]]]`
   :param valid_pds_dir: Path to directory where valid structures will be output.
   :type valid_pds_dir: :py:class:`str`
   :param goal_ensemble_size: If the number of structures in valid pdbs directory is ever greater than this value
                              do not write any more structures into the directory.
   :type goal_ensemble_size: :py:class:`int`

   :returns:     True if the sampled .pdb structure has steric clashes, False otherwise.
                 None if an error occured.
   :rtype: bool | None


