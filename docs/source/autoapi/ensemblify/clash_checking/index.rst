ensemblify.clash_checking
=========================

.. py:module:: ensemblify.clash_checking

.. autoapi-nested-parse::

   Simple module to re-check previously generated ensembles for steric clashes.



Functions
---------

.. autoapisummary::

   ensemblify.clash_checking.process_pulchra_output
   ensemblify.clash_checking.check_report_pdb_clashes
   ensemblify.clash_checking.check_steric_clashes


Module Contents
---------------

.. py:function:: process_pulchra_output(sampled_pdb: str, pulchra_output_buffer: str, sampling_targets: dict[str, tuple[tuple[str, tuple[int, Ellipsis], str, str]]] | None = None, input_clashes: list[tuple[str, str]] | None = None) -> list[str]

   Check if there are recorded steric clashes in given PULCHRA output.

   Clashes present in input structure (if provided) are not ignored.
   Clashes are only considered when at least one residue belongs to a sampled region (if those
   regions are provided).

   :param sampled_pdb: Path to .pdb file output from conformational sampling.
   :type sampled_pdb: :py:class:`str`
   :param pulchra_output_buffer: Stdout from applying PULCHRA to the sampled .pdb structure.
   :type pulchra_output_buffer: :py:class:`str`
   :param sampling_targets: Mapping of chain identifiers to sampled residue numbers.
   :type sampling_targets: :py:class:`dict[str,tuple[tuple[str,tuple[int,...],str,str]]] | None`
   :param input_clashes: Clashes present in the sampling input structure, that will be ignored if
                         present in the given PULCHRA output.
   :type input_clashes: :py:class:`list[tuple[str,str]] | None`

   :returns:     List of clashes present in sampled .pdb file.
   :rtype: list[str]


.. py:function:: check_report_pdb_clashes(pdb2check: str, sampling_targets: dict[str, tuple[tuple[str, tuple[int, Ellipsis], str, str]]] | None = None, input_clashes: list[tuple[str, str]] | None = None) -> tuple[str, list[str] | None]

   Check for steric clashes in a .pdb file, optionally considering sampling targets and clashes
   in the input structure.

   A steric clash is reported when the distance between any two non bonded atoms is less than two
   angstrom.

   :param pdb2check: Path to .pdb file to check for steric clashes.
   :type pdb2check: :py:class:`str`
   :param sampling_targets: Mapping of chains to sampled regions following Ensemblify parameters style. If
                            provided, clashes are only checked for in these regions. Defaults to None.
   :type sampling_targets: :py:class:`dict[str,tuple[tuple[str,tuple[int,...],str,str]]] | None`
   :param input_clashes: List of clashes detected in the ensemble generation input structure. If provided,
                         clashes detailed here will be ignored if found in the .pdb to check (only in sampled
                         regions, if sampling targets is provided). Defaults to None.
   :type input_clashes: :py:class:`list[tuple[str,str]] | None`

   :returns:

                 pdb2check (str):
                     Path to the sampled .pdb to check for clashes.
                 steric_clashes (list[str] | None):
                     List of clashes present in sampled .pdb file or None, if PULCHRA erred.
   :rtype: tuple[str,list[str] | None]


.. py:function:: check_steric_clashes(ensemble_dir: str, sampling_targets: str | dict[str, tuple[tuple[str, tuple[int, Ellipsis], str, str]]] | None = None, input_structure: str | None = None) -> tuple[str, str]

   Check a generated ensemble for steric clashes, outputting clash reports.

   A directory is created inside the ensemble directory where clash reports (simple and detailed)
   will be stored, as well as any files output by processing the input structure (if provided).

   :param ensemble_dir: Path to directory where ensemble .pdb structures are stored.
   :type ensemble_dir: :py:class:`str`
   :param sampling_targets: Mapping of chains to sampled regions following Ensemblify parameters style or path to
                            this mapping in .yaml format. If a path is provided, it can either be solely the
                            mapping or a full Ensemblify parameters file. If provided, clashes are only checked
                            for in these regions. Defaults to None.
   :type sampling_targets: :py:class:`str | dict[str,tuple[tuple[str,tuple[int,...],str,str]]] | None`
   :param input_structure: Path to input structure used to generate the ensemble. If provided, steric clashes
                           present in this structure (only in sampled regions, if sampling targets is provided)
                           are ignored if they are detected in any of the sampled structures. Defaults to None.
   :type input_structure: :py:class:`str | None`

   :returns:

                 clash_report (str):
                     Path to file with simplified ensemble clash report, i.e. total number of clashed
                     structures and how many clashes were detected in each structure.
                 clash_report_detailed (str):
                     Path to file with detailed ensemble clash report, i.e. how many clashes were
                     detected in each structure and the atoms involved in the detected clash.
   :rtype: tuple[str,str]


