ensemblify.generation.ensemble_utils.inputs_processing
======================================================

.. py:module:: ensemblify.generation.ensemble_utils.inputs_processing

.. autoapi-nested-parse::

   Auxiliary functions to process input parameters file and the input structure.



Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.inputs_processing.process_input_pdb
   ensemblify.generation.ensemble_utils.inputs_processing.register_input_clashes
   ensemblify.generation.ensemble_utils.inputs_processing.get_protein_info
   ensemblify.generation.ensemble_utils.inputs_processing.setup_ensemble_gen_params
   ensemblify.generation.ensemble_utils.inputs_processing.read_input_parameters


Module Contents
---------------

.. py:function:: process_input_pdb(faspr_path: str, pulchra_path: str, inputs_dir: str, input_pdb: str) -> tuple[str, str]

   Process input structure prior to sampling.

   Apply FASPR and Pulchra to the .pdb of our sampling input structure, outputting the resulting
   .pdb file. Allows us to start the sampling from a pdb file with less clashes, while taking note
   of the clashes still present in the input pdb so those are ignored later.

   :param faspr_path: Path to the FASPR executable.
   :type faspr_path: :py:class:`str`
   :param pulchra_path: Path to the PULCHRA executable.
   :type pulchra_path: :py:class:`str`
   :param inputs_dir: Path to directory where files and .logs resulting from input processing will be stored.
   :type inputs_dir: :py:class:`str`
   :param input_pdb: Path to the .pdb input structure to process.
   :type input_pdb: :py:class:`str`

   :returns:

                 clashes_file (str):
                     Path to the .log file with the output of applying PULCHRA to our input structure.
                 processed_pdb (str):
                     Path to the .pdb file resulting from applying FASPR and PULCHRA to our input structure.
   :rtype: tuple[str,str]


.. py:function:: register_input_clashes(input_clashes_file: str | None) -> list[tuple[str, str]]

   Register clashes in input structure to be ignored later.

   :param input_clashes_file: Path to the PULCHRA output file for the input structure. If None, no clashes are
                              registered.
   :type input_clashes_file: :py:class:`str, optional`

   :returns:     List of clashes in the input structure. Can be empty list if input_clashes_file
                 is None. For example:

                 [ ('ARG[277]', 'GLY[287]'),('ARG[706]', 'GLY[716]'), ... ]
   :rtype: list[tuple[str,str]]


.. py:function:: get_protein_info(uniprot_accession: str) -> dict

   Get information about a protein from the AlphaFold Database using a given UniProt accession.

   :param uniprot_accession: UniProt accession to use in request for AlphaFold's Database API.
   :type uniprot_accession: :py:class:`str`

   :returns:     Information about the protein identified by the given UniProt accession, including
                 links to its .pdb structure and .json PAE matrix.
   :rtype: dict

   Adapted from:
       https://github.com/PDBeurope/afdb-notebooks/blob/main/AFDB_API.ipynb


.. py:function:: setup_ensemble_gen_params(input_params: dict, inputs_dir: str) -> tuple[str, str | None]

   Update sampling input parameters, store steric clashes present in input structure.

   If a UniProt accession is given in either the sequence or PAE fields, replace it with the
   corresponding downloaded .pdb or .json file.

   - 'sequence' field is updated with the path to the processed input structure. If a UniProt
     accession is provided, replace it with the corresponding downloaded .pdb file.
   - 'pae' field is updated with the path to the downloaded .json PAE matrix file, if a UniProt
     accession is provided.
   - 'output_path' field is updated with the ensemble directory inside the created directory named
     'job_name'.
   - File with the updated parameters is saved to the inputs directory.
   - Input structure is processed and any steric clashes present after processing are stored in a
     file so they can later be ignored on non-sampled regions of structures resulting from the
     sampling process.


   :param input_params: Parameters following the Ensemblify template.
   :type input_params: :py:class:`dict`
   :param inputs_dir: Path to directory where files and .logs resulting from input processing will be stored.
   :type inputs_dir: :py:class:`str`

   :returns:

                 processed_parameters_path (str):
                     Path to file where updated parameters are stored.
                 input_clashes (str):
                     Path to file with PULCHRA output from processing input .pdb.
   :rtype: tuple[str,str] | None


.. py:function:: read_input_parameters(parameter_path: str) -> dict

   Read input parameters and assert the validity of its contents.

   :param parameter_path: Path to the parameter .yaml file.
   :type parameter_path: :py:class:`str`

   :returns:     Dictionary with all the parameters validated for correct Python types.
   :rtype: dict


