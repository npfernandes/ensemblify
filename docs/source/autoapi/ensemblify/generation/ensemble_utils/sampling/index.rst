ensemblify.generation.ensemble_utils.sampling
=============================================

.. py:module:: ensemblify.generation.ensemble_utils.sampling

.. autoapi-nested-parse::

   Generate conformational ensembles by sampling regions in an input structure.



Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.sampling.run_sampling


Module Contents
---------------

.. py:function:: run_sampling(input_parameters: str, input_clashes_file: str | None, valid_pdbs_dir: str, sampling_log: str)

   Perform conformational sampling according to input parameters.

   The Pulchra clashes file for the input structure must be provided, along with the path to
   the directory where sampled structures will be stored and a sampling log file.
   Additional log files will be created in the same directory as the provided log file.

   :param input_parameters: Path to parameters file following the Ensemblify template.
   :type input_parameters: :py:class:`str`
   :param input_clashes_file: Path to Pulchra log file for the input structure.
   :type input_clashes_file: :py:class:`str`
   :param valid_pdbs_dir: Path to the directory where sampled structures will be stored.
   :type valid_pdbs_dir: :py:class:`str`
   :param sampling_log: Path to the .log file for the sampling process.
   :type sampling_log: :py:class:`str`


