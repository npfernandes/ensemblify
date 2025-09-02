"""
Modelling - ``ensemblify.modelling``
======================================

:Author(s): Nuno P. Fernandes
:Year: 2024
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains functions for creating and manipulating PyRosetta Pose objects and associated PDB files.

If desired, these can then be used as starting structures for conformational ensemble generation.

Example applications
--------------------

Create a PyRosetta Pose object from a PDB file or sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ensemblify.modelling.setup_pose`` function can be used to initialize a Pose object from
a .pdb file, a .txt/.fa file containing a sequence, or directly from a provided sequence string.

For example, we can create a Pose object for Histatin5 (Hst5), an intrinsically disordered peptide
with 24 aminoacid residues.

Assuming the path to your .pdb file is assigned to a variable named HST5_PDB, you should run:

>>> import ensemblify.modelling as em
>>> em.setup_pose(HST5_PDB)

Create a full-length PDB structure from ordered and disordered fragments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ensemblify.modelling.fuse_structures`` function can be used to fuse PDB structures of
folded protein domains with the sequences of connecting disordered regions.
For this you need: a list of FASTA files, containing the sequences of all protein domains (folded +
disordered) to be included in the final fused structure; a list of PDB files, containing the
structure of folded domains to be fused.

Assuming the paths to your (Multi-)FASTA and (Multi-Model)PDB files are assigned to variables named
FASTAS and PDBS, respectively, you should run:

>>> import ensemblify.modelling as em
>>> em.fuse_structures(FASTAS, PDBS, 'myprotein')

Available Functions
-------------------

``apply_basic_idealize``
    Idealize bond geometry of a range of residues in a Pose object.
``apply_constraints``
    Apply energy constraints to desired regions of a Pose object.
``apply_pae_constraints``
    Apply energy constraints to a Pose created from an AlphaFold structure based on its Predicted
    Aligned Error (PAE) matrix.
``create_fusion_pose``
    Create a full-length fused Pose object from provided sequences and PDB structures.
``fuse_structures``
    Fuse sequences and PDB structures into a full-length PDB structure.
``pack_sidechains``
    Add side-chains to a Pose object using a PackRotamersMover with the given score function.
``process_pdb_structure``
    Apply FASPR and PULCHRA to PDB structure.
``setup_minmover``
    Setup a PyRosetta MinMover object given the necessary parameters.
``setup_pose``
    Initialize a Pose object from a sequence, a .txt file containing the sequence or a PDB file.

"""

from ensemblify.modelling.constraints import apply_constraints, apply_pae_constraints
from ensemblify.modelling.fusion import apply_basic_idealize, create_fusion_pose, fuse_structures, pack_sidechains
from ensemblify.modelling.objects import setup_minmover, setup_pose
from ensemblify.modelling.pdb_processing import process_pdb_structure

__all__ = [
    'apply_basic_idealize',
    'apply_constraints',
    'apply_pae_constraints',
    'create_fusion_pose',
    'fuse_structures',
    'pack_sidechains',
    'process_pdb_structure',
    'setup_minmover',
    'setup_pose'
]
