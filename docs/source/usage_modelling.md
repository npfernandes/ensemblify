# The `modelling` module
  
With the `modelling` module, you can create a full-length structure of your protein of interest, 
by fusing together the structures of folded domains with disordered tails/inter-domain linkers.
This structure can then be used as an input starting structure for conformational ensemble generation.

## Create your full-length protein structure

To fuse disordered sequences with PDB structures into a full-length PDB structure, provide Ensemblify with:

- the path(s) to FASTA file(s) containing the sequences of all (folded + disordered) protein domains, from N- to C-terminal. Can be provided either in FASTA or Multi-FASTA format.

- the path(s) to PDB file(s) containing the structures of all folded protein domains, from N- to C-terminal. Can be provided as a single PDB file with multiple MODEL entries.

- the name for the created full-length fused structure (without .pdb extension).

````{tabs}

   ```{code-tab} console CLI
   (ensemblify_env) $ ensemblify modelling -f protein_name_all.fa -p protein_name_all.pdb -id protein_name
   ```

   ```{code-tab} python Python
   from ensemblify.modelling import fuse_structures
   fuse_structures(['protein_name_all.fa'], ['protein_name_all.pdb'], 'protein_name')
   ```
````
