"""Auxiliary functions for reading and aligning sequences from FASTA or PDB files."""

# IMPORTS
## Standard Library Imports
import os

## Third Party Imports
from Bio import Align

# FUNCTIONS
def get_sequence_from_file(fasta_path: str) -> str:
    """Get the single sequence contained in a FASTA file.
    
    Args:
        fasta_path (str):
            Path to the FASTA file.
    
    Returns:
        str:
            The sequence as a single string.
    """
    with open(fasta_path, 'r') as fasta_file:
        lines = fasta_file.readlines()
        sequence = ''
        for line in lines:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence


def get_sequences_from_file(fasta_path: str) -> list[str]:
    """Get the multiple sequences contained in a FASTA file.
    
    Args:
        fasta_path (str):
            Path to the FASTA file.
    
    Returns:
        list[str]:
            A list of sequences, each as a single string.
    """
    sequences = []
    sequence = ''
    with open(fasta_path, 'r') as fasta_file:
        lines = fasta_file.readlines()
        for line in lines:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences


def split_by_model(pdb_path: str) -> list[str]:
    """Split a PDB file into separate files, one for each model.
    
    Args:
        pdb_path (str):
            Path to the multi-model PDB file.
    Returns:
        list[str]:
            A list of file paths, one for each model.
    """
    model_filepaths = []
    with open(pdb_path, 'r') as pdb_file:
        lines = pdb_file.readlines()
        model_lines = []
        model_number = 0
        for line in lines:
            if line.startswith('MODEL'):
                if model_lines:
                    output_file = os.path.join(
                        os.path.dirname(pdb_path),
                        f'{os.path.splitext(os.path.basename(pdb_path))[0]}_mdl{model_number}.pdb'
                    )
                    with open(output_file, 'w') as model_file:
                        model_file.writelines(model_lines[1:-1])
                    model_filepaths.append(output_file)
                model_lines = [line]
                model_number += 1
            else:
                model_lines.append(line)
        if model_number == 0:
            # If no MODEL lines were found, treat the whole file as a single model
            return [pdb_path]
        if model_lines:
            output_file = os.path.join(
                os.path.dirname(pdb_path),
                f'{os.path.splitext(os.path.basename(pdb_path))[0]}_mdl{model_number}.pdb'
            )
            with open(output_file, 'w') as model_file:
                model_file.writelines(model_lines[1:-1])
            model_filepaths.append(output_file)
    return model_filepaths


def setup_sequences_structures(
    input_fastas: list[str],
    input_pdbs: list[str],
    ) -> tuple[list[str], list[str]]:
    """Setup sequence strings and structure filepaths from the provided FASTA and PDB files.
    
    Parameters:
        input_fastas (list):
            List of FASTA file paths.
        input_pdbs (list):
            List of PDB file paths.
    
    Returns:
        tuple[list,list]:
            list:
                A list of sequence strings extracted from the provided FASTA files.
            list:
                A list of PDB filepaths, one for each model extracted from the provided PDB files.
    """
    # So that using len() below does not trick us
    assert type(input_fastas) == list, "input_fastas must be a list of FASTA file paths."
    assert type(input_pdbs) == list, "input_pdbs must be a list of PDB file paths."

    # Setup the FASTA sequences and PDB structures
    sequences = [] # List of sequences in string format
    pdb_structures = [] # List of PDB file paths

    # Get sequences
    if len(input_fastas) == 1:
        sequences = get_sequences_from_file(input_fastas[0]) # Multiple sequences in a single file
    else:	
        for sequence in input_fastas:
            sequences.append(get_sequence_from_file(sequence)) # Single sequence in each file

    # Get structures
    if len(input_pdbs) == 1:
        pdb_structures = split_by_model(input_pdbs[0]) # Multi-model pdb file
    else:
        for pdb in input_pdbs:
            pdb_structures.append(pdb) # Single structure in each file

    return sequences, pdb_structures


def perform_global_alignment(
    target: str,
    query: str,
    ) -> Align.Alignment:
    """Perform global alignment of two sequences.
    
    Args:
        target (str):
            The target sequence to align against.
        query (str):
            The query sequence to align.
    
    Returns:
        Bio.Align.Alignment:
            The alignment object containing the aligned sequences.
    """
    # Setup the aligner
    aligner = Align.PairwiseAligner(mode='global')
    ## defaults: match = 1.0, mismatch = 0.0, gap score = 0.0
    
    # Set the penalty for opening gaps, enough for our goal here
    aligner.open_gap_score = -10 # automatically sets the global alignment algorithm
    ## between the Needleman-Wunsch, Gotoh and Waterman-Smith-Beyer

    # Perform the alignment
    alignments = aligner.align(target, query)

    # Return the alignment with the highest score
    return alignments[0]


def parse_alignment(alignment: Align.Alignment) -> tuple[str, str, str]:
    """Parse the alignment to extract aligned sequences and encoding.
    
    Args:
        alignment (Bio.Align.Alignment):
            The alignment object containing the aligned sequences.
    
    Returns:
        tuple[str, str, str]:
            str:
                The aligned target sequence.
            str:
                The aligned query sequence.
            str:
                The encoding of the alignment, where:
                - '|' indicates a match,
                - '-' indicates a gap,
                - '.' indicates a mutation.
    """
    # Split the alignment into lines
    alignment_lines = str(alignment).split('\n')

    # Valid characters include the 20 standard amino acids and
    # the alignment symbols for match, gap, mutation
    valid_chars = [
        'G', 'A', 'V', 'L', 'I', 'T', 'S', 'M', 'C', 'P', 'F',
        'Y', 'W', 'H', 'K', 'R','D', 'E', 'N', 'Q', '|', '-', '.'
        ]

    # Build full strings for the aligned sequences and encoding
    aligned_target = ''
    aligned_query = ''
    alignment_encoding = ''
    for line in alignment_lines:
        if 'target' in line:
            aligned_target += ''.join([x for x in line.strip() if x in valid_chars])
        elif 'query' in line:
            aligned_query += ''.join([x for x in line.strip() if x in valid_chars])
        else:
            alignment_encoding += ''.join([x for x in line.strip() if x in valid_chars])
   
    return aligned_target, aligned_query, alignment_encoding


def extract_indels_mutations(
    aligned_target: str,
    aligned_query: str,
    alignment_encoding: str,
    ) -> tuple[list, list]:
    """Extract insertions, deletions, and mutations from the aligned sequences.
    
    Args:
        aligned_target (str):
            The aligned target sequence.
        aligned_query (str):
            The aligned query sequence.
        alignment_encoding (str):
            The encoding of the alignment.
    
    Returns:
        tuple[list, list]:
            list:
                A list of insertions and deletions, where each entry is a list with the
                insertion/deletion segment, start position, end position, and type ('I'
                for insertion, 'D' for deletion).
            list:
                A list of mutations, where each entry is a list with the position, query
                amino acid, and target amino acid.
    
    Reference:
        https://github.com/jferrie3/FusionProteinEnsemble/blob/main/FusionProteinModeler.py
    """
    # Identify positions where the query is aligned to the target
    query_seq_numbering = []
    current_query_seg_numbering = 0
    for encoding_idx, encoding in enumerate(alignment_encoding):       
        
        if aligned_query[encoding_idx] != '-': # Query is aligned in this position
            # Numbering in the current continuous aligned query segment
            current_query_seg_numbering += 1

            # Numbering in the overall aligned query sequence
            # (possibly non-continuous in the alignment)
            query_seq_numbering.append(current_query_seg_numbering)
        else:
            # If the query is not aligned in this position, add a placeholder 0 to
            # ensure query_seq_numbering has same length as the alignment encoding,
            # guaranteeing correct indexing later
            query_seq_numbering.append(0)

    # Parse the alignment encodings, store positions of insertions, deletions, and mutations
    indels = []
    mutations = []
    temp_insert = ''
    insert_start = 0
    insert_end = 0
    deletion_start = 0
    deletion_end = 0

    for encoding_idx, encoding in enumerate(alignment_encoding):
        # A gap in the alignment could mean either the start or continuation
        # of either a deletion or insertion
        if encoding == '-':
            if aligned_target[encoding_idx] == '-': # Deletion
                if insert_start != 0: # Insertion ends, reset start and end
                    insert_end = query_seq_numbering[encoding_idx]
                    indel_code = 'I'
                    indels.append([temp_insert, insert_start, insert_end, indel_code])
                    temp_insert = ''
                    insert_start = 0
                    insert_end = 0
                if deletion_start == 0: # Deletion starts
                    deletion_start = query_seq_numbering[encoding_idx]
                    continue
                else:
                    continue
            else: # Insertion
                temp_insert += aligned_target[encoding_idx] # Build insertion segment
                if deletion_start != 0: # Deletion ends, reset start and end
                    deletion_end = query_seq_numbering[encoding_idx] - 1
                    indel_code = 'D'
                    indels.append(['X', deletion_start, deletion_end, indel_code])
                    deletion_start = 0
                    deletion_end = 0
                if insert_start == 0: # Insertion starts
                    if encoding_idx == 0:
                        insert_start = 1
                    else:
                        insert_start = query_seq_numbering[encoding_idx-1]
                    continue
        # If there is no gap in the alignment, it can be either the end of an insertion or
        # deletion, a mutation, or a match
        else:
            if insert_start != 0: # Insertion ends
                insert_end = query_seq_numbering[encoding_idx]
                indel_code = 'I'
                indels.append([temp_insert, insert_start, insert_end, indel_code])
                temp_insert = '' # Reset insertion segment
                insert_start = 0
                insert_end = 0
            elif deletion_start != 0: # Deletion ends
                deletion_end = encoding_idx
                indel_code = 'D'
                indels.append(['X', deletion_start, deletion_end, indel_code])
                deletion_start = 0
                deletion_end = 0
            if encoding == '.': # Mutation occurs
                mutations.append([query_seq_numbering[encoding_idx],
                                  aligned_query[encoding_idx],
                                  aligned_target[encoding_idx]]) 
            if encoding == '|': # Match occurs
                continue
    
    # If there is an open insertion or deletion at the end of the alignment, close it
    if deletion_start != 0:
        deletion_end = max(query_seq_numbering) # Last position where query is aligned
        indel_code = 'D'
        indels.append(['X', deletion_start, deletion_end, indel_code])			
    if insert_start != 0:
        insert_end = max(query_seq_numbering) # Last position where query is aligned
        indel_code = 'I'
        indels.append([temp_insert, insert_start, insert_end, indel_code])

    return indels, mutations
