from io import StringIO
from typing import Dict, List, Tuple

import foldcomp
import numpy as np
from Bio.PDB import PDBParser
from pysam import FastaFile, FastxFile
from scipy.spatial.distance import pdist, squareform


def load_fasta_as_dict(fasta_file):
    """Load FASTA file as dict"""
    with FastxFile(fasta_file) as fasta:
        fasta_dict = {entry.name: entry.sequence for entry in fasta}
    return fasta_dict


def retrieve_fasta_entries_as_dict(fasta_file, entries):
    """Retrieve selected FASTA entries as dict"""
    fasta_dict = dict()
    with FastaFile(fasta_file) as fasta:
        for name in entries:
            fasta_dict[name] = fasta.fetch(name)
    return fasta_dict


# copied from Biopython to remove dependency
protein_letters_1to3 = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}
protein_letters_1to3_extended = {
    **protein_letters_1to3,
    **{
        "B": "Asx",
        "X": "Xaa",
        "Z": "Glx",
        "J": "Xle",
        "U": "Sec",
        "O": "Pyl"
    },
}

protein_letters_3to1_extended = {
    value: key
    for key, value in protein_letters_1to3_extended.items()
}

PROTEIN_LETTERS = dict()
for k, v in protein_letters_3to1_extended.items():
    PROTEIN_LETTERS[str.upper(k)] = v
PROTEIN_LETTERS["UNK"] = "X"


def seq2onehot(seq):
    """Create 26-dim embedding"""
    chars = [
        '-', 'D', 'G', 'U', 'L', 'N', 'T', 'K', 'H', 'Y', 'W', 'C', 'P', 'V',
        'S', 'O', 'I', 'E', 'F', 'X', 'Q', 'A', 'B', 'Z', 'R', 'M'
    ]
    vocab_size = len(chars)
    vocab_embed = dict(zip(chars, range(vocab_size)))

    # Convert vocab to one-hot
    vocab_one_hot = np.zeros((vocab_size, vocab_size), int)
    for _, val in vocab_embed.items():
        vocab_one_hot[val, val] = 1

    embed_x = [vocab_embed[v] for v in seq]
    seqs_x = np.array([vocab_one_hot[j, :] for j in embed_x])

    return seqs_x


def calculate_cmap(pdb_string, max_seq_len=1000, threshold=6.0, mode="matrix"):

    parser = PDBParser()
    structure = parser.get_structure("", StringIO(pdb_string))

    residues = [x for x in structure.get_residues()][:max_seq_len]

    coords = np.array([residue["CA"].get_coord() for residue in residues])
    distances = squareform(pdist(coords))

    cmap = (distances < threshold).astype(np.int32)

    if mode == "sparse":
        cmap = np.argwhere(cmap == 1).astype(np.uint32)
    else:
        pass

    return cmap


def retrieve_structure(idx, db):
    with foldcomp.open(db, ids=[idx]) as db:
        for _, pdb in db:
            structure = pdb

        return structure


def align_contact_map(query_alignment: str,
                      target_alignment: str,
                      sparse_target_contact_map: List[Tuple[int, int]],
                      generated_contacts: int = 2) -> np.ndarray:
    """
    Aligns a contact map based on the alignments of query and target sequences.

    Args:
        query_alignment: The alignment of the query sequence.
        target_alignment: The alignment of the target sequence.
        sparse_target_contact_map: The sparse contact map of the target
                                   sequence represented as a list of tuples (i, j)
                                   indicating contacts between residues iand j.
        generated_contacts: The number of generated contacts to add for gapped
                            regions in the query alignment. Defaults to 2.

    Returns:
        The aligned contact map as a numpy array.

    Algorithm:
    1. Initialize an empty list `sparse_query_contact_map` to store the contacts in the aligned contact map.
    2. Initialize variables `target_index` and `query_index` to track the indices of residues in the target
    and query proteins, respectively.
    3. Initialize an empty dictionary `target_to_query_indices` to map target residues to query residues
    using shift resulting from the alignments.
    4. Iterate over each position in the query alignment:
        - If the query residue is '-', increment the `target_index` and do not add a contact
        to the aligned contact map.
        - If the query residue is not '-', check the target residue:
            - If the target residue is '-', add contacts for the generated region in the query alignment:
                - For each generated contact, add the contact (query_index + j, query_index)
                and (query_index - j, query_index) to the `sparse_query_contact_map` according to generated_contacts.
                - Increment the `query_index`.
            - If the target residue is not '-', map the target residue to the query residue by adding
            an entry in the `target_to_query_indices` dictionary.
                - Increment both the `query_index` and `target_index`.
    5. Translate the target residue indices to query residue indices
    in the `sparse_target_contact_map` by using the `target_to_query_indices` dictionary.
    6. Filter out the contacts that are not present in the query alignment by removing contacts
    with '-1' indices from the `sparse_target_contact_map`.
    7. Add the filtered contacts from the filtered `sparse_target_contact_map` to the `sparse_query_contact_map`.
    8. Build the output contact map with dimensions (query_index, query_index) initialized as all zeros.
    Query index is the number of residues in the query sequence.
    9. Set the diagonal elements of the output contact map to 1.
    10. For each contact (i, j) in the `sparse_query_contact_map`:
        - If i is less than 0 or greater than or equal to `query_index`, skip the contact.
        - Otherwise, set the corresponding elements in the output contact map to 1 symmetrically.
    11. Return the aligned contact map as a numpy array.
    """
    # The sparse contact map of the query sequence will contain all contacts. Will be used to create dense contact map
    sparse_query_contact_map: List[Tuple[int, int]] = []

    # The index of the residues in sequences
    target_index: int = 0
    query_index: int = 0

    # Map target residues to query residues based on the alignments
    target_to_query_indices: Dict[int, int] = {}

    # Map target residues to query residues based on the alignments
    for i in range(len(query_alignment)):
        # If the query residue is a gap, skip target residue
        if query_alignment[i] == "-":
            target_to_query_indices[target_index] = -1
            target_index += 1
        else:
            # If the target residue is a gap, add contacts to the query residue
            # connected to generated_contacts nearest residues
            if target_alignment[i] == "-":
                for j in range(1, generated_contacts + 1):
                    sparse_query_contact_map.append(
                        (query_index + j, query_index))
                    sparse_query_contact_map.append(
                        (query_index - j, query_index))
                query_index += 1
            else:
                # If there is an alignment match, map target residue to query residue
                target_to_query_indices[target_index] = query_index
                query_index += 1
                target_index += 1

    # Translate the target residues index to query residues index
    sparse_map = list(
        map(
            lambda x:
            (target_to_query_indices[x[0]], target_to_query_indices[x[1]]),
            sparse_target_contact_map))
    # Filter out the contacts that are not in the query alignment by removing columns and rows from sparse contact map
    sparse_map = list(filter(lambda x: x[0] != -1 and x[1] != -1, sparse_map))
    # Add the contacts to the output contact map
    sparse_query_contact_map.extend(sparse_map)

    # Build the output contact map
    output_contact_map = np.zeros((query_index, query_index))
    # Fill the diagonal
    for i in range(query_index):
        output_contact_map[i, i] = 1
    # Fill the contacts from the sparse query contact map
    for i, j in sparse_query_contact_map:
        if i < 0:
            continue
        if i >= query_index:
            continue
        # Apply symmetry
        output_contact_map[i, j] = 1
        output_contact_map[j, i] = 1

    return output_contact_map
