from io import StringIO
from typing import Dict, List, Tuple

import foldcomp
import numpy as np
import requests
from biotite.sequence import ProteinSequence
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import PDBxFile, get_structure
from pysam import FastaFile, FastxFile

from mDeepFRI.alignment_utils import alignment_identity, pairwise_sqeuclidean


class AlignmentResult:
    """
    Class for storing pairwise alignment results.

    Attributes:
        query_name (str): Name of the query sequence.
        query_sequence (str): Query sequence.
        target_name (str): Name of the target sequence.
        target_sequence (str): Target sequence.
        alignment (str): Alignment string.
        gapped_sequence (str): Query sequence with gaps.
        gapped_target (str): Target sequence with gaps.
        identity (float): Identity between two sequences.

    Methods:
        insert_gaps: Inserts gaps into query and target sequences.
        calculate_identity: Calculates identity between query and target.
    """
    def __init__(self, query_name, query_sequence, target_name,
                 target_sequence, alignment):

        self.query_name = query_name
        self.query_sequence = query_sequence
        self.target_name = target_name
        self.target_sequence = target_sequence
        self.alignment = alignment
        self.insert_gaps()
        self.calculate_identity()
        self.db_name = None

    def __str__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"identity={self.identity})"

    def __repr__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"identity={self.identity})"

    def insert_gaps(self):
        """
        Inserts gaps into query and target sequences.
        """

        self.gapped_sequence, self.gapped_target = insert_gaps(
            self.query_sequence, self.target_sequence, self.alignment)

        return self

    def calculate_identity(self):
        """
        Calculates identity between query and target.

        Args:
            None
        Returns:
            self
        """

        self.identity = alignment_identity(self.gapped_sequence,
                                           self.gapped_target)

        return self


def seq2onehot(seq: str) -> np.ndarray:
    """Create 26-dim one-hot encoding of a protein sequence.

    Args:
        seq (str): Protein sequence.

    Returns:
        np.ndarray: One-hot encoding of protein sequence.
    """

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


def insert_gaps(sequence: str, reference: str,
                alignment_string: str) -> Tuple[str, str]:
    """
    Inserts gaps into query and target sequences.

    Args:
        sequence (str): Query sequence.
        reference (str): Target sequence.
        alignment_string (str): Alignment string.

    Returns:
        gapped_sequence (str): Query sequence with gaps.
        gapped_target (str): Target sequence with gaps.
    """

    sequence = list(sequence)
    reference = list(reference)
    alignment_string = list(alignment_string)

    for i, a in enumerate(alignment_string):
        if a == "I":
            sequence.insert(i, "-")
        elif a == "D":
            reference.insert(i, "-")
    return "".join(sequence), "".join(reference)


def load_fasta_as_dict(fasta_file: str) -> Dict[str, str]:
    """
    Load FASTA file as dict of headers to sequences

    Args:
        fasta_file (str): Path to FASTA file. Can be compressed.

    Returns:
        Dict[str, str]: Dictionary of FASTA entries sorted by length.
    """

    with FastxFile(fasta_file) as fasta:
        fasta_dict = {entry.name: entry.sequence for entry in fasta}

    return fasta_dict


def retrieve_fasta_entries_as_dict(fasta_file: str,
                                   entries: List[str]) -> Dict[str, str]:
    """
    Retrieve selected FASTA entries as dict

    Args:
        fasta_file (str): Path to FASTA file. Can be compressed.
        entries (List[str]): List of entries to retrieve.

    Returns:
        Dict[str, str]: Dictionary of FASTA entries.
    """

    fasta_dict = dict()
    with FastaFile(fasta_file) as fasta:
        for name in entries:
            fasta_dict[name] = fasta.fetch(name)
    return fasta_dict


def calculate_contact_map(coordinates: np.ndarray,
                          threshold=6.0,
                          mode="matrix") -> np.ndarray:
    """
    Calculate contact map from PDB string.

    Args:
        pdb_string (str): PDB file read into string.
        max_seq_len (int): Maximum sequence length.
        threshold (float): Distance threshold for contact map.
        mode (str): Output mode. Either "matrix" or "sparse".

    Returns:
        np.ndarray: Contact map.
    """
    # squared euclidean is used for efficiency
    threshold = threshold**2

    distances = pairwise_sqeuclidean(coordinates)
    cmap = (distances < threshold).astype(np.int32)

    if mode == "sparse":
        cmap = np.argwhere(cmap == 1).astype(np.uint32)
    else:
        pass

    return cmap


def extract_residues_coordinates(structure_string: str,
                                 chain: str = None) -> Tuple[List, np.ndarray]:
    """
    Extracts residues and coordinates from structural file.
    Automatically processes PDB and mmCIF files.

    Args:
        structure_string (str): Structure file read into string.
        max_seq_len (int): Maximum sequence length.

    Returns:
        Tuple[List, np.ndarray]: Tuple of residues and coordinates.
    """
    try:
        mmcif = PDBxFile.read(StringIO(structure_string))
        structure = get_structure(mmcif, model=1)
    except UnboundLocalError:
        pdb = PDBFile(StringIO(structure_string))
        structure = pdb.get_structure()[0]

    protein_chain = structure[structure.chain_id == chain]
    # extract CA atoms coordinates
    ca_atoms = protein_chain[protein_chain.atom_name == "CA"]
    residues = str(ProteinSequence(ca_atoms.res_name))
    coords = ca_atoms.coord

    return (residues, coords)


def retrieve_structure_features(idx: str,
                                database_path: str = "pdb100"
                                ) -> Tuple[List, np.ndarray]:
    """
    Retrieves structure either from PDB or supplied FoldComp database.
    Extracts sequence and coordinate infromaton.

    Args:
        idx (str): Index of structure.
        database_path (str): Path to FoldComp database. If empty, the structure will be retrieved from the PDB.
        max_seq_len (int): Maximum sequence length.

    Returns:
        Tuple[List, np.ndarray]: Tuple of residues and coordinates.
    """
    pdb_http = "https://files.rcsb.org/view/{pdb_id}.cif"

    if database_path == "pdb100":
        pdb_id, chain = idx.split("_")
        pdb_id = pdb_id.lower()
        url = pdb_http.format(pdb_id=pdb_id)
        structure = requests.get(url).text

    else:
        chain = None
        with foldcomp.open(database_path, ids=[idx]) as db:
            for _, pdb in db:
                structure = pdb

        # issue with FoldComp inconsistency
        # https://github.com/steineggerlab/foldcomp/issues/45
        if "structure" not in locals():
            idx = idx + ".pdb"
            with foldcomp.open(database_path, ids=[idx]) as db:
                for _, pdb in db:
                    structure = pdb

    residues, coords = extract_residues_coordinates(structure, chain=chain)

    return (residues, coords)


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


def retrieve_align_contact_map(
        alignment: AlignmentResult,
        database: str = "pdb100",
        max_seq_len: int = 1000,
        threshold: float = 6,
        generated_contacts: int = 2) -> Tuple[AlignmentResult, np.ndarray]:
    """
    Retrieve contact map for aligned sequences.

    Args:
        alignment (AlignmentResult): Alignment of query and target sequences.
        database (str): Path to FoldComp database. If empty, the structure will be retrieved from the PDB.
        max_seq_len (int): Maximum sequence length.
        threshold (float): Distance threshold for contact map.
        generated_contacts (int): Number of generated contacts to add for gapped regions in the query alignment.

    Returns:
        Tuple[AlignmentResult, np.ndarray]: Tuple of alignment and contact map.
    """

    if "pdb" in str(database):
        database = "pdb100"

    idx = alignment.target_name.rsplit(".", 1)[0]
    coordinates = retrieve_structure_features(idx, database_path=database)[1]

    cmap = calculate_contact_map(coordinates,
                                 threshold=threshold,
                                 mode="sparse")

    aligned_cmap = align_contact_map(alignment.gapped_sequence,
                                     alignment.gapped_target, cmap,
                                     generated_contacts)

    return (alignment, aligned_cmap)
