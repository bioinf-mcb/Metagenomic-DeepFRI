import logging
from io import StringIO
from typing import Dict, List, Tuple

import foldcomp
import numpy as np
import pysam
import requests
from biotite.sequence import ProteinSequence
from biotite.structure import get_chains
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import PDBxFile, get_structure
from pysam import FastaFile, FastxFile

from mDeepFRI.alignment_utils import (align_contact_map, alignment_identity,
                                      pairwise_sqeuclidean)

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


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
    # silence pysam warnings for duplicate sequences
    verb = pysam.set_verbosity(0)

    with FastaFile(fasta_file) as fasta:
        for name in entries:
            fasta_dict[name] = fasta.fetch(name)

    # reset verbosity
    pysam.set_verbosity(verb)

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
        cmap = np.argwhere(cmap == 1).astype(np.int32)
    else:
        pass

    return cmap


def get_residues_coordinates(structure: np.ndarray, chain: str = "A"):
    """
    Retrieves residues and coordinates from biotite structure.

    Args:
        structure (np.ndarray): Structure file read into string.

    Returns:
        Tuple[str, np.ndarray]: Tuple of residues and coordinates.
    """

    chains = get_chains(structure)
    if chain not in chains:
        raise ValueError(f"Chain {chain} not found in structure.")

    protein_chain = structure[structure.chain_id == chain]
    # extract CA atoms coordinates
    ca_atoms = protein_chain[protein_chain.atom_name == "CA"]

    residues = str(ProteinSequence(ca_atoms.res_name))
    coords = ca_atoms.coord

    return (residues, coords)


def extract_residues_coordinates(structure_string: str,
                                 chain: str = "A") -> Tuple[str, np.ndarray]:
    """
    Extracts residues and coordinates from structural string.
    Automatically processes PDB and mmCIF files.

    Args:
        structure_string (str): Structure file read into string.
        max_seq_len (int): Maximum sequence length.

    Returns:
        Tuple[str, np.ndarray]: Tuple of residues and coordinates.
    """
    try:
        mmcif = PDBxFile.read(StringIO(structure_string))
        structure = get_structure(mmcif, model=1)
    except (UnboundLocalError, ValueError):
        pdb = PDBFile.read(StringIO(structure_string))
        structure = pdb.get_structure()[0]

    residues, coords = get_residues_coordinates(structure, chain=chain)

    return (residues, coords)


def retrieve_structure_features(idx: str,
                                database_path: str = "pdb100"
                                ) -> Tuple[str, np.ndarray]:
    """
    Retrieves structure either from PDB or supplied FoldComp database.
    Extracts sequence and coordinate infromaton.

    Args:
        idx (str): Index of structure.
        database_path (str): Path to FoldComp database. If empty, the structure will be retrieved from the PDB.
        max_seq_len (int): Maximum sequence length.

    Returns:
        Tuple[str, np.ndarray]: Tuple of residues and coordinates.
    """
    pdb_http = "https://files.rcsb.org/view/{pdb_id}.cif"

    if database_path == "pdb100":
        pdb_id, chain = idx.split("_")
        pdb_id = pdb_id.lower()
        url = pdb_http.format(pdb_id=pdb_id)
        structure = requests.get(url).text

    else:
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

    residues, coords = extract_residues_coordinates(structure)

    return (residues, coords)


def retrieve_align_contact_map(
        alignment: AlignmentResult,
        database: str = "pdb100",
        threshold: float = 6,
        generated_contacts: int = 2) -> Tuple[AlignmentResult, np.ndarray]:
    """
    Retrieve contact map for aligned sequences.

    Args:
        alignment (AlignmentResult): Alignment of query and target sequences.
        database (str): Path to FoldComp database. If empty, the structure will be retrieved from the PDB.
        threshold (float): Distance threshold for contact map.
        generated_contacts (int): Number of generated contacts to add for gapped regions in the query alignment.

    Returns:
        Tuple[AlignmentResult, np.ndarray]: Tuple of alignment and contact map.
    """

    if "pdb" in str(database):
        database = "pdb100"

    idx = alignment.target_name.rsplit(".", 1)[0]
    try:
        coordinates = retrieve_structure_features(idx,
                                                  database_path=database)[1]

    # catch error for non-standard aminoacids in database
    except KeyError as e:
        coordinates = None
        pdb_id, chain = idx.upper().split("_")

        logger.warning(
            f"Error extracting residues and coordinates for PDB ID {pdb_id}[Chain {chain}] - "
            f"non-standard residue {str(e)} present; {alignment.query_name} alignment skipped."
        )

    # output of the function might None
    if coordinates is not None:

        cmap = calculate_contact_map(coordinates,
                                     threshold=threshold,
                                     mode="sparse")
        logger.debug("Aligning contact map for %s against %s.", idx,
                     alignment.query_name)
        try:
            aligned_cmap = align_contact_map(alignment.gapped_sequence,
                                             alignment.gapped_target, cmap,
                                             generated_contacts)
        except KeyError:
            pdb_id, chain = idx.upper().split("_")
            logger.warning(
                f"Error aligning contact map for PDB ID {pdb_id}[Chain {chain}] "
                f"against {alignment.query_name}.")
            aligned_cmap = None

    else:
        aligned_cmap = None

    return (alignment, aligned_cmap)
