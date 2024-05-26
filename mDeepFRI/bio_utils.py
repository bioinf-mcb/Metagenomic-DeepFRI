import logging
from io import StringIO
from typing import List, Literal, Tuple

import foldcomp
import numpy as np
from biotite.sequence import ProteinSequence
from biotite.structure import get_chains
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import PDBxFile, get_structure

from mDeepFRI.alignment import AlignmentResult
from mDeepFRI.alignment_utils import align_contact_map, pairwise_sqeuclidean

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


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


def extract_residues_coordinates(
        structure_string: str,
        chain: str = "A",
        filetype: Literal["mmcif", "pdb"] = "mmcif") -> Tuple[str, np.ndarray]:
    """
    Extracts residues and coordinates from structural string.
    Automatically processes PDB and mmCIF files.

    Args:
        structure_string (str): Structure file read into string.
        max_seq_len (int): Maximum sequence length.

    Returns:
        Tuple[str, np.ndarray]: Tuple of residues and coordinates.
    """
    if filetype == "mmcif":
        mmcif = PDBxFile.read(StringIO(structure_string))
        structure = get_structure(mmcif, model=1)
    elif filetype == "pdb":
        pdb = PDBFile.read(StringIO(structure_string))
        structure = pdb.get_structure()[0]
    else:
        raise NotImplementedError(f"Filetype {filetype} not supported.")

    residues, coords = get_residues_coordinates(structure, chain=chain)

    return (residues, coords)


def foldcomp_sniff_suffix(idx: str, database_path: str) -> str:
    """
    Sniff suffix for FoldComp database.

    Args:
        idx (str): Protein ID.
        database_path (str): Path to FoldComp database.

    Returns:
        str: Suffix for the database ids.
    """
    with foldcomp.open(database_path, ids=[idx]) as db:
        for _, structure in db:
            suffix = None
    if "structure" not in locals():
        idx = idx + ".pdb"
        with foldcomp.open(database_path, ids=[idx]) as db:
            for _, structure in db:
                suffix = ".pdb"

    return suffix


def get_foldcomp_structures(ids: List[str], database_path: str) -> List[str]:
    """
    Retrieves structure either from PDB or supplied FoldComp database.
    Extracts sequence and coordinate infromaton.

    Args:
        ids (List[str]): List of protein ids.
        database_path (str): Path to FoldComp database. If empty, the structure will be retrieved from the PDB.

    Returns:
        List[str]: List of structures.
    """
    structures = []
    with foldcomp.open(database_path, ids=ids) as db:
        for _, pdb in db:
            structures.append(pdb)

    return structures


def build_align_contact_map(
        alignment: AlignmentResult,
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
    idx = alignment.target_name.rsplit(".", 1)[0]
    coordinates = alignment.coords
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
        except IndexError:
            pdb_id, chain = idx.upper().split("_")
            logger.warning(
                f"Error aligning contact map for PDB ID {pdb_id}[Chain {chain}] "
                f"against {alignment.query_name}.")
            aligned_cmap = None

    else:
        aligned_cmap = None

    return (alignment, aligned_cmap)
