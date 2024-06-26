import logging
import sys
from io import StringIO
from typing import List, Literal, Tuple

import foldcomp
import numpy as np
from biotite.sequence import ProteinSequence
from biotite.structure import get_chains
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import PDBxFile, get_structure

from mDeepFRI.alignment import AlignmentResult
from mDeepFRI.contact_map_utils import align_contact_map, pairwise_sqeuclidean

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    '[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# https://github.com/openmm/pdbfixer/blob/master/pdbfixer/pdbfixer.py
substitutions = {
    '2AS': 'ASP',
    '3AH': 'HIS',
    '5HP': 'GLU',
    '5OW': 'LYS',
    'ACL': 'ARG',
    'AGM': 'ARG',
    'AIB': 'ALA',
    'ALM': 'ALA',
    'ALO': 'THR',
    'ALY': 'LYS',
    'ARM': 'ARG',
    'ASA': 'ASP',
    'ASB': 'ASP',
    'ASK': 'ASP',
    'ASL': 'ASP',
    'ASQ': 'ASP',
    'AYA': 'ALA',
    'BCS': 'CYS',
    'BHD': 'ASP',
    'BMT': 'THR',
    'BNN': 'ALA',
    'BUC': 'CYS',
    'BUG': 'LEU',
    'C5C': 'CYS',
    'C6C': 'CYS',
    'CAS': 'CYS',
    'CCS': 'CYS',
    'CEA': 'CYS',
    'CGU': 'GLU',
    'CHG': 'ALA',
    'CLE': 'LEU',
    'CME': 'CYS',
    'CSD': 'ALA',
    'CSO': 'CYS',
    'CSP': 'CYS',
    'CSS': 'CYS',
    'CSW': 'CYS',
    'CSX': 'CYS',
    'CXM': 'MET',
    'CY1': 'CYS',
    'CY3': 'CYS',
    'CYG': 'CYS',
    'CYM': 'CYS',
    'CYQ': 'CYS',
    'DAH': 'PHE',
    'DAL': 'ALA',
    'DAR': 'ARG',
    'DAS': 'ASP',
    'DCY': 'CYS',
    'DGL': 'GLU',
    'DGN': 'GLN',
    'DHA': 'ALA',
    'DHI': 'HIS',
    'DIL': 'ILE',
    'DIV': 'VAL',
    'DLE': 'LEU',
    'DLY': 'LYS',
    'DNP': 'ALA',
    'DPN': 'PHE',
    'DPR': 'PRO',
    'DSN': 'SER',
    'DSP': 'ASP',
    'DTH': 'THR',
    'DTR': 'TRP',
    'DTY': 'TYR',
    'DVA': 'VAL',
    'EFC': 'CYS',
    'FLA': 'ALA',
    'FME': 'MET',
    'GGL': 'GLU',
    'GL3': 'GLY',
    'GLZ': 'GLY',
    'GMA': 'GLU',
    'GSC': 'GLY',
    'HAC': 'ALA',
    'HAR': 'ARG',
    'HIC': 'HIS',
    'HIP': 'HIS',
    'HMR': 'ARG',
    'HPQ': 'PHE',
    'HTR': 'TRP',
    'HYP': 'PRO',
    'IAS': 'ASP',
    'IIL': 'ILE',
    'IYR': 'TYR',
    'KCX': 'LYS',
    'LLP': 'LYS',
    'LLY': 'LYS',
    'LTR': 'TRP',
    'LYM': 'LYS',
    'LYZ': 'LYS',
    'MAA': 'ALA',
    'MEN': 'ASN',
    'MHS': 'HIS',
    'MIS': 'SER',
    'MK8': 'LEU',
    'MLE': 'LEU',
    'MPQ': 'GLY',
    'MSA': 'GLY',
    'MSE': 'MET',
    'MVA': 'VAL',
    'NEM': 'HIS',
    'NEP': 'HIS',
    'NLE': 'LEU',
    'NLN': 'LEU',
    'NLP': 'LEU',
    'NMC': 'GLY',
    'OAS': 'SER',
    'OCS': 'CYS',
    'OMT': 'MET',
    'PAQ': 'TYR',
    'PCA': 'GLU',
    'PEC': 'CYS',
    'PHI': 'PHE',
    'PHL': 'PHE',
    'PR3': 'CYS',
    'PRR': 'ALA',
    'PTR': 'TYR',
    'PYX': 'CYS',
    'SAC': 'SER',
    'SAR': 'GLY',
    'SCH': 'CYS',
    'SCS': 'CYS',
    'SCY': 'CYS',
    'SEL': 'SER',
    'SEP': 'SER',
    'SET': 'SER',
    'SHC': 'CYS',
    'SHR': 'LYS',
    'SMC': 'CYS',
    'SOC': 'CYS',
    'STY': 'TYR',
    'SVA': 'SER',
    'TIH': 'ALA',
    'TPL': 'TRP',
    'TPO': 'THR',
    'TPQ': 'ALA',
    'TRG': 'LYS',
    'TRO': 'TRP',
    'TYB': 'TYR',
    'TYI': 'TYR',
    'TYQ': 'TYR',
    'TYS': 'TYR',
    'TYY': 'TYR'
}


def calculate_contact_map(coordinates: np.ndarray,
                          threshold=6.0,
                          distance="sqeuclidean",
                          mode="matrix") -> np.ndarray:
    """
    Calculate contact map from PDB string.

    Args:
        pdb_string (str): PDB file read into string.
        threshold (float): Distance threshold for contact map.
        mode (str): Output mode. Either "matrix" or "sparse".

    Returns:
        np.ndarray: Contact map.
    """
    # squared euclidean is used for efficiency
    distance_functions = {"sqeuclidean": pairwise_sqeuclidean}

    if distance == "sqeuclidean":
        threshold = threshold**2

    distance_func = distance_functions[distance]

    distances = distance_func(coordinates)
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
    ca_atoms = protein_chain[(protein_chain.atom_name == "CA")
                             & (protein_chain.hetero == False)]  # noqa

    residues = str(
        ProteinSequence(
            [substitutions.get(res, res) for res in ca_atoms.res_name]))
    coords = ca_atoms.coord

    return (residues, coords)


def load_structure(structure_string: str,
                   filetype: Literal["mmcif", "pdb"] = "mmcif") -> np.ndarray:
    """
    Load structure from string.

    Args:
        structure_string (str): Structure file read into string.
        filetype (str): Filetype of the structure.

    Returns:
        np.ndarray: Structure.
    """
    if filetype == "mmcif":
        mmcif = PDBxFile.read(StringIO(structure_string))
        structure = get_structure(mmcif, model=1)
    elif filetype == "pdb":
        pdb = PDBFile.read(StringIO(structure_string))
        structure = pdb.get_structure()[0]
    else:
        raise NotImplementedError(f"Filetype {filetype} not supported.")

    return structure


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

    structure = load_structure(structure_string, filetype=filetype)
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
        logger.warning(f"No coordinates found for {alignment.target_name}.")
        aligned_cmap = None

    return (alignment, aligned_cmap)
