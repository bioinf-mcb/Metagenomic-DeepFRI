import dataclasses
import gzip
import logging
import pathlib
import traceback

import numpy as np

from mDeepFRI import ATOMS, SEQUENCES
from mDeepFRI.CPP_lib import libAtomDistanceIO  # type: ignore[attr-defined]
from mDeepFRI.structure_files.parsers import parse_mmcif, parse_pdb
from mDeepFRI.utils import bio_utils

# need to parse different type of files? Add a file name pattern with a parser in this dict.
# read structure_files_parsers/README.md for more information about how to create new parser.
PARSERS = {
    '.pdb': parse_pdb,
    '.pdb.gz': parse_pdb,
    '.cif': parse_mmcif,
    '.cif.gz': parse_mmcif,
    '.ent': parse_pdb,
    '.ent.gz': parse_pdb
}


@dataclasses.dataclass
class SeqAtoms:
    protein_id: str
    atom_amino_group: np.ndarray
    positions: np.ndarray
    groups: np.ndarray


def search_structure_files(input_paths: list):
    """

    :param input_paths:
    :return:
    """
    structure_files_paths = dict()
    # Iterate input_paths and check if file extension matches any key in structure_files.PARSERS
    for input_path in input_paths:
        print(f"{str(input_path)}")
        if input_path.is_file():
            for pattern in PARSERS:
                if str(input_path).endswith(pattern):
                    print(f"\tFile {input_path} is {pattern}")
                    structure_file_id = input_path.name[:-len(pattern)]
                    structure_files_paths[structure_file_id] = input_path
                    continue
        elif input_path.is_dir():
            files_inside_input_path = list(input_path.glob("**/*"))
            for pattern in PARSERS:
                structure_file_paths = list(
                    filter(lambda x: str(x).endswith(pattern),
                           files_inside_input_path))
                if len(structure_file_paths) == 0:
                    continue
                print(
                    f"\tFound {len(structure_file_paths)} {pattern} files in {input_path}"
                )
                structure_file_ids = [
                    x.name[:-len(pattern)] for x in structure_file_paths
                ]
                structure_files_paths.update(
                    zip(structure_file_ids, structure_file_paths))
        else:
            print(f"\tUnable to find {str(input_path)}")
    return structure_files_paths


def read_structure_file(file_path: pathlib.Path) -> SeqAtoms:
    """
    Extract sequence and atom positions from structure file
    :param file_path:
    :return:
    """
    for pattern, structure_parsing_func in PARSERS.items():
        if file_path.name.endswith(pattern):
            if file_path.name.endswith('.gz'):
                f = gzip.open(file_path, 'rt')
            else:
                f = open(file_path, 'r')
            atom_amino_group, positions, groups = structure_parsing_func(f)
            f.close()

            protein_id = file_path.name.replace(pattern, '')

    return SeqAtoms(protein_id, atom_amino_group, positions, groups)


def save_sequence_and_atoms(seq_atoms: SeqAtoms, sequence_path: pathlib.Path,
                            atoms_path: pathlib.Path,
                            max_target_chain_length: int) -> str:
    """
    Reads the structure file extracting sequence and atom positions.
    If the sequence is too short, it is skipped as it is considered to be corrupted or DNA sequence.
    If the sequence is too long, it is truncated based on the max_target_chain_length defined in target_db_config.json.

    Saves the extracted data:
        sequence - protein_id.faa file containing single sequence f.write(f">{protein_id}\n{sequence}\n")
            SEQ_ATOMS_DATASET_PATH / SEQUENCES / (protein_id + ".faa")
        atom positions - binary file containing positions of all atoms and correlated amino acid chain index
            SEQ_ATOMS_DATASET_PATH / ATOMS / (protein_id + ".bin")
            For more information on how those binary are saved check out source code at CPP_lib/atoms_file_io.h

    Args:
        seq_atoms: An object of SeqAtoms containing the sequence and atom positions.
        sequence_path: A pathlib.Path object containing the path to the sequence file to be saved.
        atoms_path: A pathlib.Path object containing the path to the atom positions file to be saved.
        max_target_chain_length: An integer value of the maximum length of the protein chain.

    Returns:
        If the protein chain is too short, it returns "Sequence is too short, probably DNA or corrupted".
        If the protein chain is truncated, it returns the string with the new length of the chain.
        Otherwise, it returns "SUCCEED".
    """
    _, groups_index = np.unique(seq_atoms.groups, return_index=True)
    groups_index.sort()

    if len(groups_index) < 9:
        # Files containing protein chains shorter than 9 amino acids might be corrupted or contain DNA
        return "FAIL - Sequence is too short, probably DNA or corrupted"

    truncated = False
    if len(groups_index) > max_target_chain_length:
        # Protein chains with more than {max_target_chain_length} are truncated.
        # Change CONFIG/RUNTIME_PARAMETERS.py :max_target_chain_length
        # and rerun this script with --overwrite to process this file again
        truncated = True
        group_indexes = groups_index[:max_target_chain_length]
        group_indexes = np.append(
            group_indexes,
            groups_index[max_target_chain_length]).astype(np.int32)
    else:
        group_indexes = np.append(
            groups_index, seq_atoms.positions.shape[0]).astype(np.int32)

    sequence = ''.join([
        bio_utils.PROTEIN_LETTERS[seq_atoms.atom_amino_group[i]]
        for i in group_indexes[:-1]
    ])
    with open(sequence_path, "w", encoding="utf-8") as f:
        f.write(f">{seq_atoms.protein_id}\n{sequence}\n")

    libAtomDistanceIO.save_atoms(seq_atoms.positions, group_indexes,
                                 str(atoms_path))

    if truncated:
        return f"SUCCESS, but sequences and contact maps got truncated to {max_target_chain_length}"
    else:
        return "SUCCESS"


def process_structure_file(structure_file, save_path, max_target_chain_length):
    """

    :param structure_file:
    :param save_path:
    :param max_target_chain_length:
    :return:
    """
    try:
        seq_atoms = read_structure_file(structure_file)
    except Exception:
        print("EXCEPTION WHILE READING FILE ", str(structure_file))
        logging.error(traceback.format_exc())
        return "file reading exceptions"

    # create directories for saving if they don't exist
    sequences_folder = save_path / SEQUENCES
    atoms_path = save_path / ATOMS
    sequences_folder.mkdir(parents=True, exist_ok=True)
    atoms_path.mkdir(parents=True, exist_ok=True)

    # process single protein
    prot_sequence_path = save_path / SEQUENCES / (seq_atoms.protein_id +
                                                  ".faa")
    prot_atoms_path = save_path / ATOMS / (seq_atoms.protein_id + ".bin")

    try:
        return save_sequence_and_atoms(seq_atoms, prot_sequence_path,
                                       prot_atoms_path,
                                       max_target_chain_length)
    except Exception:
        print("EXCEPTION DURING FILE PROCESSING ", str(structure_file))
        logging.error(traceback.format_exc())
        prot_sequence_path.unlink(missing_ok=True)
        prot_atoms_path.unlink(missing_ok=True)
        return "file processing exceptions"
