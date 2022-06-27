import dataclasses
import gzip
import logging
import pathlib
import traceback

import numpy as np

from meta_deepFRI import CPP_lib
from meta_deepFRI import structure_files
from CONFIG.FOLDER_STRUCTURE import SEQUENCES, ATOMS

from meta_deepFRI.utils import bio_utils

# need to parse different type of files? Add a pattern with a parser in this dict.
# read structure_files_parsers/README.md for more information about how to create new parser.
PARSERS = {
    '.pdb': structure_files.parse_pdb,
    '.pdb.gz': structure_files.parse_pdb,
    '.cif': structure_files.parse_mmcif,
    '.cif.gz': structure_files.parse_mmcif,
    '.ent': structure_files.parse_pdb,
    '.ent.gz': structure_files.parse_pdb
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
            for pattern in structure_files.PARSERS.keys():
                if str(input_path).endswith(pattern):
                    print(f"\tFile {input_path} is {pattern}")
                    structure_file_id = input_path.name[:-len(pattern)]
                    structure_files_paths[structure_file_id] = input_path
                    continue
        elif input_path.is_dir():
            files_inside_input_path = list(input_path.glob("**/*"))
            for pattern in structure_files.PARSERS.keys():
                structure_file_paths = list(filter(lambda x: str(x).endswith(pattern), files_inside_input_path))
                if len(structure_file_paths) == 0:
                    continue
                print(f"\tFound {len(structure_file_paths)} {pattern} files in {input_path}")
                structure_file_ids = [x.name[:-len(pattern)] for x in structure_file_paths]
                structure_files_paths.update(zip(structure_file_ids, structure_file_paths))
        else:
            print(f"\tUnable to find {str(input_path)}")
    return structure_files_paths


def read_structure_file(file_path: pathlib.Path) -> SeqAtoms:
    """
    Extract sequence and atom positions from structure file
    :param file_path:
    :return:
    """
    for pattern in PARSERS.keys():
        if file_path.name.endswith(pattern):
            if file_path.name.endswith('.gz'):
                f = gzip.open(file_path, 'rt')
            else:
                f = open(file_path, 'r')
            atom_amino_group, positions, groups = PARSERS[pattern](f)
            f.close()

            protein_id = file_path.name.replace(pattern, '')
            return SeqAtoms(protein_id, atom_amino_group, positions, groups)


def save_sequence_and_atoms(seq_atoms: SeqAtoms, sequence_path: pathlib.Path, atoms_path: pathlib.Path,
                            max_target_chain_length: int) -> str:
    """

    :param seq_atoms:
    :param sequence_path:
    :param atoms_path:
    :param max_target_chain_length:
    :return:
    """
    _, groups = np.unique(seq_atoms.groups, return_index=True)
    groups.sort()

    if len(groups) < 9:
        # Files containing protein chains shorter than 9 amino acids might be corrupted or contain DNA
        return "sequences too short, probably DNA or corrupted"

    truncated = False
    if len(groups) > max_target_chain_length:
        # Protein chains with more than {max_target_chain_length} are truncated.
        # Change CONFIG/RUNTIME_PARAMETERS.py :max_target_chain_length
        # and rerun this script with --overwrite to process this file again
        truncated = True
        group_indexes = groups[:max_target_chain_length]
        group_indexes = np.append(group_indexes, groups[max_target_chain_length]).astype(np.int32)
    else:
        group_indexes = np.append(groups, seq_atoms.positions.shape[0]).astype(np.int32)

    sequence = ''.join([bio_utils.PROTEIN_LETTERS[seq_atoms.atom_amino_group[i]] for i in group_indexes[:-1]])
    with open(sequence_path, "w") as f:
        f.write(f">{seq_atoms.protein_id}\n{sequence}\n")
    CPP_lib.save_atoms(seq_atoms.positions, group_indexes, str(atoms_path))

    if truncated:
        return f"SUCCEED, but sequences and contact maps got truncated to {max_target_chain_length}"
    else:
        return "SUCCEED"


def process_structure_file(structure_file, save_path, max_target_chain_length):
    """

    :param structure_file:
    :param save_path:
    :param max_target_chain_length:
    :return:
    """
    try:
        seq_atoms = structure_files.read_structure_file(structure_file)
    except Exception:
        print("EXCEPTION WHILE READING FILE ", str(structure_file))
        logging.error(traceback.format_exc())
        return "file reading exceptions"

    # process and store sequence and atom positions in SEQ_ATOMS_DATASET_PATH / output_name
    sequence_path = save_path / SEQUENCES / (seq_atoms.protein_id + ".faa")
    atoms_path = save_path / ATOMS / (seq_atoms.protein_id + ".bin")

    try:
        return structure_files.save_sequence_and_atoms(seq_atoms, sequence_path, atoms_path, max_target_chain_length)
    except Exception:
        print("EXCEPTION DURING FILE PROCESSING ", str(structure_file))
        logging.error(traceback.format_exc())
        sequence_path.unlink(missing_ok=True)
        atoms_path.unlink(missing_ok=True)
        return "file processing exceptions"
