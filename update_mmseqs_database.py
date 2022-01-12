import argparse
import gzip
import logging
import multiprocessing
import os
import pathlib
import traceback

from itertools import repeat
import numpy as np

from CONFIG.FOLDER_STRUCTURE import TARGET_DB_NAME, ATOMS, SEQUENCES, STRUCTURE_FILES_PATH, SEQ_ATOMS_DATASET_PATH, MMSEQS_DATABASES_PATH
from CONFIG.RUNTIME_PARAMETERS import CPU_COUNT, MAX_CHAIN_LENGTH
from utils.bio_utils import PROTEIN_LETTERS, STRUCTURE_FILES_PATTERNS

from CPP_lib.libAtomDistanceIO import save_atoms
from CPP_lib.libAtomDistanceIO import initialize as initialize_CPP_LIB
from utils.structure_files_parsers.parse_mmcif import parse_mmcif
from utils.structure_files_parsers.parse_pdb import parse_pdb
from utils.mmseqs_utils import mmseqs_createdb
from utils.mmseqs_utils import mmseqs_createindex
from utils.utils import create_unix_time_folder


def parse_args():
    parser = argparse.ArgumentParser(description="Read structure files from -i to extract sequence and atom positions. "
                                                 "Save them in -o as .faa and .bin files. "
                                                 "Create and index new MMSEQS2 database in -db")

    parser.add_argument("-i", "--input", required=False, default=STRUCTURE_FILES_PATH, help="Path to folder with structure files")
    parser.add_argument("-o", "--output", required=False, default=SEQ_ATOMS_DATASET_PATH, help="Path to folder with sequences and atom positions")
    parser.add_argument("-db", "--database", required=False, default=MMSEQS_DATABASES_PATH, help="Path to create new MMSEQS2 database")
    parser.add_argument("--overwrite", action="store_true", help="Flag to override existing sequence and atom positions")
    return parser.parse_args()


def parse_structure_file(structure_file, save_path):
    file_id = structure_file.name
    sequence_path = pathlib.Path(save_path) / SEQUENCES / (file_id + ".faa")
    atoms_path = pathlib.Path(save_path) / ATOMS / (file_id + ".bin")

    try:
        if file_id.endswith('.pdb'):
            file_id = file_id.replace('.pdb', '')
            with open(structure_file, 'r') as f:
                atom_amino_group, positions, groups = parse_pdb(f)
        elif file_id.endswith('.pdb.gz'):
            file_id = file_id.replace('.pdb.gz', '')
            with gzip.open(structure_file, 'rt') as f:
                atom_amino_group, positions, groups = parse_pdb(f)
        elif file_id.endswith('.cif'):
            file_id = file_id.replace('.cif', '')
            with open(structure_file, 'r') as f:
                atom_amino_group, positions, groups = parse_mmcif(f)
        elif file_id.endswith('.cif.gz'):
            file_id = file_id.replace('.cif.gz', '')
            with gzip.open(structure_file, 'rt') as f:
                atom_amino_group, positions, groups = parse_mmcif(f)
        else:
            print("Unsupported file format of file " + str(structure_file))
            return

    except Exception:
        print("EXCEPTION WHILE READING FILE ", str(structure_file))
        logging.error(traceback.format_exc())
        return

    try:
        _, groups = np.unique(groups, return_index=True)
        groups.sort()

        if len(groups) < 9:
            # print("Files containing protein chains shorter than 9 amino acids might be corrupted or contain DNA ", file)
            return

        if len(groups) > MAX_CHAIN_LENGTH:
            # print(f"Protein chains with more than {MAX_CHAIN_LENGTH} are truncated. {str(structure_file)}\n"
            #       f"Change CONFIG/RUNTIME_PARAMETERS.py :MAX_CHAIN_LENGTH and rerun this script with --overwrite to apply process this file again")
            group_indexes = groups[:MAX_CHAIN_LENGTH]
            group_indexes = np.append(group_indexes, groups[MAX_CHAIN_LENGTH]).astype(np.int32)
        else:
            group_indexes = np.append(groups, positions.shape[0]).astype(np.int32)

        sequence = ''.join([PROTEIN_LETTERS[atom_amino_group[i]] for i in group_indexes[:-1]])
        with open(sequence_path, "w") as f:
            f.write(">" + file_id + "\n" + sequence + "\n")
        save_atoms(positions, group_indexes, str(atoms_path))

    except Exception:
        print("EXCEPTION DURING FILE PROCESSING ", str(structure_file))
        logging.error(traceback.format_exc())
        sequence_path.unlink(missing_ok=True)
        atoms_path.unlink(missing_ok=True)
        return


def main(input_path, atoms_path, db_path, overwrite):
    atoms_path.mkdir(exist_ok=True, parents=True)
    save_path = str(atoms_path.absolute())
    print("Sequences and Atoms positions will be stored in: ", save_path)
    (atoms_path / SEQUENCES).mkdir(exist_ok=True)
    (atoms_path / ATOMS).mkdir(exist_ok=True)

    print("Searching for structure files in: ", input_path)
    structure_files = dict()
    for pattern in STRUCTURE_FILES_PATTERNS:
        pattern_structure_files = list(input_path.glob("**/*" + pattern))
        pattern_structure_ids = [x.name[:-len(pattern)] for x in pattern_structure_files]
        if len(pattern_structure_files) == 0:
            continue
        print("Found", len(pattern_structure_files), pattern, "files")
        structure_files.update(zip(pattern_structure_ids, pattern_structure_files))

    if len(structure_files) == 0:
        print("No structure files found")
        return

    if not overwrite:
        existing_structures = set([x.name[:-4] for x in (atoms_path / ATOMS).glob("**/*.bin")])
        print("Found ", len(existing_structures), " already processed structures")
        duplicated_ids_counter = 0
        for id in list(structure_files.keys()):
            if id in existing_structures:
                structure_files.pop(id)
                duplicated_ids_counter += 1
        print("Found ", duplicated_ids_counter, " duplicated IDs")

    initialize_CPP_LIB()
    print("Processing", len(structure_files), "files")
    with multiprocessing.Pool(processes=CPU_COUNT) as p:
        p.starmap(parse_structure_file, zip(structure_files.values(), repeat(save_path)))

    print("Merging sequence files for mmseqs2")
    os.system(f"cat {atoms_path / SEQUENCES}/* > {atoms_path / 'merged_sequences.faa'}")

    mmseqs2_path = create_unix_time_folder(db_path)
    print("Creating new mmseqs2 database " + str(mmseqs2_path))
    mmseqs_createdb(atoms_path / 'merged_sequences.faa', mmseqs2_path / TARGET_DB_NAME)
    print("Indexing new mmseqs2 database " + str(mmseqs2_path))
    mmseqs_createindex(mmseqs2_path / TARGET_DB_NAME)


if __name__ == '__main__':
    args = parse_args()

    input_path = pathlib.Path(args.input)
    output_path = pathlib.Path(args.output)
    db_path = pathlib.Path(args.database)
    overwrite = args.overwrite

    main(input_path, output_path, db_path, overwrite)
