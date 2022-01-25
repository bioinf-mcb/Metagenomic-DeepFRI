import argparse
import gzip
import json
import logging
import multiprocessing
import pathlib
import traceback

from itertools import repeat
import numpy as np

from CONFIG.FOLDER_STRUCTURE import TARGET_MMSEQS_DB_NAME, ATOMS, SEQUENCES, STRUCTURE_FILES_PATH, \
    DEFAULT_NAME, SEQ_ATOMS_DATASET_PATH, MMSEQS_DATABASES_PATH, MERGED_SEQUENCES
from CONFIG.RUNTIME_PARAMETERS import CPU_COUNT, MAX_TARGET_CHAIN_LENGTH
from utils.bio_utils import PROTEIN_LETTERS

from structure_files_parsers.parse_pdb import parse_pdb
from structure_files_parsers.parse_mmcif import parse_mmcif

from CPP_lib.libAtomDistanceIO import save_atoms
from CPP_lib.libAtomDistanceIO import initialize as initialize_CPP_LIB
from utils.mmseqs_utils import mmseqs_createdb
from utils.mmseqs_utils import mmseqs_createindex
from utils.utils import create_unix_timestamp_folder, merge_files_binary

STRUCTURE_FILES_PARSERS = {
    '.pdb': parse_pdb,
    '.pdb.gz': parse_pdb,
    '.cif': parse_mmcif,
    '.cif.gz': parse_mmcif,
    '.ent': parse_pdb,
    '.ent.gz': parse_pdb
}


def parse_args():
    parser = argparse.ArgumentParser(description="Read structure files from folders --input to extract sequence and atom positions. "
                                                 "Create and index new --output MMSEQS2 database")

    parser.add_argument("-t", "--target_db_name", required=False, default=DEFAULT_NAME, help="Name of the target database")
    parser.add_argument("-i", "--input", nargs='+', required=False, default=None,
                        help=f"List of folder or file paths containing structure files. Both absolute and relative to {STRUCTURE_FILES_PATH} are accepted."
                             f"If not provided pipeline will search in {STRUCTURE_FILES_PATH}/--name. "
                             f"Use '--input all' to process all files within {STRUCTURE_FILES_PATH}")  # logic described here is implemented at the bottom of this file

    parser.add_argument("--overwrite", action="store_true", help="Flag to override existing sequences and atom positions")
    return parser.parse_args()


def parse_structure_file(structure_file, save_path):
    protein_id = structure_file.name

    # extract sequence and atom positions from structure file
    try:
        for pattern in STRUCTURE_FILES_PARSERS.keys():
            if protein_id.endswith(pattern):
                if protein_id.endswith('.gz'):
                    with gzip.open(structure_file, 'rt') as f:
                        atom_amino_group, positions, groups = STRUCTURE_FILES_PARSERS[pattern](f)
                else:
                    with open(structure_file, 'r') as f:
                        atom_amino_group, positions, groups = STRUCTURE_FILES_PARSERS[pattern](f)
                protein_id = structure_file.name.replace(pattern, '')
                break

    except Exception:
        print("EXCEPTION WHILE READING FILE ", str(structure_file))
        logging.error(traceback.format_exc())
        return "file reading exceptions"

    # process and store sequence and atom positions in SEQ_ATOMS_DATASET_PATH / output_name
    sequence_path = save_path / SEQUENCES / (protein_id + ".faa")
    atoms_path = save_path / ATOMS / (protein_id + ".bin")
    try:
        _, groups = np.unique(groups, return_index=True)
        groups.sort()

        if len(groups) < 9:
            # print("Files containing protein chains shorter than 9 amino acids might be corrupted or contain DNA ", file)
            return "sequences too short, probably DNA or corrupted"

        truncated = False
        if len(groups) > MAX_TARGET_CHAIN_LENGTH:
            # print(f"Protein chains with more than {MAX_TARGET_CHAIN_LENGTH} are truncated. {str(structure_file)}\n"
            #       f"Change CONFIG/RUNTIME_PARAMETERS.py :MAX_TARGET_CHAIN_LENGTH and rerun this script with --overwrite to process this file again")
            truncated = True
            group_indexes = groups[:MAX_TARGET_CHAIN_LENGTH]
            group_indexes = np.append(group_indexes, groups[MAX_TARGET_CHAIN_LENGTH]).astype(np.int32)
        else:
            group_indexes = np.append(groups, positions.shape[0]).astype(np.int32)

        sequence = ''.join([PROTEIN_LETTERS[atom_amino_group[i]] for i in group_indexes[:-1]])
        with open(sequence_path, "w") as f:
            f.write(f">{protein_id}\n{sequence}\n")
        save_atoms(positions, group_indexes, str(atoms_path))

        if truncated:
            return f"SUCCEED, but sequences and contact maps got truncated to {MAX_TARGET_CHAIN_LENGTH}"
        else:
            return "SUCCEED"

    except Exception:
        print("EXCEPTION DURING FILE PROCESSING ", str(structure_file))
        logging.error(traceback.format_exc())
        sequence_path.unlink(missing_ok=True)
        atoms_path.unlink(missing_ok=True)
        return "file processing exceptions"


def main(input_paths, target_db_name, overwrite):
    seq_atoms_path = SEQ_ATOMS_DATASET_PATH / target_db_name
    seq_atoms_path.mkdir(exist_ok=True, parents=True)
    print("Sequences and Atoms positions will be stored in: ", seq_atoms_path)
    (seq_atoms_path / SEQUENCES).mkdir(exist_ok=True)
    (seq_atoms_path / ATOMS).mkdir(exist_ok=True)

    print("Searching for structure files in paths:")
    for input_path in input_paths:
        print(f"\t{str(input_path)}")
    # dict [ structure_file_id ] = pathlib.Path( structure_file_path )
    structure_files = dict()

    for input_path in input_paths:
        print(f"{str(input_path)}")
        if input_path.is_file():
            for pattern in STRUCTURE_FILES_PARSERS.keys():
                if str(input_path).endswith(pattern):
                    print(f"\tFile {input_path} is {pattern}")
                    structure_file_id = input_path.name[:-len(pattern)]
                    structure_files[structure_file_id] = input_path
                    continue
        elif input_path.is_dir():
            files_inside_input_path = list(input_path.glob("**/*"))
            for pattern in STRUCTURE_FILES_PARSERS.keys():
                structure_file_paths = list(filter(lambda x: str(x).endswith(pattern), files_inside_input_path))
                structure_file_ids = [x.name[:-len(pattern)] for x in structure_file_paths]
                if len(structure_file_paths) == 0:
                    continue
                print(f"\tFound {len(structure_file_paths)} {pattern} files in {input_path}")
                structure_files.update(zip(structure_file_ids, structure_file_paths))
        else:
            print(f"\tUnable to find {str(input_path)}")

    if len(structure_files) == 0:
        print("No structure files found")
        return
    print(f"Found {len(structure_files)} structure files to extract")

    # search for already processed protein_ids to skip them
    if not overwrite:
        duplicated_ids_counter = 0
        for structure_id in list(structure_files.keys()):
            if (seq_atoms_path / ATOMS / (structure_id + ".bin")).exists():
                structure_files.pop(structure_id)
                duplicated_ids_counter += 1
        print(f"Found {duplicated_ids_counter} duplicated IDs")

    print("\nProcessing", len(structure_files), "files")
    initialize_CPP_LIB()
    with multiprocessing.Pool(processes=CPU_COUNT) as p:
        results = p.starmap(parse_structure_file, zip(structure_files.values(), repeat(seq_atoms_path.absolute())))

    result_types, types_counts = np.unique(results, return_counts=True)
    for i in range(len(result_types)):
        print(f"\t{types_counts[i]} {result_types[i]}")

    freshly_added_ids = []
    structure_file_ids = list(structure_files.keys())
    for i in range(len(results)):
        if results[i].startswith("SUCCEED"):
            freshly_added_ids.append(structure_file_ids[i])

    if len(freshly_added_ids) == 0:
        print(f"\nNo new protein structures added.\nNo new target database will be created.")
        return
    freshly_added_ids = sorted(freshly_added_ids)

    sequence_files = list((seq_atoms_path / SEQUENCES).glob("**/*.faa"))
    print("\nMerging " + str(len(sequence_files)) + " sequence files for mmseqs2")
    merge_files_binary(sequence_files, seq_atoms_path / MERGED_SEQUENCES)

    mmseqs2_path = create_unix_timestamp_folder(MMSEQS_DATABASES_PATH / target_db_name)
    print("Creating new target mmseqs2 database " + str(mmseqs2_path))
    mmseqs_createdb(seq_atoms_path / MERGED_SEQUENCES, mmseqs2_path / TARGET_MMSEQS_DB_NAME)
    print("Indexing new target mmseqs2 database " + str(mmseqs2_path))
    mmseqs_createindex(mmseqs2_path / TARGET_MMSEQS_DB_NAME)
    print("Saving freshly added sequence ids to " + str(mmseqs2_path / "new_structure_ids.json"))
    json.dump(freshly_added_ids, open(mmseqs2_path / "new_structure_ids.json", "w"), indent=4)


if __name__ == '__main__':
    args = parse_args()

    target_db_name = pathlib.Path(args.target_db_name)
    overwrite = args.overwrite

    # here is the logic responsible for selecting correct structure files folders and files based on args
    if args.input is None:
        input_paths = [pathlib.Path(STRUCTURE_FILES_PATH / target_db_name)]
    else:
        if args.input == ["all"]:
            input_paths = [STRUCTURE_FILES_PATH]
        else:
            input_paths = []
            for input_path in [pathlib.Path(x) for x in args.input]:
                if input_path.is_absolute():
                    input_paths.append(input_path)
                else:
                    input_paths.append(STRUCTURE_FILES_PATH / input_path)

    main(input_paths, target_db_name, overwrite)
