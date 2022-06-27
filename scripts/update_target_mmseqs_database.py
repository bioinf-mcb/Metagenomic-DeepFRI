import argparse
import json
import multiprocessing
import pathlib

from itertools import repeat
import numpy as np

from CONFIG.FOLDER_STRUCTURE import ATOMS, SEQUENCES, STRUCTURE_FILES_PATH, \
    DEFAULT_NAME, SEQ_ATOMS_DATASET_PATH, MMSEQS_DATABASES_PATH, TARGET_DB_CONFIG
from CONFIG.RUNTIME_PARAMETERS import CPU_COUNT
from CONFIG.get_config_dict import target_db_config

from meta_deepFRI import CPP_lib

from meta_deepFRI.structure_files.parse_structure_file import process_structure_file
from meta_deepFRI.utils.mmseqs import create_target_database

###################################################################################################################
# utils.mmseqs.update_target_mmseqs_database:
#   1.  iterates over --input searching for file extensions that match the PARSERS keys.
#   2.  filter out protein_ids that already exists in SEQ_ATOMS_DATASET_PATH / target_db_name / ATOMS
#   3.  process_structure_file in parallel
#
# structure_files.parse_structure_file.process_structure_file:
#   1. reads the structure file extracting sequence and atom positions
#   2. skip short sequences and truncate ones that size is over MAX_TARGET_CHAIN_LENGTH inside target_db_config.json
#   3. save extracted data:
#           sequence - protein_id.faa file containing single sequence f.write(f">{protein_id}\n{sequence}\n")
#               SEQ_ATOMS_DATASET_PATH / target_db_name / SEQUENCES / (protein_id + ".faa")
#           atom positions - binary file containing positions of all atoms and correlated amino acid chain index
#               SEQ_ATOMS_DATASET_PATH / target_db_name / ATOMS / (protein_id + ".bin")
#               For more information on how those binary are saved check out source code at CPP_lib/atoms_file_io.h
#
# create_target_database:
#   1.  merges all sequences from SEQ_ATOMS_DATASET_PATH / target_db_name / SEQUENCES
#   2.  creates and index a new mmseqs target database inside MMSEQS_DATABASES_PATH / target_db_name / timestamp
#   3.  creates a structure_ids_added.json, so it is easier to remove unwanted ATOM and SEQUENCE after mistake
###################################################################################################################


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(description="Read structure files from folders --input to extract sequence and atom positions. "
                                                 "Create and index new --output MMSEQS2 database")

    parser.add_argument("-p", "--project_name", required=False, default=DEFAULT_NAME,
                        help="Name for the target database")

    parser.add_argument("-i", "--input", nargs='+', required=False, default=None,
                        help=f"List of folder or file paths containing structure files. Both absolute and relative to {STRUCTURE_FILES_PATH} are accepted."
                             f"If not provided pipeline will search in {STRUCTURE_FILES_PATH}/--name. "
                             f"Use '--input .' to process all files within {STRUCTURE_FILES_PATH}")  # logic described here is implemented in parse_input_paths

    # todo add max_target_chain_length argument and inform user if there is difference between this arg and existing target_db_config.json
    # parser.add_argument("-m", "--max_target_chain_length", required=False, default=None,
    #                     help="If protein chain is longer than this value, it will be truncated")
    parser.add_argument("--overwrite", action="store_true",
                        help="Flag to override existing sequences and atom positions")
    return parser.parse_args()
    # yapf: enable


def update_target_mmseqs_database(input_paths, project_name, overwrite) -> None:
    """

    :param input_paths:
    :param project_name:
    :param overwrite:
    :return:
    """
    seq_atoms_path = SEQ_ATOMS_DATASET_PATH / project_name
    print("Sequences and Atoms positions will be stored in: ", seq_atoms_path)
    seq_atoms_path.mkdir(exist_ok=True, parents=True)
    (seq_atoms_path / SEQUENCES).mkdir(exist_ok=True)
    (seq_atoms_path / ATOMS).mkdir(exist_ok=True)

    # get MAX_TARGET_CHAIN_LENGTH from existing config or create new config
    target_db_path = MMSEQS_DATABASES_PATH / project_name
    target_db_path.mkdir(exist_ok=True, parents=True)
    if (target_db_path / TARGET_DB_CONFIG).exists():
        print(f"Using existing TARGET_DB_CONFIG {target_db_path / TARGET_DB_CONFIG}")
        max_target_chain_length = json.load(open(target_db_path / TARGET_DB_CONFIG))["MAX_TARGET_CHAIN_LENGTH"]
    else:
        print(f"Creating a new PROJECT_CONFIG {target_db_path / TARGET_DB_CONFIG}")
        config = target_db_config()
        json.dump(config, open(target_db_path / TARGET_DB_CONFIG, "w"), indent=4)
        max_target_chain_length = config["MAX_TARGET_CHAIN_LENGTH"]
    print(f"\t{max_target_chain_length} - MAX_TARGET_CHAIN_LENGTH")

    print("Searching for structure files in paths:")
    for input_path in input_paths:
        print(f"\t{str(input_path)}")
    # dict [ structure_file_id ] = pathlib.Path( structure_file_path )
    structure_files_paths = structure_files.search_structure_files(input_paths)
    if len(structure_files_paths) == 0:
        print("No structure files found")
        return
    print(f"Found {len(structure_files_paths)} structure files to extract")

    # search for already processed protein_ids to skip them
    if not overwrite:
        duplicated_ids_counter = 0
        for structure_id in list(structure_files_paths.keys()):
            if (seq_atoms_path / ATOMS / (structure_id + ".bin")).exists():
                structure_files_paths.pop(structure_id)
                duplicated_ids_counter += 1
        print(f"Found {duplicated_ids_counter} duplicated IDs")

    # actual processing of structure files takes place here
    print("\nProcessing", len(structure_files_paths), "files")
    CPP_lib.initialize()
    with multiprocessing.Pool(processes=CPU_COUNT) as p:
        processing_status = p.starmap(
            process_structure_file,
            zip(structure_files_paths.values(), repeat(seq_atoms_path.absolute()), repeat(max_target_chain_length)))

    status, status_count = np.unique(processing_status, return_counts=True)
    for i in range(len(status)):
        print(f"\t{status_count[i]} - {status[i]}")

    freshly_added_ids = []
    structure_file_ids = list(structure_files_paths.keys())
    for i in range(len(processing_status)):
        if processing_status[i].startswith("SUCCEED"):
            freshly_added_ids.append(structure_file_ids[i])

    if len(freshly_added_ids) == 0:
        print(f"\n No new protein structures added.\n No new target database will be created.")
        return

    create_target_database(seq_atoms_path, project_name, freshly_added_ids)


def main():
    args = parse_args()

    project_name = pathlib.Path(args.project_name)
    overwrite = args.overwrite
    input_paths = utils.parse_input_paths(args.input, project_name, STRUCTURE_FILES_PATH)

    update_target_mmseqs_database(input_paths, project_name, overwrite)


if __name__ == '__main__':
    main()
