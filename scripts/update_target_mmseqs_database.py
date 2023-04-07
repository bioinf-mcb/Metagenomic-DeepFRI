import argparse
import json
import multiprocessing
import pathlib

from itertools import repeat
import numpy as np

from meta_deepFRI.config.folder_structure import FolderStructureConfig, load_folder_structure_config
from meta_deepFRI import CPP_lib
from meta_deepFRI.config.names import DEFAULT_NAME, ATOMS, TARGET_DB_CONFIG
from meta_deepFRI.config import CPU_COUNT

from meta_deepFRI.structure_files.parse_structure_file import process_structure_file, search_structure_files
from meta_deepFRI.utils.mmseqs import create_target_database
from meta_deepFRI.utils import create_unix_timestamp_folder, parse_input_paths

###################################################################################################################
# utils.mmseqs.update_target_mmseqs_database:

#
# structure_files.parse_structure_file.process_structure_file:
#   1. reads the structure file extracting sequence and atom positions
#   2. skip short sequences and truncate ones that size is over MAX_TARGET_CHAIN_LENGTH inside target_db_config.json
#   3. save extracted data:
#           sequence - protein_id.faa file containing single sequence f.write(f">{protein_id}\n{sequence}\n")
#               SEQ_ATOMS_DATASET_PATH / project_name / SEQUENCES / (protein_id + ".faa")
#           atom positions - binary file containing positions of all atoms and correlated amino acid chain index
#               SEQ_ATOMS_DATASET_PATH / project_name / ATOMS / (protein_id + ".bin")
#               For more information on how those binary are saved check out source code at CPP_lib/atoms_file_io.h
#
# create_target_database:
#   1.  merges all sequences from SEQ_ATOMS_DATASET_PATH / project_name / SEQUENCES
#   2.  creates and index a new mmseqs target database inside MMSEQS_DATABASES_PATH / project_name / timestamp
#   3.  creates a structure_ids_added.json, so it is easier to remove unwanted ATOM and SEQUENCE after mistake
###################################################################################################################


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(description="Read structure files from folders --input to extract sequence and atom positions. "
                                                 "Create and index new --output MMSEQS2 database")

    parser.add_argument("-r", "--data_root", required=False, default="default_config/data_root_path.json",
                        help="Path to json file containing DATA_ROOT or path to folder that will be used as DATA_ROOT")

    parser.add_argument("-p", "--project_name", required=False, default=DEFAULT_NAME,
                        help="Name for the target database")

    # logic described here is implemented in parse_input_paths
    parser.add_argument("-i", "--input", nargs='+', required=False, default=None,
                        help=f"List of folder or file paths containing structure files. Both absolute and relative to FolderStructureConfig.STRUCTURE_FILES_PATH are accepted."
                             f"If not provided pipeline will search in FolderStructureConfig.STRUCTURE_FILES_PATH/--name. "
                             f"Use '--input .' to process all files within FolderStructureConfig.STRUCTURE_FILES_PATH")

    # todo add max_target_chain_length argument and inform user if there is difference between this arg and existing target_db_config.json
    # parser.add_argument("-m", "--max_target_chain_length", required=False, default=None,
    #                     help="If protein chain is longer than this value, it will be truncated")
    parser.add_argument("--overwrite", action="store_true",
                        help="Flag to override existing sequences and atom positions")
    return parser.parse_args()
    # yapf: enable


def update_target_mmseqs_database(fsc: FolderStructureConfig, input_paths, project_name, overwrite) -> None:
    """
    1.  iterates over --input searching for file extensions that match the PARSERS keys.
    2.  filter out protein_ids that already exists in SEQ_ATOMS_DATASET_PATH / project_name / ATOMS
    3.  process_structure_file in parallel
    :param fsc:
    :param input_paths:
    :param project_name:
    :param overwrite:
    :return:
    """
    seq_atoms_path = fsc.SEQ_ATOMS_DATASET_PATH / project_name
    if not seq_atoms_path.exists():
        raise RuntimeError(f"No project found. Please run: python create_project --project_name {project_name} first")
    print("Sequences and Atoms positions will be stored in: ", seq_atoms_path)

    # get MAX_TARGET_CHAIN_LENGTH from existing config or load from project config
    target_db_path = fsc.MMSEQS_DATABASES_PATH / project_name
    print(f"Using TARGET_DB_CONFIG {target_db_path / TARGET_DB_CONFIG}")
    max_target_chain_length = json.load(open(target_db_path / TARGET_DB_CONFIG))
    print(f"\t{max_target_chain_length} - MAX_TARGET_CHAIN_LENGTH")

    print("Searching for structure files in paths:")
    for input_path in input_paths:
        print(f"\t{str(input_path)}")
    # dict [ structure_file_id ] = pathlib.Path( structure_file_path )
    structure_files_paths = search_structure_files(input_paths)
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
        print("\n No new protein structures added.\n No new target database will be created.")
        return

    new_mmseqs2_db_path = create_unix_timestamp_folder(fsc.MMSEQS_DATABASES_PATH / project_name)
    create_target_database(seq_atoms_path, new_mmseqs2_db_path, freshly_added_ids)


def main():
    args = parse_args()

    data_root = pathlib.Path(args.data_root)
    if data_root.is_dir():
        fsc = FolderStructureConfig(data_root)
    else:
        fsc = load_folder_structure_config(data_root)

    project_name = pathlib.Path(args.project_name)
    overwrite = args.overwrite
    input_paths = parse_input_paths(args.input, project_name, fsc.STRUCTURE_FILES_PATH)

    update_target_mmseqs_database(fsc, input_paths, project_name, overwrite)


if __name__ == '__main__':
    main()
