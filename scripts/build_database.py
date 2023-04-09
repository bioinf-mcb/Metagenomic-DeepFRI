import argparse
import multiprocessing
import pathlib

from itertools import repeat

from typing import List
# Create logger
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)

import numpy as np

from meta_deepFRI.CPP_lib import libAtomDistanceIO
from meta_deepFRI.config.names import ATOMS

from meta_deepFRI.structure_files.parse_structure_file import process_structure_file, search_structure_files
from meta_deepFRI.utils.mmseqs import create_target_database
from meta_deepFRI.utils.utils import shutdown

###################################################################################################################
# utils.mmseqs.update_target_mmseqs_database:

#
# structure_files.parse_structure_file.process_structure_file:
#   1. reads the structure file extracting sequence and atom positions
#   2. skip short sequences and truncate ones that size is over max_protein_length inside target_db_config.json
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

## TODO update documentation
## TODO write tests
## TODO speed up globbing
## TODO replace merging of aminoacid sequence files with pysam.libcfaidx


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(description="Read structure files from folders --input to extract sequence and atom positions. "
                                                 "Create and index new --output MMSEQS2 database")

    # logic described here is implemented in parse_input_paths
    parser.add_argument("-i", "--input", nargs='+', required=True, default=None,
                        help="List of folder or file paths containing structure files. Both absolute and relative to FolderStructureConfig.STRUCTURE_FILES_PATH are accepted."
                             "If not provided pipeline will search in FolderStructureConfig.STRUCTURE_FILES_PATH/--name. "
                             "Use '--input .' to process all files within FolderStructureConfig.STRUCTURE_FILES_PATH")

    parser.add_argument("-o", "--output", required=True, default=None,
                        help="Path to folder where new MMSeqs2 database will be created.")

    parser.add_argument("-t", "--threads", required=False, default=1, type=int,
                        help="Number of threads to use. If not provided, the program will use single core.")
    parser.add_argument("-max_len", "--max_protein_length", required=False, default=1_000, type=int,
                        help="If protein chain is longer than this value, it will be truncated")

    parser.add_argument("--overwrite", action="store_true",
                        help="Flag to override existing sequences and atom positions")
    return parser.parse_args()
    # yapf: enable


def build_database(
    input_paths: List[str],
    output_path: str,
    overwrite: bool,
    threads: int,
    max_protein_length: int = 1_000,
) -> None:
    """
    Searches for structure files in input (both files & folders are accepted, folders are globbed),
    parses them, creates a database of available structures & sequences.
    From sequences creates and indexes MMSeqs 2 database.
    Additionally, creates a .json of added structures.

    Args:
        input_paths (list): list of paths to structure files or folders containing structure files.
        output_path (str): path to folder where the database will be created.
        overwrite (bool): flag to override existing sequences and atom positions.
        threads (int): number of threads to use.
        max_protein_length (int): if protein is longer than this value, it will be truncated. Default is 1_000.

    Returns:
        None
    """
    output_path = pathlib.Path(output_path)
    seq_atoms_path = output_path / "seq_atom_db"
    if not seq_atoms_path.exists():
        seq_atoms_path.mkdir(parents=True)
    target_db_path = output_path / "mmseqs_db"

    # report max_protein_length
    logger.info("MAX_PROTEIN_LENGTH: %s" % max_protein_length)

    logger.info("Searching for structure files in following paths:")
    for input_path in input_paths:
        logger.info("%s" % str(input_path))
    # dict [ structure_file_id ] = pathlib.Path( structure_file_path )
    structure_files_paths = search_structure_files(input_paths)

    n_structures = len(structure_files_paths)
    if n_structures == 0:
        shutdown("No structure files found")
    else:
        logger.info("Found %s structure files to extract" % n_structures)

    # search for already processed protein_ids to skip them
    if not overwrite:
        duplicated_ids_counter = 0
        # remove duplicates from structure_files_paths
        for structure_id in structure_files_paths.keys():
            if (seq_atoms_path / structure_id).exists():
                logger.info(f"Skipping {structure_id} - already in database")
                del structure_files_paths[structure_id]
                duplicated_ids_counter += 1
        logger.info("Found %s duplicated IDs" % duplicated_ids_counter)

    # actual processing of structure files takes place here
    logger.info("Processing %s files" % len(structure_files_paths))
    libAtomDistanceIO.initialize()
    with multiprocessing.Pool(processes=threads) as p:
        processing_status = p.starmap(
            process_structure_file,
            zip(structure_files_paths.values(), repeat(seq_atoms_path.absolute()), repeat(max_protein_length)))

    for struct, status in zip(structure_files_paths.values(), processing_status):
        logging.info(f"{struct} - {status}")

    freshly_added_ids = []
    structure_file_ids = list(structure_files_paths.keys())
    for status, structure_id in zip(processing_status, structure_file_ids):
        if status.startswith("SUCCEED"):
            freshly_added_ids.append(structure_id)

    if len(freshly_added_ids) == 0:
        message = "\nNo new protein structures added.\n No new target database will be created."
        shutdown(message)

    # create a folder for target db
    if not target_db_path.exists():
        target_db_path.mkdir(parents=True)

    create_target_database(seq_atoms_path, target_db_path, freshly_added_ids)


def main() -> None:
    args = parse_args()
    input_seqs = [pathlib.Path(seqs) for seqs in args.input]
    output_path = pathlib.Path(args.output)

    build_database(
        input_paths=input_seqs,
        output_path=output_path,
        overwrite=args.overwrite,
        threads=args.threads,
        max_protein_length=args.max_protein_length)


if __name__ == '__main__':
    main()
