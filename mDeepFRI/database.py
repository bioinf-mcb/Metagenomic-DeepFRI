import json
# Create logger
import logging
import multiprocessing
import pathlib
from itertools import repeat
from typing import List

import numpy as np

from mDeepFRI.config.names import SEQ_ATOMS_DATASET_PATH
from mDeepFRI.CPP_lib import libAtomDistanceIO  # type: ignore[attr-defined]
from mDeepFRI.structure_files.parse_structure_file import (
    process_structure_file, search_structure_files)
from mDeepFRI.utils.mmseqs import create_target_database
from mDeepFRI.utils.utils import shutdown

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


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

    seq_atoms_path = output_path / SEQ_ATOMS_DATASET_PATH
    if not seq_atoms_path.exists():
        seq_atoms_path.mkdir(parents=True)

    # report max_protein_length
    logger.info("MAX_PROTEIN_LENGTH: %s", max_protein_length)

    logger.info("Searching for structure files in following paths:")
    for input_path in input_paths:
        logger.info("%s", str(input_path))
    # dict [ structure_file_id ] = pathlib.Path( structure_file_path )
    structure_files_paths = search_structure_files(input_paths)

    n_structures = len(structure_files_paths)
    if n_structures == 0:
        shutdown("No structure files found")
    else:
        logger.info("Found %s structure files to extract", n_structures)

    # search for already processed protein_ids to skip them
    if not overwrite:
        duplicated_ids_counter = 0
        # remove duplicates from structure_files_paths
        for structure_id in structure_files_paths:
            if (seq_atoms_path / structure_id).exists():
                logger.info("Skipping %s - already in database", structure_id)
                del structure_files_paths[structure_id]
                duplicated_ids_counter += 1
        logger.info("Found %s duplicated IDs", duplicated_ids_counter)

    # actual processing of structure files takes place here
    logger.info("Processing %s files", len(structure_files_paths))
    libAtomDistanceIO.initialize()
    with multiprocessing.Pool(processes=threads) as p:
        processing_status = p.starmap(
            process_structure_file,
            zip(structure_files_paths.values(),
                repeat(seq_atoms_path.absolute()), repeat(max_protein_length)))

    # General stats over database
    status, status_count = np.unique(processing_status, return_counts=True)

    # reversing the list, so that FAIL is the last one printed
    for stat, count in list(zip(status, status_count))[::-1]:
        logger.info("%s - %i", stat, count)

    # Print out failed structures
    for struct, status in zip(structure_files_paths.values(),
                              processing_status):
        if status.startswith("FAIL"):
            logger.error("%s - %s", struct, status)

    freshly_added_ids = []
    structure_file_ids = list(structure_files_paths.keys())
    for status, structure_id in zip(processing_status, structure_file_ids):
        if status.startswith("SUCCESS"):
            freshly_added_ids.append(structure_id)

    if len(freshly_added_ids) == 0:
        message = "\nNo new protein structures added.\n No new target database will be created."
        shutdown(message)

    # save db params
    param_dict = {}
    param_dict["sequences"] = sorted(freshly_added_ids)
    param_dict["MAX_PROTEIN_LENGTH"] = max_protein_length
    param_dict["input_structures_path"] = [
        str(path) for path in sorted(input_paths)
    ]

    json.dump(param_dict,
              open(output_path / "db_params.json", "w", encoding="utf-8"),
              indent=4)

    create_target_database(seq_atoms_path, output_path)
