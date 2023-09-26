# Create logger
import logging
from pathlib import Path

from mDeepFRI import MERGED_SEQUENCES, TARGET_MMSEQS_DB_NAME
from mDeepFRI.mmseqs import (check_mmseqs_database, create_target_database,
                             extract_fasta_foldcomp)

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def build_database(
    input_path: str,
    output_path: str,
    overwrite: bool,
    threads: int,
) -> None:
    """
    Extracts FASTA file from FoldComp database. Creates MMSeqs2 database and index.

    Args:
        input_path (str): path to a FoldComp database with compressed structures.
        output_path (str): path to folder where the database for DeepFRI will be created.
        overwrite (bool): overwrite existing database.
        threads (int): number of threads to use.

    Returns:
        None
    """

    logging.info("Building MMSeqs2 database from %s", input_path)
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    output_sequences = output_path / MERGED_SEQUENCES
    # check if files exist in output directory
    if output_sequences.exists() and not overwrite:
        logging.info("Found %s in %s", MERGED_SEQUENCES, output_path)
        logging.info(
            "Skipping extraction of FASTA file from FoldComp database.")
    else:
        logging.info("Extracting FASTA file from FoldComp database.")
        extract_fasta_foldcomp(input_path, output_sequences, threads)
        logging.info("FASTA file extracted to %s", output_sequences)

    # create mmseqs db
    target_db = check_mmseqs_database(output_path / TARGET_MMSEQS_DB_NAME)
    if not target_db and not overwrite:
        logging.info("Found %s in %s", TARGET_MMSEQS_DB_NAME, output_path)
        logging.info("Skipping creation of MMSeqs2 database.")
    else:
        logging.info("Creating MMSeqs2 database.")
        create_target_database(output_sequences,
                               output_path / TARGET_MMSEQS_DB_NAME)
        logging.info("Database created at %s", output_path)

    return (output_sequences, output_path / TARGET_MMSEQS_DB_NAME)
