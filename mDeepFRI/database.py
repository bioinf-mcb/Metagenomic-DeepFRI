# Create logger
import logging
from pathlib import Path

from mDeepFRI.mmseqs import (create_target_database, extract_fasta_foldcomp,
                             validate_mmseqs_database)

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def build_database(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    threads: int = 1,
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
    output_sequences = output_path / Path(input_path.stem + ".fasta.gz")
    unzipped_sequences = output_sequences.with_suffix("")
    needs_new_mmseqs = False

    # check if files exist in output directory
    if output_sequences.exists() and not overwrite:
        logging.info("Found %s in %s", output_sequences, output_path)
        logging.info(
            "Skipping extraction of FASTA file from FoldComp database.")
    else:
        logging.info("Extracting FASTA file from FoldComp database.")
        output_sequences = extract_fasta_foldcomp(input_path,
                                                  unzipped_sequences, threads)
        logging.info("FASTA file extracted to %s", output_sequences)
        needs_new_mmseqs = True

    # create mmseqs db
    mmseqs_path = output_path / Path(input_path.stem + ".mmseqsDB")
    mmseqs_valid = validate_mmseqs_database(mmseqs_path)
    if not mmseqs_valid:
        logging.info("Creating MMSeqs2 database.")
        create_target_database(output_sequences, mmseqs_path)
    elif overwrite or needs_new_mmseqs:
        logging.info("Creating MMSeqs2 database.")
        create_target_database(output_sequences, mmseqs_path)
    else:
        logging.info("Database created at %s", output_path)
        logging.info("Found %s in %s", mmseqs_path, output_path)
        logging.info("Skipping creation of MMSeqs2 database.")

    return (output_sequences, mmseqs_path)
