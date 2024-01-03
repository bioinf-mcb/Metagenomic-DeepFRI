# Create logger
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mDeepFRI.mmseqs import (create_target_database, extract_fasta_foldcomp,
                             validate_mmseqs_database)

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


@dataclass
class Database:
    """
    Class for storing database paths.
    """
    foldcomp_db: Path
    sequence_db: Path
    mmseqs_db: Path

    def __post_init__(self):
        self.foldcomp_db = Path(self.foldcomp_db)
        self.sequence_db = Path(self.sequence_db)
        self.mmseqs_db = Path(self.mmseqs_db)
        self.name = self.sequence_db.stem.rsplit(".", 1)[0]


@dataclass
class AlignedQuery:
    name: str
    best_hit_name: str
    db_name: str
    identity: float
    aligned_contact_map: np.ndarray


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

    logger.info("Building MMSeqs2 database from %s", input_path)
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    output_sequences = output_path / Path(input_path.stem + ".fasta.gz")
    unzipped_sequences = output_sequences.with_suffix("")
    needs_new_mmseqs = False

    # check if files exist in output directory
    if output_sequences.exists() and not overwrite:
        logger.info("Found %s in %s", output_sequences, output_path)
        logger.info(
            "Skipping extraction of FASTA file from FoldComp database.")
    else:
        logger.info("Extracting FASTA file from FoldComp database.")
        output_sequences = extract_fasta_foldcomp(input_path,
                                                  unzipped_sequences, threads)
        logger.info("FASTA file extracted to %s", output_sequences)
        needs_new_mmseqs = True

    # create mmseqs db
    mmseqs_path = output_path / Path(input_path.stem + ".mmseqsDB")
    mmseqs_valid = validate_mmseqs_database(mmseqs_path)
    if not mmseqs_valid:
        logger.info("Creating MMSeqs2 database.")
        create_target_database(output_sequences, mmseqs_path)
    elif overwrite or needs_new_mmseqs:
        logger.info("Creating MMSeqs2 database.")
        create_target_database(output_sequences, mmseqs_path)
    else:
        logger.info("Database created at %s", output_path)
        logger.info("Found %s in %s", mmseqs_path, output_path)
        logger.info("Skipping creation of MMSeqs2 database.")

    database = Database(foldcomp_db=input_path,
                        sequence_db=output_sequences,
                        mmseqs_db=mmseqs_path)

    return database
