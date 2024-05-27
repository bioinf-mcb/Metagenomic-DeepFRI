import logging
import sys
from dataclasses import dataclass
from pathlib import Path

from mDeepFRI.mmseqs import _createdb, _createindex, extract_fasta_foldcomp

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    '[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


@dataclass
class Database:
    """
    Class for storing database paths.
    """
    foldcomp_db: Path
    sequence_db: Path
    mmseqs_db: Path
    mmseqs_result: Path = None

    def __post_init__(self):
        self.foldcomp_db = Path(self.foldcomp_db)
        self.sequence_db = Path(self.sequence_db)
        self.mmseqs_db = Path(self.mmseqs_db)
        if self.mmseqs_result:
            self.mmseqs_result = Path(self.mmseqs_result)
        self.name = self.sequence_db.stem.rsplit(".", 1)[0]


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
    if overwrite or needs_new_mmseqs:
        logger.info("Creating and indexing MMSeqs2 daStabase.")
        _createdb(output_sequences, mmseqs_path)
        _createindex(mmseqs_path, threads)
        logger.info("Database created at %s", output_path)
    else:
        logger.info("Found %s in %s", mmseqs_path, output_path)
        logger.info("Skipping creation of MMSeqs2 database.")

    database = Database(foldcomp_db=input_path,
                        sequence_db=output_sequences,
                        mmseqs_db=mmseqs_path)

    return database
