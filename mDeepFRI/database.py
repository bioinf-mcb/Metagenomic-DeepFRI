import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from mDeepFRI.mmseqs import (QueryFile, _createdb, _createindex,
                             extract_fasta_foldcomp)

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
    else:
        logger.info("Database created at %s", output_path)
        logger.info("Found %s in %s", mmseqs_path, output_path)
        logger.info("Skipping creation of MMSeqs2 database.")

    database = Database(foldcomp_db=input_path,
                        sequence_db=output_sequences,
                        mmseqs_db=mmseqs_path)

    return database


def search_database(query_file: str,
                    database: str,
                    ids: Iterable[str] = None,
                    min_seq_len: int = None,
                    max_seq_len: int = None,
                    max_eval: float = 10e-5,
                    min_bits: float = 0,
                    min_ident: float = 0.5,
                    min_coverage: float = 0.9,
                    top_k: int = 5,
                    threads: int = 1):
    """
    Search databases for sequences similar to the query sequences.

    Args:
        query_file (str): path to the query file in FASTA format.
        database (str): path to MMSeqs database.
        min_seq_len (int): minimum sequence length.
        max_seq_len (int): maximum sequence length.
        max_eval (float): maximum e-value.
        min_bits (float): minimum bitscore.
        min_ident (float): minimum identity.
        top_k (int): number of top hits to return.
        skip_pdb (bool): skip PDB100 database.
        threads (int): number of threads to use.

    """

    query_file = QueryFile(query_file)
    if min_seq_len or max_seq_len or ids:
        query_file.load_sequences(ids=ids)
        query_file.filter_sequences(min_seq_len, max_seq_len)

    results = query_file.search(database, eval=max_eval, threads=threads)
    results.apply_filters(min_cov=min_coverage,
                          min_bits=min_bits,
                          min_ident=min_ident)
    best_hits = results.select_best_matches(top_k, threads=threads)

    return best_hits
