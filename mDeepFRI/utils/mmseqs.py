import logging
import tempfile
from pathlib import Path

import pandas as pd

import mDeepFRI
from mDeepFRI import MMSEQS_SEARCH_RESULTS
from mDeepFRI.utils.utils import run_command

MMSEQS_COLUMN_NAMES = [
    "query", "target", "identity", "alignment_length", "mismatches",
    "gap_openings", "query_start", "query_end", "target_start", "target_end",
    "e_value", "bit_score"
]

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def createdb(sequences_file, db_path):
    """
    Converts FASTA file to a DB format needed for MMseqs2.
    This should generate five files,
    e.g. queryDB, queryDB_h and its corresponding index file queryDB.index,
    queryDB_h.index and queryDB.lookup from the FASTA QUERY.fasta input sequences.

    sequence_file (str): path to FASTA file.
    db_path (str): path to output db file.

    Returns:
        None
    """
    run_command(f"mmseqs createdb {sequences_file} {db_path} --dbtype 1")


def createindex(db_path):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(f"mmseqs createindex {db_path} {tmp_path}")


def search(query_db, target_db, result_db):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(
            f"mmseqs search {query_db} {target_db} {result_db} {tmp_path}")


def convertalis(query_db, target_db, result_db, output_file):
    run_command(
        f"mmseqs convertalis {query_db} {target_db} {result_db} {output_file}")


def extract_fasta_foldcomp(foldcomp_db: str,
                           output_file: str,
                           threads: int = 1):
    """
    Extracts FASTA from database
    """
    foldcomp_bin = Path(mDeepFRI.__path__[0]).parent / "foldcomp"

    # run command
    run_command(
        f"{foldcomp_bin} extract --fasta -t {threads} {foldcomp_db} {output_file}"
    )


def create_target_database(foldcomp_fasta_path: Path,
                           mmseqs_db_path: Path) -> None:
    """
    Extracts sequences from compressed FoldComp database.

    Args:
        foldcomp_db_path (pathlib.Path): Path to FoldComp database.
        new_db_path (pathlib.Path): Path to new MMSeqs database.

    Returns:
        None
    """
    createdb(foldcomp_fasta_path, mmseqs_db_path)
    logging.info("Indexing new target mmseqs2 database %s",
                 str(mmseqs_db_path))
    createindex(mmseqs_db_path)


def run_mmseqs_search(query_file: Path, target_db: Path, output_path: Path,
                      min_bit_score: float, max_evalue: float,
                      min_identity: float) -> pd.DataFrame:
    """Creates a database from query sequences and runs mmseqs2 search against database.

    Args:
        query_file (pathlib.Path): Path to query FASTA file.
        target_db (pathlib.Path): Path to target MMSeqs2 database.
        output_path (pathlib.Path):

    Returns:
        pd.DataFrame: Pandas DataFrame with MMSeqs2 search results.
    """
    output_file = output_path / MMSEQS_SEARCH_RESULTS

    query_db = output_path / 'queryDB'
    createdb(query_file, query_db)

    with tempfile.TemporaryDirectory() as tmp_path:

        result_db = Path(tmp_path) / 'search_resultDB'
        search(query_db, target_db, result_db)

        # Convert results to tabular format
        convertalis(query_db, target_db, result_db, output_file)

    output = pd.read_csv(output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)

    # MMSeqs2 alginment filters
    if min_identity:
        output = output.query(f"identity >= {min_identity}")
    if min_bit_score:
        output = output.query(f"bitscore >= {min_bit_score}")
    if max_evalue:
        output = output.query(f"e_value <= {max_evalue}")

    output.to_csv(output_file, sep="\t", index=False)

    return output
