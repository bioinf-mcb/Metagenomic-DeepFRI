from pathlib import Path
import tempfile
import pandas as pd

from meta_deepFRI.config.names import MERGED_SEQUENCES, TARGET_MMSEQS_DB_NAME, MMSEQS_SEARCH_RESULTS
from meta_deepFRI.utils.utils import run_command, merge_files_binary

MMSEQS_COLUMN_NAMES = [
    "query", "target", "identity", "alignment_length", "mismatches", "gap_openings", "query_start", "query_end",
    "target_start", "target_end", "e_value", "bit_score"
]


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
        run_command(f"mmseqs search {query_db} {target_db} {result_db} {tmp_path}")


def convertalis(query_db, target_db, result_db, output_file):
    run_command(f"mmseqs convertalis {query_db} {target_db} {result_db} {output_file}")


def create_target_database(seq_atoms_path: Path, new_db_path: Path) -> None:
    """

    :param seq_atoms_path:
    :param new_db_path:
    :param freshly_added_ids:
    :return:
    """
    sequence_files = list((seq_atoms_path).glob("**/*.faa"))
    print("\nMerging " + str(len(sequence_files)) + " sequence files for mmseqs2")
    merge_files_binary(sequence_files, seq_atoms_path / MERGED_SEQUENCES)

    print("Creating new target mmseqs2 database " + str(new_db_path))
    createdb(seq_atoms_path / MERGED_SEQUENCES, new_db_path / TARGET_MMSEQS_DB_NAME)
    print("Indexing new target mmseqs2 database " + str(new_db_path))
    createindex(new_db_path / TARGET_MMSEQS_DB_NAME)


def run_mmseqs_search(query_file: Path, target_db: Path, output_path: Path) -> pd.DataFrame:
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

    output.to_csv(output_file, sep="\t", index=False)

    return output
