import json
import pathlib
import tempfile

import pandas as pd

from meta_deepFRI.config.names import MERGED_SEQUENCES, TARGET_MMSEQS_DB_NAME, MMSEQS_SEARCH_RESULTS
from meta_deepFRI.utils.utils import run_command, merge_files_binary
from .fasta_file_io import encode_faa_ids

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


def create_target_database(seq_atoms_path: pathlib.Path, new_db_path: pathlib.Path, freshly_added_ids: list) -> None:
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

    print("Saving freshly added sequence ids to " + str(new_db_path / "structure_ids_added.json"))
    json.dump(sorted(freshly_added_ids), open(new_db_path / "structure_ids_added.json", "w"), indent=4)


def run_mmseqs_search(query_file, target_db, job_path):
    output_file = job_path / MMSEQS_SEARCH_RESULTS

    if output_file.exists():
        return pd.read_csv(output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)

    faa_hashed_ids_path, hash_lookup_dict = encode_faa_ids(query_file)

    query_db = job_path / 'queryDB'
    createdb(faa_hashed_ids_path, query_db)

    result_db = job_path / 'search_resultDB'
    search(query_db, target_db, result_db)

    hashed_output_file = str(output_file) + ".hashed_ids"
    convertalis(query_db, target_db, result_db, hashed_output_file)

    output = pd.read_csv(hashed_output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)
    output["query"] = output["query"].map(hash_lookup_dict)

    output.to_csv(output_file, sep="\t", index=False, header=False)
    return output
