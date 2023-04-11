import hashlib
import pathlib
import tempfile

from typing import Tuple

import pandas as pd

from meta_deepFRI.config.names import MERGED_SEQUENCES, TARGET_MMSEQS_DB_NAME, MMSEQS_SEARCH_RESULTS
from meta_deepFRI.utils.utils import run_command, merge_files_binary
from .fasta_file_io import load_fasta_file, write_fasta_file

MMSEQS_COLUMN_NAMES = [
    "query", "target", "identity", "alignment_length", "mismatches", "gap_openings", "query_start", "query_end",
    "target_start", "target_end", "e_value", "bit_score"
]


def hash_sequence_id(sequence: str) -> str:
    """Return the SHA256 encoding of protein sequence

    Args:
        sequence (str): Aminoacid sequence of the protein.

    Returns:
        SHA256 encoding of the protein sequence.
    """
    return hashlib.sha256(bytes(sequence, encoding='utf-8')).hexdigest()


def encode_faa_ids(input_file: pathlib.Path, output_file: pathlib.Path) -> Tuple[str, dict]:
    """
    Encodes multiline FASTA file IDs with SHA256 encoding.

    Args:
        input_file (pathlib.Path): Path to a fasta file.
        output_file (pathlib.Path): Path to the output file.

    Returns:
        Path to the new file with ID and respective SHA256-encoded sequences and
        a dictionary where keys represent sequence ID, and values SHA256 encoding.
    """

    seq_records = load_fasta_file(input_file)
    hash_lookup_dict = {}

    for fasta_entry in seq_records:
        seq_id_hash = hash_sequence_id(fasta_entry.id)
        while seq_id_hash in hash_lookup_dict.keys():
            seq_id_hash = hash_sequence_id(fasta_entry.id + "1")
        hash_lookup_dict[seq_id_hash] = fasta_entry.id
        fasta_entry.id = seq_id_hash

    write_fasta_file(seq_records, output_file)
    return output_file, hash_lookup_dict


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


def create_target_database(seq_atoms_path: pathlib.Path, new_db_path: pathlib.Path) -> None:
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


def run_mmseqs_search(query_file: pathlib.Path, target_db: pathlib.Path, output_path: pathlib.Path) -> pd.DataFrame:
    """Creates a database from query sequences and runs mmseqs2 search against database.

    Args:
        query_file (pathlib.Path): Path to query FASTA file.
        target_db (pathlib.Path): Path to target MMSeqs2 database.
        output_path (pathlib.Path):

    Returns:
        pd.DataFrame: Pandas DataFrame with MMSeqs2 search results.
    """
    output_file = output_path / MMSEQS_SEARCH_RESULTS

    # Check if search results already exist
    if output_file.exists():

        return pd.read_csv(output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)

    else:
        # Create hashed IDs for query sequences
        hashed_file = output_path / query_file.stem + ".hashed_ids.faa"
        faa_hashed_ids_path, hash_lookup_dict = encode_faa_ids(query_file, hashed_file)
        # Create MMSeqs2 database from query sequences.
        query_db = output_path / 'queryDB'
        createdb(faa_hashed_ids_path, query_db)

        result_db = output_path / 'search_resultDB'
        search(query_db, target_db, result_db)

        # Convert results to tabular format
        hashed_output_file = str(output_file) + ".hashed_ids"
        convertalis(query_db, target_db, result_db, hashed_output_file)

        # Decode hashed IDs
        output = pd.read_csv(hashed_output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)
        output["query"] = output["query"].map(hash_lookup_dict)

        output.to_csv(output_file, sep="\t", index=False, header=False)

        return output
