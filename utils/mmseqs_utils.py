import tempfile

from utils.utils import run_command


MMSEQS_COLUMN_NAMES = ["query", "target", "identity", "alignment_length", "mismatches", "gap_openings", "query_start",
                "query_end", "target_start", "target_end", "e_value", "bit_score"]


def mmseqs_createdb(sequences_file, db_path):
    run_command(f"mmseqs createdb {sequences_file} {db_path} --dbtype 1")


def mmseqs_createindex(db_path):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(f"mmseqs createindex {db_path} {tmp_path}")


def mmseqs_search(query_db, target_db, result_db):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(f"mmseqs search {query_db} {target_db} {result_db} {tmp_path}")


def mmseqs_convertalis(query_db, target_db, result_db, output_file):
    run_command(f"mmseqs convertalis {query_db} {target_db} {result_db} {output_file}")
