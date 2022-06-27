import json
import pathlib
import tempfile

from CONFIG.FOLDER_STRUCTURE import SEQUENCES, MERGED_SEQUENCES, MMSEQS_DATABASES_PATH, TARGET_MMSEQS_DB_NAME
from meta_deepFRI.utils.utils import run_command, merge_files_binary, create_unix_timestamp_folder


def createdb(sequences_file, db_path):
    run_command(f"mmseqs createdb {sequences_file} {db_path} --dbtype 1")


def createindex(db_path):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(f"mmseqs createindex {db_path} {tmp_path}")


def search(query_db, target_db, result_db):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(f"mmseqs search {query_db} {target_db} {result_db} {tmp_path}")


def convertalis(query_db, target_db, result_db, output_file):
    run_command(f"mmseqs convertalis {query_db} {target_db} {result_db} {output_file}")


def create_target_database(seq_atoms_path: pathlib.Path, project_name: str, freshly_added_ids: list) -> None:
    """

    :param seq_atoms_path:
    :param project_name:
    :param freshly_added_ids:
    :return:
    """
    sequence_files = list((seq_atoms_path / SEQUENCES).glob("**/*.faa"))
    print("\nMerging " + str(len(sequence_files)) + " sequence files for mmseqs2")
    merge_files_binary(sequence_files, seq_atoms_path / MERGED_SEQUENCES)

    mmseqs2_db_path = create_unix_timestamp_folder(MMSEQS_DATABASES_PATH / project_name)
    print("Creating new target mmseqs2 database " + str(mmseqs2_db_path))
    createdb(seq_atoms_path / MERGED_SEQUENCES, mmseqs2_db_path / TARGET_MMSEQS_DB_NAME)
    print("Indexing new target mmseqs2 database " + str(mmseqs2_db_path))
    createindex(mmseqs2_db_path / TARGET_MMSEQS_DB_NAME)

    print("Saving freshly added sequence ids to " + str(mmseqs2_db_path / "structure_ids_added.json"))
    json.dump(sorted(freshly_added_ids), open(mmseqs2_db_path / "structure_ids_added.json", "w"), indent=4)
