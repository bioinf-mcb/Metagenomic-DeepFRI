import pandas as pd
from CONFIG import *
from utils import run_command, add_path_to_env


def run_foldseek_search(query_file, work_path):
    target_databases = sorted(list(FOLDSEEK_DATABASES_PATH.iterdir()))[-1]

    add_path_to_env(FOLDSEEK_BIN_PATH)
    run_command(f"foldseek createdb {query_file} {work_path / 'queryDB'}")
    run_command(f"foldseek search {work_path / 'queryDB'} {target_databases / TARGET_DB_NAME} {work_path / 'resultDB'} {TMP_PATH}")
    run_command(f"foldseek convertalis {work_path / 'queryDB'} {target_databases / TARGET_DB_NAME} {work_path / 'resultDB'} {work_path / 'resultDB.m8'}")

    column_names = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits"]
    output_data = pd.read_csv(work_path / 'resultDB.m8', sep="\t", names=column_names)

    # work around pandas Automatic exclusion of “nuisance” columns
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#other-useful-features
    output_data["original_index"] = output_data.index
    grouped = output_data.groupby(["query"], sort=False)
    maximums = grouped.max('fident')
    output_data = output_data.loc[maximums["original_index"]]

    output_data = output_data[output_data.identity > 0.95]
    return output_data