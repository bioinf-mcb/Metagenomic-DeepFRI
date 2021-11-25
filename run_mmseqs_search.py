import pandas as pd
from CONFIG import *
from utils import run_command


def run_mmseqs_search(query_file, work_path):
    target_databases = sorted(list(MMSEQS_DATABASES_PATH.iterdir()))[-1]

    run_command(f"mmseqs createdb {query_file} {work_path / 'queryDB'} --dbtype 1")
    run_command(f"mmseqs search {work_path / 'queryDB'} {target_databases / TARGET_DB_NAME} {work_path / 'resultDB'} {TMP_PATH}")
    run_command(f"mmseqs convertalis {work_path / 'queryDB'} {target_databases / TARGET_DB_NAME} {work_path / 'resultDB'} {work_path / 'resultDB.m8'}")

    column_names = ["query", "target", "identity", "alignment_length", "mismatches", "gap_openings", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score"]
    output_data = pd.read_csv(work_path / 'resultDB.m8', sep="\t", names=column_names)

    # work around pandas Automatic exclusion of “nuisance” columns
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#other-useful-features
    output_data["original_index"] = output_data.index
    grouped = output_data.groupby(["query"], sort=False)
    maximums = grouped.max('identity')
    output_data = output_data.loc[maximums["original_index"]]

    # output_data = output_data[output_data.identity > 0.95]
    return output_data