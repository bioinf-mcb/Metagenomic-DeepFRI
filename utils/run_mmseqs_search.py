import pandas as pd

from utils.mmseqs_utils import *


def run_mmseqs_search(query_file, target_db, work_path):
    output_file = work_path / 'mmseqs2_search_results.m8'

    if not output_file.exists():
        query_db = work_path / 'queryDB'
        mmseqs_createdb(query_file, query_db)
        
        result_db = work_path / 'search_resultDB'
        mmseqs_search(query_db, target_db, result_db)
        mmseqs_convertalis(query_db, target_db, result_db, output_file)

    return pd.read_csv(output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)
