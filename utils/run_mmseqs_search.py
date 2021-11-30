import pandas as pd
from CONFIG import *
from utils.mmseqs_utils import mmseqs_createdb
from utils.mmseqs_utils import mmseqs_search
from utils.mmseqs_utils import mmseqs_convertalis
from utils.mmseqs_utils import MMSEQS_COLUMN_NAMES


def run_mmseqs_search(query_file, work_path):
    target_database_path = sorted(list(MMSEQS_DATABASES_PATH.iterdir()))[-1]

    query_db = work_path / 'queryDB'
    target_db = target_database_path / TARGET_DB_NAME
    result_db = work_path / 'resultDB'
    output_file = work_path / 'resultDB.m8'

    mmseqs_createdb(query_file, query_db)
    mmseqs_search(query_db, target_db, result_db)
    mmseqs_convertalis(query_db, target_db, result_db, output_file)

    return pd.read_csv(output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)
