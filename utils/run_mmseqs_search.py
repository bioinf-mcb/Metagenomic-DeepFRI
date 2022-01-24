import pandas as pd

from CONFIG.FOLDER_STRUCTURE import MMSEQS_SEARCH_RESULTS
from utils.mmseqs_utils import *


MMSEQS_COLUMN_NAMES = ["query", "target", "identity", "alignment_length", "mismatches", "gap_openings", "query_start",
                       "query_end", "target_start", "target_end", "e_value", "bit_score"]


def run_mmseqs_search(query_file, target_db, job_path):
    output_file = job_path / MMSEQS_SEARCH_RESULTS

    if not output_file.exists():
        query_db = job_path / 'queryDB'
        mmseqs_createdb(query_file, query_db)
        
        result_db = job_path / 'search_resultDB'
        mmseqs_search(query_db, target_db, result_db)
        mmseqs_convertalis(query_db, target_db, result_db, output_file)

    return pd.read_csv(output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)
