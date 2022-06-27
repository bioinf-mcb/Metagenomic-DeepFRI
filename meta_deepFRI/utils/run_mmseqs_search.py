import pandas as pd

from CONFIG.FOLDER_STRUCTURE import MMSEQS_SEARCH_RESULTS
from meta_deepFRI.utils.encode_sequence_ids import encode_faa_ids
from meta_deepFRI.utils import mmseqs

MMSEQS_COLUMN_NAMES = [
    "query", "target", "identity", "alignment_length", "mismatches", "gap_openings", "query_start", "query_end",
    "target_start", "target_end", "e_value", "bit_score"
]


def run_mmseqs_search(query_file, target_db, job_path):
    output_file = job_path / MMSEQS_SEARCH_RESULTS

    if output_file.exists():
        return pd.read_csv(output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)

    faa_hashed_ids_path, hash_lookup_dict = encode_faa_ids(query_file)

    query_db = job_path / 'queryDB'
    mmseqs.createdb(faa_hashed_ids_path, query_db)

    result_db = job_path / 'search_resultDB'
    mmseqs.search(query_db, target_db, result_db)

    hashed_output_file = str(output_file) + ".hashed_ids"
    mmseqs.convertalis(query_db, target_db, result_db, hashed_output_file)

    output = pd.read_csv(hashed_output_file, sep="\t", names=MMSEQS_COLUMN_NAMES)
    output["query"] = output["query"].map(hash_lookup_dict)

    output.to_csv(output_file, sep="\t", index=False, header=False)
    return output
