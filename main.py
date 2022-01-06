import shutil
import time

from CONFIG import *
from metagenomic_deepfri_pipeline import metagenomic_deepfri_pipeline

if __name__ == '__main__':
    if not DEEPFRI_MODEL_WEIGHTS_JSON_PATH.exists():
        print("Please run post_setup.py script to download/unzip model weights")
        exit(1)

    query_faa_files = list(QUERY_PATH.glob("**/*.faa"))
    if len(query_faa_files) == 0:
        print("No query protein files found terminating")
        exit(0)
    print("Query files found: ", query_faa_files)

    target_database_path = sorted(list(MMSEQS_DATABASES_PATH.iterdir()))[-1]
    target_db = target_database_path / TARGET_DB_NAME
    print("Target database: ", target_database_path)

    work_start = str(time.time())
    work_path = (WORK_PATH / work_start)
    while work_path.exists():
        time.sleep(1)
        work_start = str(time.time())
        work_path = (WORK_PATH / work_start)
    work_path.mkdir(parents=True)
    print("Work path: ", work_path)

    query_file = work_path / 'merged_query_sequences.faa'
    with open(query_file, 'wb') as writer:
        for seq_file in query_faa_files:
            with open(seq_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)
            # todo remove query_file

    metagenomic_deepfri_pipeline(query_file, target_db, work_path)
