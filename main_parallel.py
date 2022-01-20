import argparse
import multiprocessing
import os
import shutil

from itertools import repeat

from CONFIG.FOLDER_STRUCTURE import QUERY_PATH, WORK_PATH, FINISHED_PATH, \
    DEEPFRI_MODEL_WEIGHTS_JSON_FILE, DEFAULT_TARGET_DB_NAME
from CONFIG.RUNTIME_PARAMETERS import ANGSTROM_CONTACT_THRESHOLD, GENERATE_CONTACTS, CPU_COUNT

from metagenomic_deepfri_pipeline import metagenomic_deepfri_pipeline
from utils.elapsed_time_handler import ElapsedTimeHandler
from utils.utils import create_unix_time_folder


def parse_args():
    # todo add description
    parser = argparse.ArgumentParser(description="main pipeline")

    parser.add_argument("-t", "--target_db", required=False, default=DEFAULT_TARGET_DB_NAME, help="Target database name")
    return parser.parse_args()


def job(work_path, target_db_name):
    metagenomic_deepfri_pipeline(target_db_name, work_path, ANGSTROM_CONTACT_THRESHOLD, GENERATE_CONTACTS)
    finished_path = FINISHED_PATH / work_path.name
    print("Finished! Saving output files to ", finished_path)
    finished_path.mkdir(parents=True, exist_ok=True)
    os.system(f"cp {work_path}/*.faa {finished_path}")
    os.system(f"cp {work_path}/results* {finished_path}")
    os.system(f"cp {work_path}/alignments.json {finished_path}")
    os.system(f"cp {work_path}/mmseqs2_search_results.m8 {finished_path}")
    os.system(f"cp {work_path}/metadata* {finished_path}")

    shutil.rmtree(work_path)


def main():
    args = parse_args()
    timer = ElapsedTimeHandler()
    if not DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists():
        print("Please run post_setup.py script to download and unzip model weights")
        exit(1)

    # todo fancier way of selecting query files
    query_faa_files = list(QUERY_PATH.glob("**/*.faa"))
    if len(query_faa_files) == 0:
        print(f"No protein sequences files found in {QUERY_PATH}. Terminating")
        exit(0)
    print(f"Query files to be processed: {len(query_faa_files)}")
    for file in query_faa_files:
        print(f"\t{file}")

    work_paths = []
    for query_file in query_faa_files:
        work_path = create_unix_time_folder(WORK_PATH)
        print("Work path: ", work_path)
        os.system(f"cp {query_file} {work_path}")
        work_paths.append(work_path)
    # todo remove query_faa_files from QUERY_PATH so they don't get processed again later

    with multiprocessing.Pool(min(len(query_faa_files), CPU_COUNT)) as p:
        p.starmap(job, zip(work_paths, repeat(args.target_db)))
    timer.log_total_time()


if __name__ == '__main__':
    main()
