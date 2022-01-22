import argparse
import multiprocessing
import os
import pathlib
import shutil
from itertools import repeat

from Bio import SeqIO

from CONFIG.FOLDER_STRUCTURE import QUERY_PATH, WORK_PATH, FINISHED_PATH, \
    DEEPFRI_MODEL_WEIGHTS_JSON_FILE, DEFAULT_NAME, MMSEQS_DATABASES_PATH
from CONFIG.RUNTIME_PARAMETERS import ANGSTROM_CONTACT_THRESHOLD, GENERATE_CONTACTS, CPU_COUNT

from metagenomic_deepfri_pipeline import metagenomic_deepfri_pipeline
from utils.elapsed_time_handler import ElapsedTimeHandler
from utils.pipeline_utils import split_query_into_jobs, merge_result_files
from utils.utils import create_unix_time_folder, merge_files_binary, search_files_in_paths


def parse_args():
    # todo add description
    parser = argparse.ArgumentParser(description="main pipeline")
    parser.add_argument("-n", "--task_name", required=False, default=DEFAULT_NAME, help="Task name")
    parser.add_argument("-q", "--query_paths", nargs='+', required=False, default=None,
                        help=f"Folders paths containing query .faa files and/or paths to .faa files. "
                             f"If not provided pipeline will search in {QUERY_PATH}/task_name. "
                             f"Use '-q all' to process all files within {QUERY_PATH}")
    parser.add_argument("-t", "--target_db_name", required=False, default=DEFAULT_NAME, help="Target database name")
    parser.add_argument("-d", "--delete_query", action="store_true",
                        help="Use this flag so that query files are deleted from --query_paths after being copied to task workspace")
    parser.add_argument("-p", "--parallel_jobs", required=False, default=1, type=int, help="Number of parallel jobs")
    return parser.parse_args()


def main(task_name: str, query_paths: list, target_db_name: str, delete_query: bool, parallel_jobs: int):
    # check if deepfri model weights exists
    assert DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists(), \
        f"No DeepFri model weights json file found at {DEEPFRI_MODEL_WEIGHTS_JSON_FILE}" \
        f"Please run post_setup.py script to download and unzip model weights and json file"
    # check if mmseqs database exists
    assert len(list((MMSEQS_DATABASES_PATH / target_db_name).iterdir())) > 0, \
        f"No target database found {MMSEQS_DATABASES_PATH / target_db_name}. " \
        f"Create one using update_mmseqs_database.py {'' if target_db_name == DEFAULT_NAME else '--name ' + target_db_name} script"

    query_faa_files = search_files_in_paths(query_paths, ".faa")
    assert len(query_faa_files) > 0, f"No protein sequences .faa files found inside {[str(x) for x in query_paths]}"

    print(f"Query files to be processed: {len(query_faa_files)}")
    for file in query_faa_files:
        print(f"\t{file}")

    # create a new work_path for this task
    task_work_path = WORK_PATH / task_name
    task_work_path.mkdir(parents=True, exist_ok=True)
    work_path = create_unix_time_folder(task_work_path)
    print("Work path: ", work_path)
    timer = ElapsedTimeHandler(work_path / "metadata_total_pipeline_time.csv")

    # merge sequences from all the query files
    merged_queries_file = work_path / 'merged_query_sequences.faa'
    merge_files_binary(query_faa_files, merged_queries_file)
    # copy query files into work_path to have them saved in the results directory
    (work_path / "query_files").mkdir()
    for query_faa_file in query_faa_files:
        os.system(f"cp {query_faa_file} {work_path / 'query_files'}")
    # delete query files from query_paths if specified
    if delete_query:
        for query_path in query_faa_files:
            query_path.unlink()

    # split query sequences across parallel_jobs
    with open(merged_queries_file, "r") as f:
        query_records = [record for record in SeqIO.parse(f, "fasta")]
    job_paths = split_query_into_jobs(query_records, work_path, parallel_jobs)
    jobs_args = zip(repeat(target_db_name), job_paths, repeat(ANGSTROM_CONTACT_THRESHOLD), repeat(GENERATE_CONTACTS))
    timer.log("data_preparation")

    # run metagenomic_deepfri_pipeline in parallel
    with multiprocessing.Pool(min(parallel_jobs, CPU_COUNT)) as p:
        p.starmap(metagenomic_deepfri_pipeline, jobs_args)
    timer.log("metagenomic_deepfri_pipeline")

    # merge jobs results and store them in finished_path
    finished_path = FINISHED_PATH / task_name / work_path.name
    print("Finished! Saving output files to ", finished_path)
    merge_result_files(work_path, finished_path)
    timer.log_total_time()
    os.system(f"cp {work_path}/metadata* {finished_path}")

    shutil.rmtree(work_path)


if __name__ == '__main__':
    args = parse_args()

    task_name = args.task_name
    target_db_name = args.target_db_name
    delete_query = args.delete_query
    parallel_jobs = args.parallel_jobs

    if args.query_paths is None:
        query_paths = [pathlib.Path(QUERY_PATH / task_name)]
    else:
        if args.query_paths == ["all"]:
            query_paths = [QUERY_PATH]
        else:
            query_paths = [pathlib.Path(x) for x in args.query_paths]

    main(task_name, query_paths, target_db_name, delete_query, parallel_jobs)
