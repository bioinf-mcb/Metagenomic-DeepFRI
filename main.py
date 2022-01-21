import argparse
import json
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
from utils.utils import create_unix_time_folder, merge_files_binary, split_dict


def parse_args():
    # todo add description
    parser = argparse.ArgumentParser(description="main pipeline")
    parser.add_argument("-n", "--task_name", required=False, default=DEFAULT_NAME, help="Task name")
    parser.add_argument("-q", "--query_paths", required=False, default=None,
                        help=f"Folders paths containing query .faa files and/or paths to .faa files. "
                             f"If not provided pipeline will search in {QUERY_PATH}/task_name")
    parser.add_argument("-t", "--target_db_name", required=False, default=DEFAULT_NAME, help="Target database name")
    parser.add_argument("-s", "--store_query", action="store_true",
                        help="Use this flag so that query files are not deleted from --query_paths after copied to task workspace")
    parser.add_argument("-p", "--parallel_jobs", required=False, default=1, type=int, help="Number of parallel jobs")
    parser.add_argument("-i", "--ignore_duplicated_ids", action="store_true",
                        help="Use this flag to force pipeline execution that will process only one of duplicated protein ids")
    return parser.parse_args()


def locate_faa_files(query_paths):
    query_faa_files = []
    for query_path in query_paths:
        if not query_path.exists():
            print(f"Unable to locate {query_path}.")
            continue
        if query_path.is_dir():
            query_faa_files.extend(list(query_path.glob("**/*.faa")))
        else:
            if not query_path.name.endswith(".faa"):
                print(f"{query_path} is not an .faa file which is preferred format.")
            query_faa_files.append(query_path)
    return query_faa_files


def main(task_name: str, query_paths: list, target_db_name: str, store_query: bool, parallel_jobs: int, ignore_duplicated_ids: bool):
    assert DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists(), \
        f"No DeepFri model weights json file found at {DEEPFRI_MODEL_WEIGHTS_JSON_FILE}" \
        f"Please run post_setup.py script to download and unzip model weights and json file"

    assert len(list((MMSEQS_DATABASES_PATH / target_db_name).iterdir())) > 0, \
        f"No target database found {MMSEQS_DATABASES_PATH / target_db_name}. " \
        f"Create one using update_mmseqs_database.py script"

    query_faa_files = locate_faa_files(query_paths)
    assert len(query_faa_files) > 0, f"No protein sequences .faa files found inside {[str(x) for x in query_paths]}"

    print(f"Query files to be processed: {len(query_faa_files)}")
    for file in query_faa_files:
        print(f"\t{file}")

    task_work_path = WORK_PATH / task_name
    task_work_path.mkdir(parents=True, exist_ok=True)
    work_path = create_unix_time_folder(task_work_path)
    print("Work path: ", work_path)
    timer = ElapsedTimeHandler(work_path / "metadata_total_pipeline_time.csv")

    merged_queries_file = work_path / 'merged_query_sequences.faa'
    merge_files_binary(query_faa_files, merged_queries_file)
    json.dump([str(x) for x in query_faa_files], open(work_path / "metadata_source_query_files.json", "w"))

    with open(merged_queries_file, "r") as f:
        query_seqs = [record for record in SeqIO.parse(f, "fasta")]
    unique_seq_id = {record.id: record.seq for record in query_seqs}

    if len(query_seqs) > len(unique_seq_id):
        assert ignore_duplicated_ids, \
            f"There are {len(query_seqs) - len(unique_seq_id)} duplicated IDs across query files. Use --ignore_duplicated_ids"
        with open(merged_queries_file, "w") as f:
            for protein_id in unique_seq_id.keys():
                f.write(f">{protein_id}\n{unique_seq_id[protein_id]}\n")

    jobs_paths = []
    jobs_sequences = list(split_dict(unique_seq_id, max(1, int(len(unique_seq_id) / parallel_jobs))))
    for i in range(len(jobs_sequences)):
        (work_path / str(i)).mkdir()
        jobs_paths.append(work_path / str(i))
        with open(jobs_paths[i] / "jobs_sequences.faa", "w") as f:
            for protein_id in jobs_sequences[i].keys():
                f.write(f">{protein_id}\n{jobs_sequences[i][protein_id]}\n")

    if not store_query:
        for query_path in query_faa_files:
            query_path.unlink()

    jobs_args = zip(repeat(target_db_name), jobs_paths, repeat(ANGSTROM_CONTACT_THRESHOLD), repeat(GENERATE_CONTACTS))
    with multiprocessing.Pool(min(parallel_jobs, CPU_COUNT)) as p:
        p.starmap(metagenomic_deepfri_pipeline, jobs_args)

    timer.log_total_time()

    finished_path = FINISHED_PATH / task_name / work_path.name
    print("Finished! Saving output files to ", finished_path)
    finished_path.mkdir(parents=True)
    os.system(f"cp {work_path}/merged_query_sequences.faa {finished_path}")
    for job_path in jobs_paths:
        job_results_path = finished_path / job_path.name
        job_results_path.mkdir()
        os.system(f"cp {job_path}/jobs_sequences.faa {job_results_path}")
        os.system(f"cp {job_path}/results* {job_results_path}")
        os.system(f"cp {job_path}/alignments.json {job_results_path}")
        os.system(f"cp {job_path}/mmseqs2_search_results.m8 {job_results_path}")
        os.system(f"cp {job_path}/metadata* {job_results_path}")

    shutil.rmtree(work_path)


if __name__ == '__main__':
    args = parse_args()

    task_name = args.task_name
    if args.query_paths is None:
        query_paths = [pathlib.Path(QUERY_PATH / task_name)]
    else:
        query_paths = [pathlib.Path(x) for x in args.query_paths]

    target_db_name = args.target_db_name
    store_query = args.store_query
    parallel_jobs = args.parallel_jobs
    ignore_duplicated_ids = args.ignore_duplicated_ids

    main(task_name, query_paths, target_db_name, store_query, parallel_jobs, ignore_duplicated_ids)
