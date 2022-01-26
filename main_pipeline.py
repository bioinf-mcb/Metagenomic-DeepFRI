import argparse
import json
import multiprocessing
import os
import shutil

from Bio import SeqIO

from CONFIG.FOLDER_STRUCTURE import QUERY_PATH, WORK_PATH, FINISHED_PATH, DEFAULT_NAME, MERGED_SEQUENCES, TASK_CONFIG, \
    JOB_CONFIG, ALIGNMENTS, MMSEQS_SEARCH_RESULTS, PROJECT_CONFIG

from CONFIG.RUNTIME_PARAMETERS import CPU_COUNT

from CONFIG.get_config_dict import runtime_config

from metagenomic_deepfri import metagenomic_deepfri
from utils.elapsed_time_logger import ElapsedTimeLogger
from utils.pipeline_utils import find_target_database, load_deepfri_config
from utils.utils import create_unix_timestamp_folder, merge_files_binary, search_files_in_paths, chunks, parse_input_paths


def parse_args():
    parser = argparse.ArgumentParser(description="This script contains all the logic to run DeepFri's GCN or CNN experiments.")

    parser.add_argument("-p", "--project_name", required=False, default=DEFAULT_NAME, help="Task name")

    parser.add_argument("-i", "--input", nargs='+', required=False, default=None,
                        help=f"List of folder or file paths containing query .faa files. Both absolute and relative to {QUERY_PATH} are accepted."
                             f"If not provided pipeline will search in {QUERY_PATH}/--project_name. "
                             f"Use '--input .' to process all files within {QUERY_PATH}")  # logic described here is implemented in parse_input_paths
    # todo add possibility to use specific target_database path and timestamp instead of name only
    parser.add_argument("-t", "--target_db_name", required=False, default=None,
                        help="Target database name or relative path. Will use --project_name if not provided or DEFAULT_NAME if --project_name db is missing.")

    parser.add_argument("-d", "--delete_query", action="store_true", help="Use this flag so that query files are deleted from --input after being copied to project workspace")
    parser.add_argument("-n", "--n_parallel_jobs", required=False, default=1, type=int, help="Number of parallel jobs")
    return parser.parse_args()


# prepare_task validates input data and runtime data used by pipeline
# it creates a new task in WORK_PATH / project_name / timestamp. Timestamp is the task
# task_work_path contains 2 important files MERGED_SEQUENCES and TASK_CONFIG
# It also contains folder with copies of all the query files used in the processing
def prepare_task(project_name, input_paths, target_db_name, delete_query, n_parallel_jobs):
    # verify if deepfri model weights exists
    _ = load_deepfri_config()
    # find and verify target mmseqs database
    target_db = find_target_database(target_db_name)
    print(f"Using target database {target_db}")

    # search for query .faa files and remove empty ones
    query_faa_files = search_files_in_paths(input_paths, ".faa")
    empty_query_files = list(filter(lambda x: x.stat().st_size == 0, query_faa_files))
    if len(empty_query_files) > 0:
        print(f"Found {len(empty_query_files)} empty files.")
        for empty_file in empty_query_files:
            print(f"\t{empty_file}")
            query_faa_files.remove(empty_file)

    # check if there is anything to process
    assert len(query_faa_files) > 0, f"No query .faa files found inside {[str(x) for x in input_paths]}"
    print(f"Query files to be processed: {len(query_faa_files)}")
    for file in query_faa_files:
        print(f"\t{file}")

    # At this point the input data is verified. Create a new task_work_path
    project_work_directory = WORK_PATH / project_name
    project_work_directory.mkdir(parents=True, exist_ok=True)
    task_work_path = create_unix_timestamp_folder(project_work_directory)
    print("Task work path: ", task_work_path)

    # merge sequences from all the query files
    merged_queries_file = task_work_path / MERGED_SEQUENCES
    merge_files_binary(query_faa_files, merged_queries_file)
    # copy query files into task_work_path to have them saved in the results directory
    (task_work_path / "query_files").mkdir()
    for query_faa_file in query_faa_files:
        os.system(f"cp {query_faa_file} {task_work_path / 'query_files'}")
    # delete query files from input if specified
    if delete_query:
        for query_path in query_faa_files:
            query_path.unlink()

    # check if PROJECT_CONFIG exists in  WORK_PATH / project_name
    if (project_work_directory / PROJECT_CONFIG).exists():
        print(f"Using existing PROJECT_CONFIG {project_work_directory / PROJECT_CONFIG}")
        config = json.load(open(project_work_directory / PROJECT_CONFIG))
    else:
        print(f"Creating a new PROJECT_CONFIG {project_work_directory / PROJECT_CONFIG}")
        config = runtime_config()
        json.dump(config, open(project_work_directory / PROJECT_CONFIG, "w"), indent=4)
    for k, v in config.items():
        print(f"\t{v} - {k}")

    # add task specific keys to the task config.
    config.update({
        "project_name": project_name,
        "target_db": str(target_db),
        "target_db_name": target_db_name,
        "timestamp": task_work_path.name,
        "n_parallel_jobs": n_parallel_jobs,
        "query_files": [str(x) for x in query_faa_files]
    })
    # TASK_CONFIG is later copied into job subdirectories as JOB_CONFIG.
    # this file is also used to find tasks inside WORK_PATH to resume them in case of interruption like power outage
    json.dump(config, open(task_work_path / TASK_CONFIG, "w"), indent=4)
    print("Task preparation finished")
    return task_work_path


# jobs are the subdirectories of the task_work_path named 1,2,3...parallel_jobs. 
# metagenomic_deepfri.py runs inside those directories.
# it clones task config as a job config which is used to find jobs inside the task
# query sequences are divided evenly across all jobs
def split_task_into_jobs(task_work_path):
    assert task_work_path / MERGED_SEQUENCES, f"Missing {task_work_path / MERGED_SEQUENCES}"
    assert task_work_path / TASK_CONFIG, f"Missing {task_work_path / TASK_CONFIG}"
    parallel_jobs = json.load(open(task_work_path / TASK_CONFIG))["n_parallel_jobs"]

    with open(task_work_path / MERGED_SEQUENCES, "r") as f:
        query_records = [record for record in SeqIO.parse(f, "fasta")]
    jobs_records = chunks(query_records, parallel_jobs)

    print(f"Dividing {len(query_records)} records across {parallel_jobs}. {len(jobs_records[0])} per job")
    for i in range(len(jobs_records)):
        if len(jobs_records[i]) == 0:
            continue
        job_path = (task_work_path / str(i))
        if job_path.exists():
            shutil.rmtree(job_path)
        job_path.mkdir()
        os.system(f"cp {task_work_path / TASK_CONFIG} {job_path / JOB_CONFIG}")
        with open(job_path / "job_sequences.faa", "w") as f:
            for record in jobs_records[i]:
                f.write(f">{record.id}\n{record.seq}\n")


# search for job config files and start metagenomic_deepfri in parallel
def run_parallel_pipelines(task_work_path):
    timer = ElapsedTimeLogger(task_work_path / "metadata_total_task_time.csv")
    job_paths = [job.parent for job in task_work_path.glob(f"**/{JOB_CONFIG}")]
    print(f"Running {len(job_paths)} jobs on {min(len(job_paths), CPU_COUNT)} parallel threads\n")
    with multiprocessing.Pool(min(len(job_paths), CPU_COUNT)) as p:
        p.map(metagenomic_deepfri, job_paths)
    timer.log("metagenomic_deepfri")


# after completing jobs combine all interesting files into task results
def merge_completed_task_results(task_work_path):
    assert task_work_path / TASK_CONFIG, f"Missing {task_work_path / TASK_CONFIG}"
    task_config = json.load(open(task_work_path / TASK_CONFIG))

    finished_path = FINISHED_PATH / task_config['project_name'] / task_config["timestamp"]
    print("Saving output files to ", finished_path)
    finished_path.mkdir(parents=True)

    os.system(f"cp {task_work_path / TASK_CONFIG} {finished_path}")
    os.system(f"cp {task_work_path / MERGED_SEQUENCES} {finished_path}")
    os.system(f"cp -r {task_work_path}/query_files {finished_path}")

    merge_files_binary(list(task_work_path.glob(f"*/{MMSEQS_SEARCH_RESULTS}")), finished_path / MMSEQS_SEARCH_RESULTS)

    alignments = {}
    for alignment_file in list(task_work_path.glob(f"*/{ALIGNMENTS}")):
        alignments.update(json.load(open(alignment_file, "r")))
    json.dump(alignments, open(finished_path / ALIGNMENTS, "w"), indent=4, sort_keys=True)

    # deepfri result files contain 2 header lines
    for deepfri_result in list(task_work_path.glob("*/results*")):
        if not (finished_path / deepfri_result.name).exists():
            os.system(f"cp {deepfri_result} {finished_path}")
        else:
            with open(finished_path / deepfri_result.name, "a") as dst:
                with open(deepfri_result, "r") as source:
                    dst.writelines(source.readlines()[2:])

    # jobs_metadata_file.parent.name is job index - 1,2,3...parallel_jobs
    for jobs_metadata_file in list(task_work_path.glob("*/metadata*")):
        metadata_store_path = finished_path / "jobs_metadata" / jobs_metadata_file.parent.name
        metadata_store_path.mkdir(parents=True, exist_ok=True)
        os.system(f"cp {jobs_metadata_file} {metadata_store_path}")

    gcn_and_cnn_count = {"GCN": 0, "CNN": 0}
    for jobs_metadata_file in list(task_work_path.glob("*/metadata_cnn_gcn_counts.json")):
        job_metadata = json.load(open(jobs_metadata_file))
        for model in ["GCN", "CNN"]:
            gcn_and_cnn_count[model] += job_metadata[model]
    json.dump(gcn_and_cnn_count, open(finished_path / "metadata_cnn_gcn_counts.json", "w"), indent=4)

    os.system(f"cp {task_work_path}/metadata* {finished_path}")
    print(f"FINISHED {task_config['project_name']}/{task_config['timestamp']}")


# main() simply executes functions implemented above step by step
def main(project_name: str, input: list, target_db_name: str, delete_query: bool, parallel_jobs: int):
    task_work_path = prepare_task(project_name, input, target_db_name, delete_query, parallel_jobs)

    split_task_into_jobs(task_work_path)

    run_parallel_pipelines(task_work_path)

    merge_completed_task_results(task_work_path)

    shutil.rmtree(task_work_path)


if __name__ == '__main__':
    args = parse_args()

    project_name = args.project_name
    delete_query = args.delete_query
    n_parallel_jobs = args.n_parallel_jobs

    if args.target_db_name is None:
        try:
            _ = find_target_database(project_name)
            target_db_name = project_name
        except AssertionError:
            print(f"No {project_name} target database found. Using {DEFAULT_NAME} instead.")
            target_db_name = DEFAULT_NAME
    else:
        # todo add possibility to use specific target_database path and timestamp instead of name only
        target_db_name = args.target_db_name

    input_paths = parse_input_paths(args.input, project_name, QUERY_PATH)

    main(project_name, input_paths, target_db_name, delete_query, n_parallel_jobs)
