import argparse
import json
import multiprocessing
import os
import pathlib
import shutil

import dataclasses
from itertools import repeat
from typing import List, Union

##TODO: add threading option
from meta_deepFRI.config.folder_structure import FolderStructureConfig, load_folder_structure_config
from meta_deepFRI.config.names import MERGED_SEQUENCES, PROJECT_CONFIG, JOB_CONFIG, TASK_CONFIG, ALIGNMENTS, \
    MMSEQS_SEARCH_RESULTS, \
    DEFAULT_NAME, TASK_SEQUENCES
from meta_deepFRI.config.runtime_config import load_runtime_config
from meta_deepFRI.config.job_config import JobConfig, load_job_config

from meta_deepFRI.utils.elapsed_time_logger import ElapsedTimeLogger
from meta_deepFRI.utils.fasta_file_io import load_fasta_file, write_fasta_file
from meta_deepFRI.utils.pipeline_utils import find_target_database, load_deepfri_config
from meta_deepFRI.utils.utils import create_unix_timestamp_folder, merge_files_binary, search_files_in_paths, chunks, parse_input_paths
from meta_deepFRI import metagenomic_deepfri


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(description="This script contains all the logic to run DeepFri's GCN or CNN experiments.")

    parser.add_argument("-r", "--data_root", required=False, default="default_config/data_root_path.json",
                        help="Path to json file containing DATA_ROOT or path to folder that will be used as DATA_ROOT")

    # logic described here is implemented in parse_input_paths
    parser.add_argument("-i", "--input", nargs='+', required=False, default=None,
                        help="List of folder or file paths containing query .faa files. Both absolute and relative to FolderStructureConfig.QUERY_PATH are accepted."
                             "If not provided pipeline will search in FolderStructureConfig.QUERY_PATH/project_name. "
                             "Use '--input .' to process all files within FolderStructureConfig.QUERY_PATH")

    parser.add_argument("-t", "--threads", required=False, default=1, type=int,
                        help="Number of threads to use for parallel processing.")

    # todo add possibility to use specific target_database path and timestamp instead of name only
    parser.add_argument("-t", "--target_db_name", required=False, default=None,
                        help="Target database name or relative path. Will use --project_name if not provided or DEFAULT_NAME if --project_name db is missing.")

    parser.add_argument("-n", "--n_parallel_tasks", required=False, default=1, type=int,
                        help="Number of parallel tasks")
    # yapf: enable
    return parser.parse_args()


def find_query_faa_files(input_paths: List[Union[str, pathlib.Path]]) -> List[pathlib.Path]:
    """

    :param input_paths:
    :return:
    """
    # search for query .faa files and remove empty ones
    query_faa_files = search_files_in_paths(input_paths, ".faa")
    empty_query_files = list(filter(lambda x: x.stat().st_size == 0, query_faa_files))
    if len(empty_query_files) > 0:
        print(f"Found {len(empty_query_files)} empty files.")
        for empty_file in empty_query_files:
            print(f"\t{empty_file}")
            query_faa_files.remove(empty_file)

    # check if there is anything to process
    assert len(query_faa_files) > 0, f"No .faa files found inside {[str(x) for x in input_paths]}"
    print(f"Query files to be processed: {len(query_faa_files)}")
    for file in query_faa_files:
        print(f"\t{file}")

    return query_faa_files


def prepare_job(fsc: FolderStructureConfig, project_name: str, input_paths: list, target_db_name,
                n_parallel_tasks) -> pathlib.Path:
    """
    Validates input data and runtime data used by pipeline
    it creates a new job in WORK_PATH / project_name / timestamp. Timestamp is the job directory
    job_work_path contains 2 important files MERGED_SEQUENCES and JOB_CONFIG
    It also contains folder with copies of all the query files used in the processing
    Later in the pipeline additional tasks folders will be created underneath

    :param fsc:
    :param project_name:
    :param input_paths:
    :param target_db_name:
    :param n_parallel_tasks:
    :return:
    """
    # verify if deepfri model weights exists
    _ = load_deepfri_config(fsc)
    # find and verify target mmseqs database
    target_db = find_target_database(fsc, target_db_name)
    print(f"Using target database {target_db}")

    # search for query .faa files and remove empty ones
    query_faa_files = find_query_faa_files(input_paths)

    # At this point the input data is verified. Create a new job_work_path
    project_work_directory = fsc.WORK_PATH / project_name
    job_work_path = create_unix_timestamp_folder(project_work_directory)
    print("Task work path: ", job_work_path)

    # merge sequences from all the query files
    merged_queries_file = job_work_path / MERGED_SEQUENCES
    merge_files_binary(query_faa_files, merged_queries_file)
    # copy query files into job_work_path to have them saved in the results directory
    (job_work_path / "query_files").mkdir()
    for query_faa_file in query_faa_files:
        os.system(f"cp {query_faa_file} {job_work_path / 'query_files'}")

    # Load PROJECT_RUNTIME_CONFIG
    if (project_work_directory / PROJECT_CONFIG).exists():
        print(f"Using existing PROJECT_CONFIG {project_work_directory / PROJECT_CONFIG}")
        project_runtime_config = load_runtime_config(project_work_directory / PROJECT_CONFIG)
    else:
        raise RuntimeError("No project_config found. Please run create_project")

    for k, v in dataclasses.asdict(project_runtime_config).items():
        print(f"\t{v} - {k}")

    # add additional values to job config.
    job_config = JobConfig(
        project_runtime_config,
        project_name=project_name,
        target_db=str(target_db),
        target_db_name=target_db_name,
        timestamp=job_work_path.name,
        n_parallel_tasks=n_parallel_tasks,
        query_files=[str(x) for x in query_faa_files])

    # JOB_CONFIG is later copied into task subdirectories as TASK_CONFIG.
    # this file is also used to find jobs inside WORK_PATH to resume them in case of interruption
    job_config.save(job_work_path / JOB_CONFIG)
    print("Task preparation finished")
    return job_work_path


def split_job_into_tasks(job_work_path: pathlib.Path) -> None:
    """
    Tasks are the subdirectories of the job_work_path named 1,2,3...n_parallel_tasks.
    metagenomic_deepfri.py runs inside those directories.
    it clones job config as a task config which is used to find tasks inside the job
    query sequences are divided evenly across all parallel tasks
    :param job_work_path:
    :return:
    """
    assert job_work_path / MERGED_SEQUENCES, f"Missing {job_work_path / MERGED_SEQUENCES}"
    assert job_work_path / JOB_CONFIG, f"Missing {job_work_path / JOB_CONFIG}"

    job_config = load_job_config(job_work_path / JOB_CONFIG)

    job_sequences = load_fasta_file(job_work_path / MERGED_SEQUENCES)
    tasks_sequences = chunks(job_sequences, job_config.n_parallel_tasks)

    print(
        f"Dividing {len(job_sequences)} seqs across {job_config.n_parallel_tasks} parallel tasks. {len(tasks_sequences[0])} per task"
    )
    for i, task_sequences in enumerate(tasks_sequences):
        if len(task_sequences) == 0:
            continue
        task_path = (job_work_path / str(i))
        if task_path.exists():
            shutil.rmtree(task_path)
        task_path.mkdir()
        job_config.save(task_path / TASK_CONFIG)
        write_fasta_file(task_sequences, task_path / TASK_SEQUENCES)


# search for job config files and start metagenomic_deepfri in parallel
def run_parallel_pipelines(fsc: FolderStructureConfig, job_work_path: pathlib.Path, threads: int) -> None:
    timer = ElapsedTimeLogger(job_work_path / "metadata_total_task_time.csv")

    task_configs_paths = [task_config.parent for task_config in job_work_path.glob(f"**/{TASK_CONFIG}")]

    if len(task_configs_paths) == 1:
        metagenomic_deepfri.metagenomic_deepfri(fsc, task_configs_paths[0])
    else:
        print(f"Running {len(task_configs_paths)} jobs. Parallel: {min(len(task_configs_paths), threads)}\n")
        with multiprocessing.Pool(min(len(task_configs_paths), threads)) as p:
            p.starmap(metagenomic_deepfri.metagenomic_deepfri, zip(repeat(fsc), task_configs_paths))

    timer.log("metagenomic_deepfri")


# after completing jobs combine all interesting files into task results
def merge_completed_task_results(fsc: FolderStructureConfig, job_work_path: pathlib.Path) -> None:
    assert job_work_path / JOB_CONFIG, f"Missing {job_work_path / JOB_CONFIG}"
    job_conf = load_job_config(job_work_path / JOB_CONFIG)

    finished_path = fsc.FINISHED_PATH / job_conf.project_name / job_conf.timestamp
    print("Saving output files to ", finished_path)
    finished_path.mkdir(parents=True)

    os.system(f"cp {job_work_path / JOB_CONFIG} {finished_path}")
    os.system(f"cp {job_work_path / MERGED_SEQUENCES} {finished_path}")
    os.system(f"cp -r {job_work_path}/query_files {finished_path}")

    merge_files_binary(list(job_work_path.glob(f"*/{MMSEQS_SEARCH_RESULTS}")), finished_path / MMSEQS_SEARCH_RESULTS)

    alignments = {}
    for task_alignment_file in list(job_work_path.glob(f"*/{ALIGNMENTS}")):
        alignments.update(json.load(open(task_alignment_file, "r")))
    json.dump(alignments, open(finished_path / ALIGNMENTS, "w"), indent=4, sort_keys=True)

    for pattern in ["csv", "tsv"]:
        for task_deepfri_result in list(job_work_path.glob(f"*/results*{pattern}")):
            if not (finished_path / task_deepfri_result.name).exists():
                os.system(f"cp {task_deepfri_result} {finished_path}")
            else:
                with open(finished_path / task_deepfri_result.name, "a") as dst:
                    with open(task_deepfri_result, "r") as source:
                        dst.writelines(source.readlines()[1:])

    for task_deepfri_result in list(job_work_path.glob("*/results*json")):
        if not (finished_path / task_deepfri_result.name).exists():
            os.system(f"cp {task_deepfri_result} {finished_path}")
        else:
            with open(finished_path / task_deepfri_result.name, "r+") as dst:
                with open(task_deepfri_result, "r") as source:
                    existing = json.load(dst)
                    existing.extend(json.load(source))
                    json.dump(existing, dst, indent=4, sort_keys=True)

    # task_metadata_files.parent.name is job index - 1,2,3...n_parallel_tasks
    for task_metadata_files in list(job_work_path.glob("*/metadata*")):
        metadata_store_path = finished_path / "jobs_metadata" / task_metadata_files.parent.name
        metadata_store_path.mkdir(parents=True, exist_ok=True)
        os.system(f"cp {task_metadata_files} {metadata_store_path}")

    gcn_and_cnn_count = {"GCN": 0, "CNN": 0}
    for task_metadata_files in list(job_work_path.glob("*/metadata_cnn_gcn_counts.json")):
        task_metadata = json.load(open(task_metadata_files))
        for model in ["GCN", "CNN"]:
            gcn_and_cnn_count[model] += task_metadata[model]
    json.dump(gcn_and_cnn_count, open(finished_path / "metadata_cnn_gcn_counts.json", "w"), indent=4)

    os.system(f"cp {job_work_path}/metadata* {finished_path}")
    print(f"FINISHED {job_conf.project_name}/{job_conf.timestamp}")


# main_pipeline() simply executes functions implemented above step by step
def main_pipeline(fsc: FolderStructureConfig, project_name: str, input_paths: list, target_db_name: str,
                  parallel_jobs: int, threads: int):
    job_work_path = prepare_job(fsc, project_name, input_paths, target_db_name, parallel_jobs)

    split_job_into_tasks(job_work_path)
    run_parallel_pipelines(fsc, job_work_path, threads=threads)
    merge_completed_task_results(fsc, job_work_path)
    shutil.rmtree(job_work_path)


def main():
    args = parse_args()

    data_root = pathlib.Path(args.data_root)
    if data_root.is_dir():
        fsc = FolderStructureConfig(data_root)
    else:
        fsc = load_folder_structure_config(data_root)

    project_name = args.project_name

    if args.target_db_name is None:
        try:
            _ = find_target_database(fsc, project_name)
            target_db_name = project_name
        except AssertionError:
            print(f"No {project_name} target database found. Using {DEFAULT_NAME} instead.")
            target_db_name = DEFAULT_NAME
    else:
        # todo add possibility to use specific target_database path and timestamp instead of name only
        target_db_name = args.target_db_name
        _ = find_target_database(fsc, target_db_name)

    input_paths = parse_input_paths(args.input, project_name, fsc.QUERY_PATH)

    main_pipeline(fsc, project_name, input_paths, target_db_name, args.n_parallel_tasks, args.threads)


if __name__ == '__main__':
    main()
