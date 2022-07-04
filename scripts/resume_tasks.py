import argparse
import json
import pathlib
import shutil

from config.job_config import load_job_config
from meta_deepFRI.config.folder_structure import FolderStructureConfig, load_folder_structure_config
from meta_deepFRI.config.names import TASK_CONFIG, JOB_CONFIG
from main_pipeline import run_parallel_pipelines, merge_completed_task_results, split_job_into_tasks


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(description="Read structure files from folders --input to extract sequence and atom positions. "
                                                 "Create and index new --output MMSEQS2 database")

    parser.add_argument("-r", "--data_root", required=False, default="default_config/data_root_path.json",
                        help="Path to json file containing DATA_ROOT or path to folder that will be used as DATA_ROOT")

    return parser.parse_args()
    # yapf: enable


def main():
    args = parse_args()

    data_root = pathlib.Path(args.data_root)
    if data_root.is_dir():
        fsc = FolderStructureConfig(data_root)
    else:
        fsc = load_folder_structure_config(data_root)

    job_config_paths = list(fsc.WORK_PATH.glob(f"**/{JOB_CONFIG}"))
    for job_config_path in job_config_paths:
        job_work_path = job_config_path.parent
        job_config = load_job_config(job_config_path)

        # check if all job directories were created.
        tasks_configs = list(job_work_path.glob(f"**/{TASK_CONFIG}"))
        if len(tasks_configs) != job_config.n_parallel_tasks:
            split_job_into_tasks(job_work_path)

        run_parallel_pipelines(fsc, job_work_path)
        merge_completed_task_results(fsc, job_work_path)
        shutil.rmtree(job_work_path)


if __name__ == '__main__':
    main()
