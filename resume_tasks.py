import json
import shutil

from CONFIG.FOLDER_STRUCTURE import WORK_PATH, TASK_CONFIG, JOB_CONFIG

from main_pipeline import run_parallel_pipelines, merge_completed_task_results, split_task_into_jobs


def main():
    task_configs = list(WORK_PATH.glob(f"**/{TASK_CONFIG}"))
    for task_config in task_configs:
        task_work_path = task_config.parent
        config = json.load(open(task_config))

        # check if all job directories were created.
        jobs_config = list(task_work_path.glob(f"**/{JOB_CONFIG}"))
        if len(jobs_config) != config["n_parallel_jobs"]:

            split_task_into_jobs(task_work_path)

        run_parallel_pipelines(task_work_path)
        merge_completed_task_results(task_work_path)
        shutil.rmtree(task_work_path)


if __name__ == '__main__':
    main()
