import shutil

from CONFIG.FOLDER_STRUCTURE import WORK_PATH, TASK_CONFIG

from main import parallel_pipelines, merge_finalized_task_results

if __name__ == '__main__':
    task_work_paths = [task.parent for task in WORK_PATH.glob(f"**/{TASK_CONFIG}")]
    for task_work_path in task_work_paths:
        parallel_pipelines(task_work_path)
        merge_finalized_task_results(task_work_path)
        shutil.rmtree(task_work_path)
