import dataclasses
import json
from typing import List

from meta_deepFRI.config.runtime_config import RuntimeConfig


@dataclasses.dataclass
class JobConfig(RuntimeConfig):
    project_name: str = "None"
    target_db: str = "None"
    target_db_name: str = "None"
    timestamp: str = "None"
    n_parallel_tasks: int = 1
    query_files: List[str] = "None"

    def save(self, path: str) -> None:
        json.dump(dataclasses.asdict(self), open(path, "w"), indent=4)


def load_job_config(job_config_path: str) -> JobConfig:
    return JobConfig(**json.load(open(job_config_path)))
