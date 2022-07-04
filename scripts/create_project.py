import argparse
import json
import pathlib

from meta_deepFRI.config.folder_structure import FolderStructureConfig, load_folder_structure_config
from meta_deepFRI.config.names import PROJECT_CONFIG, DEFAULT_NAME, TARGET_DB_CONFIG, SEQUENCES, ATOMS
from meta_deepFRI.config.runtime_config import load_runtime_config, RuntimeConfig


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(description="This script contains all the logic to run DeepFri's GCN or CNN experiments.")

    parser.add_argument("--data_root", required=False, default="config/data_root_path.json",
                        help="Path to json file containing DATA_ROOT or path to folder that will be used as DATA_ROOT")

    parser.add_argument("--runtime_config", required=False, default="default_config/runtime_config.json",
                        help="Task name")

    parser.add_argument("--project_name", required=False, default=DEFAULT_NAME,
                        help="Task name")
    # yapf: enable
    return parser.parse_args()


def create_project(fsc: FolderStructureConfig, rc: RuntimeConfig, project_name: str) -> None:
    print("Creating new project:", project_name)
    (fsc.QUERY_PATH / project_name).mkdir(exist_ok=True, parents=True)
    (fsc.WORK_PATH / project_name).mkdir(exist_ok=True, parents=True)
    (fsc.STRUCTURE_FILES_PATH / project_name).mkdir(exist_ok=True, parents=True)
    (fsc.FINISHED_PATH / project_name).mkdir(exist_ok=True, parents=True)

    (fsc.MMSEQS_DATABASES_PATH / project_name).mkdir(exist_ok=True, parents=True)
    (fsc.SEQ_ATOMS_DATASET_PATH / project_name).mkdir(exist_ok=True, parents=True)
    (fsc.SEQ_ATOMS_DATASET_PATH / project_name / SEQUENCES).mkdir(exist_ok=True)
    (fsc.SEQ_ATOMS_DATASET_PATH / project_name / ATOMS).mkdir(exist_ok=True)

    project_config_path = fsc.WORK_PATH / project_name / PROJECT_CONFIG
    if not project_config_path.exists():
        rc.save(project_config_path)
    else:
        print(f"Project config already exists. {project_config_path}")

    target_db_config_path = fsc.MMSEQS_DATABASES_PATH / project_name / TARGET_DB_CONFIG
    if not target_db_config_path.exists():
        json.dump(rc.MAX_TARGET_CHAIN_LENGTH, open(target_db_config_path, "w"), indent=4)
    else:
        print(f"Target db config already exists. {target_db_config_path}")


def main():
    args = parse_args()

    data_root = pathlib.Path(args.data_root)
    if data_root.is_dir():
        fsc = FolderStructureConfig(data_root)
    else:
        fsc = load_folder_structure_config(data_root)

    rc = load_runtime_config(args.runtime_config)
    project_name = args.project_name

    create_project(fsc, rc, project_name)


if __name__ == '__main__':
    main()
