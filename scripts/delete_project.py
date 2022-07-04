import argparse
import pathlib

from meta_deepFRI.config.folder_structure import FolderStructureConfig, load_folder_structure_config
from meta_deepFRI.utils.utils import query_yes_no


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(
        description="This script contains all the logic to run DeepFri's GCN or CNN experiments.")

    parser.add_argument("--data_root", required=False, default="default_config/data_root_path.json",
                        help="Path to json file containing DATA_ROOT or path to folder that will be used as DATA_ROOT")

    parser.add_argument("--project_name", required=True,
                        help="Task name")
    # yapf: enable
    return parser.parse_args()


def delete_project(fsc: FolderStructureConfig, project_name: str, require_confirmation=True) -> None:
    project_paths = [
        fsc.QUERY_PATH / project_name, fsc.WORK_PATH / project_name, fsc.STRUCTURE_FILES_PATH / project_name,
        fsc.FINISHED_PATH / project_name, fsc.SEQ_ATOMS_DATASET_PATH / project_name,
        fsc.MMSEQS_DATABASES_PATH / project_name
    ]

    if require_confirmation:
        for project_path in project_paths:
            print(f"Path {project_path} contain {len(list(project_path.glob('**/*')))} files and folders")

        if query_yes_no("Delete?", "no") != "yes":
            return

    for project_path in project_paths:
        project_path.unlink()


def main():
    args = parse_args()

    data_root = pathlib.Path(args.data_root)
    if data_root.is_dir():
        fsc = FolderStructureConfig(data_root)
    else:
        fsc = load_folder_structure_config(data_root)

    delete_project(fsc, args.project_name)


if __name__ == '__main__':
    main()
