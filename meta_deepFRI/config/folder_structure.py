import json
import pathlib


class FolderStructureConfig:
    """
    DATA_ROOT - as name suggests, root directory for all the pipeline files.
    STRUCTURE_FILES_PATH - path in which one should put their protein structure files. In any structure.
    QUERY_PATH - root for all query .faa files. They should be placed in subdirectories named after project_name
    WORK_PATH - contains all currently running or interrupted tasks. Each project folder also contain runtime parameters in project_config.json
    FINISHED_PATH - root for all finished pipeline results. Will contain folders named after project_name
    """
    def __init__(self, data_root: pathlib.Path = pathlib.Path("/meta_deepfri_data")):
        self.DATA_ROOT = data_root

        # pipeline folders
        self.STRUCTURE_FILES_PATH = data_root / "structure_files"
        self.QUERY_PATH = data_root / "query"
        self.WORK_PATH = data_root / "workspace"
        self.FINISHED_PATH = data_root / "finished"

        # database folders
        self.SEQ_ATOMS_DATASET_PATH = data_root / "seq_atoms_dataset"
        self.MMSEQS_DATABASES_PATH = data_root / "mmseqs_db"

        # deepfri config file path
        self.DEEPFRI_MODEL_WEIGHTS_JSON_FILE = data_root / "trained_models/model_config.json"


def load_default_folder_structure_config() -> FolderStructureConfig:
    return FolderStructureConfig()


def load_folder_structure_config(folder_structure_config_path: pathlib.Path) -> FolderStructureConfig:
    return FolderStructureConfig(pathlib.Path(json.load(open(folder_structure_config_path))))
