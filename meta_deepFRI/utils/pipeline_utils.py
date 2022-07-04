import json
import pathlib

from meta_deepFRI.config.names import DEFAULT_NAME, TARGET_MMSEQS_DB_NAME
from meta_deepFRI.config.folder_structure import FolderStructureConfig


def find_target_database(fsc: FolderStructureConfig, project_name: str) -> pathlib.Path:
    # todo add possibility to use specific target_database path and timestamp instead of name only
    # excepted folder structure is MMSEQS_DATABASES_PATH / project_name / timestamp / TARGET_MMSEQS_DB_NAME
    assert (fsc.MMSEQS_DATABASES_PATH / project_name).exists(),         \
        f"No target databases folder found {fsc.MMSEQS_DATABASES_PATH / project_name}. " \
        f"Create one using update_mmseqs_database.py {'' if project_name == DEFAULT_NAME else '--name ' + project_name}"

    target_databases = sorted([x for x in (fsc.MMSEQS_DATABASES_PATH / project_name).iterdir() if x.is_dir()])

    assert len(target_databases) > 0,\
        f"No target databases found {fsc.MMSEQS_DATABASES_PATH / project_name}/timestamp. " \
        f"Create one using update_mmseqs_database.py {'' if project_name == DEFAULT_NAME else '--name ' + project_name}"

    assert (fsc.SEQ_ATOMS_DATASET_PATH / project_name).exists(),\
        f"There are no atom and sequences required by {project_name}. Something went wrong. " \
        f"Please, repeat target database creation update_mmseqs_database.py {'' if project_name == DEFAULT_NAME else '--name ' + project_name}"

    return target_databases[-1] / TARGET_MMSEQS_DB_NAME


def load_deepfri_config(fsc: FolderStructureConfig) -> dict:
    assert fsc.DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists(),\
        f"No DeepFri model weights json file found at {fsc.DEEPFRI_MODEL_WEIGHTS_JSON_FILE} " \
        f"Please run post_setup.py script to download and unzip model weights and json file"

    # load and replace local paths to files with absolute paths
    with open(fsc.DEEPFRI_MODEL_WEIGHTS_JSON_FILE, "r") as json_file:
        json_string = json_file.read().replace("./trained_models", f"{fsc.DEEPFRI_MODEL_WEIGHTS_JSON_FILE.parent}")
    deepfri_config = json.loads(json_string)

    return deepfri_config
