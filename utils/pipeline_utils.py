import json

from CONFIG.FOLDER_STRUCTURE import MMSEQS_DATABASES_PATH, TARGET_MMSEQS_DB_NAME, SEQ_ATOMS_DATASET_PATH, DEFAULT_NAME, \
    DEEPFRI_MODEL_WEIGHTS_JSON_FILE


def find_target_database(target_db_name):
    # todo add possibility to use specific target_database path and timestamp instead of name only
    # excepted folder structure is MMSEQS_DATABASES_PATH / target_db_name / timestamp / TARGET_MMSEQS_DB_NAME
    assert (MMSEQS_DATABASES_PATH / target_db_name).exists(),         \
        f"No target databases folder found {MMSEQS_DATABASES_PATH / target_db_name}. " \
        f"Create one using update_mmseqs_database.py {'' if target_db_name == DEFAULT_NAME else '--name ' + target_db_name}"

    target_databases = sorted([x for x in (MMSEQS_DATABASES_PATH / target_db_name).iterdir() if x.is_dir()])

    assert len(target_databases) > 0,\
        f"No target databases found {MMSEQS_DATABASES_PATH / target_db_name}/timestamp. " \
        f"Create one using update_mmseqs_database.py {'' if target_db_name == DEFAULT_NAME else '--name ' + target_db_name}"

    assert (SEQ_ATOMS_DATASET_PATH / target_db_name).exists(),\
        f"There are no atom and sequences required by {target_db_name}. Something went wrong. " \
        f"Please, repeat target database creation update_mmseqs_database.py {'' if target_db_name == DEFAULT_NAME else '--name ' + target_db_name}"

    return target_databases[-1] / TARGET_MMSEQS_DB_NAME


def load_deepfri_config():
    assert DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists(),\
        f"No DeepFri model weights json file found at {DEEPFRI_MODEL_WEIGHTS_JSON_FILE} " \
        f"Please run post_setup.py script to download and unzip model weights and json file"

    # load and replace local paths to files with absolute paths
    with open(DEEPFRI_MODEL_WEIGHTS_JSON_FILE, "r") as json_file:
        json_string = json_file.read().replace("./trained_models", f"{DEEPFRI_MODEL_WEIGHTS_JSON_FILE.parent}")
    deepfri_config = json.loads(json_string)

    return deepfri_config

