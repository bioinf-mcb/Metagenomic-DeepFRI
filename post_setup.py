import json
import os
import pathlib
import requests
import shutil

from CONFIG.FOLDER_STRUCTURE import DATA_ROOT, STRUCTURE_FILES_PATH, QUERY_PATH, DEFAULT_NAME, FINISHED_PATH, WORK_PATH, \
    SEQ_ATOMS_DATASET_PATH, MMSEQS_DATABASES_PATH, DEEPFRI_MODEL_WEIGHTS_JSON_FILE, DEEPFRI_TRAINED_MODELS_DOWNLOAD_URL, \
    PROJECT_CONFIG, TARGET_DB_CONFIG

from CONFIG.get_config_dict import target_db_config, runtime_config


def download_file(url, path):
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def main():
    print("Creating folders structure based on CONFIG/FOLDER_STRUCTURE.py")
    DATA_ROOT.mkdir(exist_ok=True, parents=True)
    STRUCTURE_FILES_PATH.mkdir(exist_ok=True, parents=True)
    QUERY_PATH.mkdir(exist_ok=True, parents=True)
    (QUERY_PATH / DEFAULT_NAME).mkdir(exist_ok=True, parents=True)
    FINISHED_PATH.mkdir(exist_ok=True, parents=True)

    WORK_PATH.mkdir(exist_ok=True, parents=True)
    (WORK_PATH / DEFAULT_NAME).mkdir(exist_ok=True, parents=True)
    if not (WORK_PATH / DEFAULT_NAME / PROJECT_CONFIG).exists():
        config = runtime_config()
        json.dump(config, open(WORK_PATH / DEFAULT_NAME / PROJECT_CONFIG, "w"), indent=4)

    SEQ_ATOMS_DATASET_PATH.mkdir(exist_ok=True, parents=True)
    MMSEQS_DATABASES_PATH.mkdir(exist_ok=True, parents=True)
    (MMSEQS_DATABASES_PATH / DEFAULT_NAME).mkdir(exist_ok=True, parents=True)
    if not (WORK_PATH / DEFAULT_NAME / TARGET_DB_CONFIG).exists():
        config = target_db_config()
        json.dump(config, open(WORK_PATH / DEFAULT_NAME / TARGET_DB_CONFIG, "w"), indent=4)

    if not DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists():
        print(f"No model config.json file found in {DATA_ROOT / 'trained_models'}.")

        if not pathlib.Path("newest_trained_models.tar.gz").exists():
            print("Downloading model weights, approx 800MB")
            download_file(DEEPFRI_TRAINED_MODELS_DOWNLOAD_URL, 'newest_trained_models.tar.gz')

        print(f"unloading models into {DATA_ROOT / 'trained_models'} directory")
        os.system(f"tar xvzf newest_trained_models.tar.gz -C {DATA_ROOT}")
    else:
        print("Found model weights")
    print("All good and ready to go!")


if __name__ == "__main__":
    main()
