import json
import os
import pathlib

import requests
from tqdm.auto import tqdm

from CONFIG.FOLDER_STRUCTURE import DATA_ROOT, STRUCTURE_FILES_PATH, QUERY_PATH, DEFAULT_NAME, FINISHED_PATH, WORK_PATH, \
    SEQ_ATOMS_DATASET_PATH, MMSEQS_DATABASES_PATH, DEEPFRI_MODEL_WEIGHTS_JSON_FILE, DEEPFRI_TRAINED_MODELS_DOWNLOAD_URL, \
    PROJECT_CONFIG, TARGET_DB_CONFIG

from CONFIG.get_config_dict import target_db_config, runtime_config


def download_file(url, destination, chunk_size=1024):
    with requests.get(url, stream=True) as r:
        total_length = int(r.headers.get("Content-Length"))
        with open(destination, "wb") as f:
            with tqdm(unit="B", unit_scale=True, unit_divisor=1024, total=total_length) as progress_bar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    datasize = f.write(chunk)
                    progress_bar.update(datasize)


def create_folder_structure():
    print("Creating folders structure based on CONFIG/FOLDER_STRUCTURE.py")
    DATA_ROOT.mkdir(exist_ok=True, parents=True)
    STRUCTURE_FILES_PATH.mkdir(exist_ok=True, parents=True)
    QUERY_PATH.mkdir(exist_ok=True, parents=True)
    (QUERY_PATH / DEFAULT_NAME).mkdir(exist_ok=True, parents=True)
    FINISHED_PATH.mkdir(exist_ok=True, parents=True)

    WORK_PATH.mkdir(exist_ok=True, parents=True)
    (WORK_PATH / DEFAULT_NAME).mkdir(exist_ok=True, parents=True)
    # todo check if creation of default config files are actually necessary
    if not (WORK_PATH / DEFAULT_NAME / PROJECT_CONFIG).exists():
        config = runtime_config()
        json.dump(config, open(WORK_PATH / DEFAULT_NAME / PROJECT_CONFIG, "w"), indent=4)

    SEQ_ATOMS_DATASET_PATH.mkdir(exist_ok=True, parents=True)
    MMSEQS_DATABASES_PATH.mkdir(exist_ok=True, parents=True)

    (MMSEQS_DATABASES_PATH / DEFAULT_NAME).mkdir(exist_ok=True, parents=True)
    # todo check if creation of default config files are actually necessary
    if not (WORK_PATH / DEFAULT_NAME / TARGET_DB_CONFIG).exists():
        config = target_db_config()
        json.dump(config, open(WORK_PATH / DEFAULT_NAME / TARGET_DB_CONFIG, "w"), indent=4)


def download_and_extract_deepfri_weights():
    if not DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists():
        print(f"No model config.json file found in {DATA_ROOT / 'trained_models'}.")

        compressed_model_weights_path = pathlib.Path("newest_trained_models.tar.gz")
        if not compressed_model_weights_path.exists():
            print("Downloading model weights")
            download_file(DEEPFRI_TRAINED_MODELS_DOWNLOAD_URL, compressed_model_weights_path)
        print(f"unloading models into {DATA_ROOT / 'trained_models'} directory")
        os.system(f"tar xvzf {compressed_model_weights_path} -C {DATA_ROOT}")
        compressed_model_weights_path.unlink()
    else:
        print("Found model weights")


def main():
    create_folder_structure()
    download_and_extract_deepfri_weights()
    print("Ready to go!")


if __name__ == "__main__":
    main()
