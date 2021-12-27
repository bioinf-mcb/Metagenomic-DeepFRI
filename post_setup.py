import requests

from CONFIG import *


def main():
    DATA_ROOT.mkdir(exist_ok=True, parents=True)
    STRUCTURE_FILES_PATH.mkdir(exist_ok=True, parents=True)
    SEQ_CMAP_DATASET_PATH.mkdir(exist_ok=True, parents=True)

    MMSEQS_DATABASES_PATH.mkdir(exist_ok=True, parents=True)
    WORK_PATH.mkdir(exist_ok=True, parents=True)
    QUERY_PATH.mkdir(exist_ok=True, parents=True)
    FINISHED_PATH.mkdir(exist_ok=True, parents=True)

    if not (DATA_ROOT / "trained_models/model_config.json").exists():
        from utils.utils import run_command
        print(f"No model config.json file found in {DATA_ROOT / 'trained_models'}.")
        if not pathlib.Path("newest_trained_models.tar.gz").exists():
            print(" Downloading models, approx 800MB")
            r = requests.get(TRAINED_MODELS_DOWNLOAD_URL, allow_redirects=True)
            open('newest_trained_models.tar.gz', 'wb').write(r.content)
        print(f"unloading models into {DATA_ROOT} directory")
        run_command(f"tar xvzf newest_trained_models.tar.gz -C {DATA_ROOT}")


if __name__ == "__main__":
    main()
