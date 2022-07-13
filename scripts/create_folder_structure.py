import argparse
import os
import pathlib

import requests
from tqdm.auto import tqdm

from meta_deepFRI.config.folder_structure import FolderStructureConfig, load_folder_structure_config
from meta_deepFRI.config.names import DEFAULT_NAME
from meta_deepFRI.config.runtime_config import load_default_runtime_config
from scripts.create_project import create_project


def parse_args():
    # yapf: disable
    parser = argparse.ArgumentParser(
        description="This script contains all the logic to run DeepFri's GCN or CNN experiments.")

    parser.add_argument("--data_root", required=False, default="default_config/data_root_path.json",
                        help="Path to json file containing DATA_ROOT or path to folder that will be used as DATA_ROOT")

    parser.add_argument("--deepfri_weights_url", required=False, default=
                        "https://users.flatironinstitute.org/vgligorijevic/public_www/DeepFRI_data/newest_trained_models.tar.gz")
    # yapf: enable
    return parser.parse_args()


def download_file(url, destination, chunk_size=1024):
    with requests.get(url, stream=True) as r:
        total_length = int(r.headers.get("Content-Length"))
        with open(destination, "wb") as f:
            with tqdm(unit="B", unit_scale=True, unit_divisor=1024, total=total_length) as progress_bar:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    datasize = f.write(chunk)
                    progress_bar.update(datasize)


def create_folder_structure(fsc: FolderStructureConfig):
    print("Creating folders structure")
    try:
        fsc.DATA_ROOT.mkdir(exist_ok=True, parents=True)
    except PermissionError:
        print("Need permissions to create this folder. \n"
              f"sudo mkdir -m 777 {fsc.DATA_ROOT}")
        exit(1)

    fsc.STRUCTURE_FILES_PATH.mkdir(exist_ok=True, parents=True)
    fsc.QUERY_PATH.mkdir(exist_ok=True, parents=True)
    fsc.FINISHED_PATH.mkdir(exist_ok=True, parents=True)

    fsc.WORK_PATH.mkdir(exist_ok=True, parents=True)

    fsc.SEQ_ATOMS_DATASET_PATH.mkdir(exist_ok=True, parents=True)
    fsc.MMSEQS_DATABASES_PATH.mkdir(exist_ok=True, parents=True)


def download_and_extract_deepfri_weights(fsc: FolderStructureConfig, url: str):
    if not fsc.DEEPFRI_MODEL_WEIGHTS_JSON_FILE.exists():
        print(f"No model config.json file found in {fsc.DATA_ROOT / 'trained_models'}.")

        compressed_model_weights_path = pathlib.Path("newest_trained_models.tar.gz")
        if not compressed_model_weights_path.exists():
            print("Downloading model weights")
            download_file(url, compressed_model_weights_path)
        print(f"unloading models into {fsc.DATA_ROOT / 'trained_models'} directory")
        os.system(f"tar xvzf {compressed_model_weights_path} -C {fsc.DATA_ROOT}")
        compressed_model_weights_path.unlink()
    else:
        print("Found model weights")


def main():
    args = parse_args()

    data_root = pathlib.Path(args.data_root)
    if data_root.is_dir():
        fsc = FolderStructureConfig(data_root)
    else:
        fsc = load_folder_structure_config(data_root)
    default_runtime_config = load_default_runtime_config()
    deepfri_weights_url = args.deepfri_weights_url

    create_folder_structure(fsc)
    download_and_extract_deepfri_weights(fsc, deepfri_weights_url)
    create_project(fsc, default_runtime_config, DEFAULT_NAME)
    print("Ready to go!")


if __name__ == "__main__":
    main()
