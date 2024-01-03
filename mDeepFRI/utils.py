import json
import logging
import re
import shlex
import shutil
import subprocess
import sys
from glob import glob
from pathlib import Path
from typing import Iterable, Literal

import requests

from mDeepFRI import cnn_model_links, gcn_model_links


def run_command(command, timeout=None):
    """
    Runs a command and returns its output.

    Args:
        command (str): Command to run.
        timeout (int): Timeout in seconds.

    Returns:
        str: Command output.
    """
    if isinstance(command, str):
        command = shlex.split(command, ' ')

    try:
        completed_process = subprocess.run(command,
                                           capture_output=True,
                                           timeout=timeout,
                                           check=True,
                                           universal_newlines=True)

    except subprocess.TimeoutExpired:
        raise TimeoutError(f"command {' '.join(command)} timed out") from None

    except subprocess.CalledProcessError as err:
        raise RuntimeError(
            f"Command '{' '.join(command)}' failed with exit code {err.returncode}\n{err.stderr}"
        ) from err

    return completed_process.stdout


def download_file(url, path):
    """
    Downloads a file from url and saves it to path.

    Args:
        url (str): URL to download.
        path (str): Path to save the file.

    Returns:
        None
    """

    with requests.get(url, stream=True) as r:
        # check response and raise error
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as err:
            raise RuntimeError(
                f"Download of {url} failed with error code {err.response.status_code}"
            ) from err

        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def download_model_weights(output_path: str,
                           version: Literal["1.0", "1.1"] = "1.1") -> None:
    """
    Downloads model weights and configs from the internet.

    Args:
        output_path (str): Path to save the files.

    Returns:
        None
    """

    output_path = Path(output_path)
    try:
        output_path.mkdir()
    except FileExistsError:
        # clean up the folder to avoid version conflicts
        shutil.rmtree(output_path)
        output_path.mkdir()

    for mode in gcn_model_links[version]:
        logging.debug("Downloading GCN %s models...", mode.upper())
        for url in gcn_model_links[version][mode].values():
            download_file(url, output_path / url.split("/")[-1])

    for mode in cnn_model_links:
        logging.debug("Downloading CNN %s models...", mode.upper())
        for url in cnn_model_links[mode].values():
            # version 1.1 does not perdict EC number
            if version == "1.1":
                if mode == "ec":
                    continue
            download_file(url, output_path / url.split("/")[-1])


## TODO: automatical generation of a config JSON
def generate_config_json(weights_path: str, version: Literal["1.0",
                                                             "1.1"]) -> None:
    """
    Generates a config json file.

    Args:
        weights_path (str): Path to the weights folder.
        version (str): Version of the model.

    Returns:
        None
    """

    weights_path = Path(weights_path)
    config = {
        "gcn": {
            "bp": None,
            "cc": None,
            "mf": None,
            "ec": None
        },
        "cnn": {
            "bp": None,
            "cc": None,
            "mf": None,
            "ec": None
        },
        "version": None
    }

    models = list(weights_path.glob("*.onnx"))
    possible_modes = "|".join(list(config["cnn"].keys()))
    for model in models:
        mode = re.search(possible_modes, model.name).group(0)
        if "CNN" in model.name:
            config["cnn"][mode] = str(model)
        elif "GraphConv" in model.name:
            config["gcn"][mode] = str(model)

    config["version"] = version
    # version 1.1 doesn't predict EC
    if version == "1.1":
        del config["cnn"]["ec"]
        del config["gcn"]["ec"]

    config_name = "model_config.json"
    with open(weights_path / config_name, "w") as f:
        json.dump(config, f, indent=4, sort_keys=True)


# def merge_files_binary(file_paths: list, output_path: pathlib.Path) -> None:
#     """
#     Merges files in binary format.

#     Args:
#         file_paths (list): List of paths to merge.
#         output_path (str): Path to save the merged file.

#     Returns:
#         None
#     """

#     with open(output_path, 'wb') as writer:
#         for input_file in file_paths:
#             with open(input_file, 'rb') as reader:
#                 shutil.copyfileobj(reader, writer)

# def search_files_in_paths(paths: list, pattern: str):
#     """
#     Searches for files in paths.

#     Args:
#         paths (list): List of paths to search.
#         pattern (str): Pattern to search for.

#     Returns:
#         list: List of files found.
#     """

#     files = []
#     for path in paths:
#         if not path.exists():
#             print(f"Unable to locate {path}.")
#             continue
#         if path.is_dir():
#             files.extend(list(path.glob("**/*" + pattern)))
#         else:
#             if not path.name.endswith(pattern):
#                 print(
#                     f"{path} is not an {pattern} file which is excepted format."
#                 )
#             else:
#                 files.append(path)
#     return files


def shutdown(message):
    """
    Terminates program execution with a reason.

    Args:
        message (str): Reason for termination.
    """
    sys.exit(message)


def remove_intermediate_files(temporary_files: Iterable):
    """
    Removes temporary files.

    Args:
        temporary_files (Iterable): List of temporary files.

    Returns:
        None
    """

    for file in temporary_files:
        extensions = glob(str(file) + "*")
        for ext in extensions:
            logging.info(f"Removing temporary file {ext}.")
            Path(ext).unlink()


def load_deepfri_config(weights: str) -> dict:
    """
    Check if DeepFRI weights are valid and load config.
    Args:
        weights: Path to DeepFRI weights.

    Returns:
        pathlib.Path: Path to DeepFRI config.
    """

    weights = Path(weights)
    assert weights.exists(), f"DeepFRI weights not found at {weights}"
    assert weights.is_dir(
    ), "DeepFRI weights should be a directory, not a file."

    config_path = weights / "model_config.json"
    assert config_path.exists(
    ), "DeepFRI weights are missing model_config.json"

    with open(config_path, "r", encoding="utf-8") as f:
        models_config = json.load(f)

    for net in ["cnn", "gcn"]:
        for model_type, model_path in models_config[net].items():
            model_name = weights / (Path(model_path).name)
            config_name = weights / (Path(model_path).stem +
                                     "_model_params.json")
            assert model_name.exists(
            ), f"DeepFRI weights are missing {model_type} model at {model_name}"
            assert config_name.exists(
            ), f"DeepFRI weights are missing {model_type} model config at {config_name}"
            # correct config
            models_config[net][model_type] = str(model_name.absolute())

    return models_config
