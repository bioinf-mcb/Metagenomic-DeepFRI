import json
import os
import pathlib
from typing import List, Union, Any
import shlex

import requests
import shutil
import subprocess
import sys

ENV_PATHS = set()


def add_path_to_env(path):
    ENV_PATHS.add(path)


def run_command(command, timeout=-1):
    if isinstance(command, str):
        command = shlex.split(command, ' ')

    my_env = os.environ.copy()
    my_env["PATH"] += os.pathsep + os.pathsep.join(ENV_PATHS)

    try:
        if timeout > 0:
            completed_process = subprocess.run(command,
                                               capture_output=True,
                                               env=my_env,
                                               timeout=timeout,
                                               check=True,
                                               universal_newlines=True)
        else:
            completed_process = subprocess.run(command,
                                               capture_output=True,
                                               env=my_env,
                                               check=True,
                                               universal_newlines=True)

    except subprocess.TimeoutExpired:
        raise TimeoutError(f"command {' '.join(command)} timed out") from None

    except subprocess.CalledProcessError as err:
        raise RuntimeError(
            f"Command '{' '.join(command)}' failed with exit code {err.returncode}: {err.stderr}") from err

    return completed_process.stdout


def download_file(url, path):
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def merge_files_binary(file_paths: list, output_path: pathlib.Path) -> None:
    with open(output_path, 'wb') as writer:
        for input_file in file_paths:
            with open(input_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)


def search_files_in_paths(paths: list, pattern: str):
    files = []
    for path in paths:
        if not path.exists():
            print(f"Unable to locate {path}.")
            continue
        if path.is_dir():
            files.extend(list(path.glob("**/*" + pattern)))
        else:
            if not path.name.endswith(pattern):
                print(f"{path} is not an {pattern} file which is excepted format.")
            else:
                files.append(path)
    return files


def shutdown(message):
    """
    Terminates program execution with a reason.

    Args:
        message (str): Reason for termination.
    """
    sys.exit(message)


def load_deepfri_config(filepath_model_config_json: str) -> dict:
    """Loads model_config.json with paths to different models.

    Args:
        filepath_model_config_json (str): path to a file within trained_models folder.
        Distributed with original DeepFRI repo.

    Returns:
        dict: a dict of different models and paths to their weights.
    """
    json_filepath = pathlib.Path(filepath_model_config_json)
    if not json_filepath.exists():
        raise FileNotFoundError(f"Config file not found at {json_filepath}.")

    # load and replace local paths to files with absolute paths
    with open(filepath_model_config_json, "r") as json_file:
        json_string = json_file.read().replace("./trained_models", f"{json_filepath.parent}")
    deepfri_config = json.loads(json_string)

    return deepfri_config
