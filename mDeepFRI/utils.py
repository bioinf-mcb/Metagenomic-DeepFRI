import logging
import pathlib
import shlex
import shutil
import subprocess
import sys
from glob import glob
from pathlib import Path
from typing import Iterable

import requests

from mDeepFRI import cnn_model_links, config_links, gcn_model_links


def run_command(command, timeout=None):
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
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def download_model_weights(output_path: pathlib.Path):

    model_links = list(cnn_model_links.values()) + list(
        gcn_model_links.values())
    total_len = len(model_links)
    for i, model_link in enumerate(model_links):
        download_file(model_link, output_path / model_link.split("/")[-1])
        logging.debug("Downloading model weights... %s/%s", i + 1, total_len)

    total_len = len(config_links)
    for i, config in enumerate(config_links):
        download_file(config, output_path / config.split("/")[-1])
        logging.debug("Downloading model configs... %s/%s", i + 1, total_len)


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
                print(
                    f"{path} is not an {pattern} file which is excepted format."
                )
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


def remove_temporary(temporary_files: Iterable):
    for file in temporary_files:
        extensions = glob(str(file) + "*")
        for ext in extensions:
            logging.info(f"Removing temporary file {ext}.")
            Path(ext).unlink()
