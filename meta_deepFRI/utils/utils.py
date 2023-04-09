import os
import pathlib
from typing import List, Union, Any
import shlex

import requests
import shutil
import subprocess
import time
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
            completed_process = subprocess.run(
                command, capture_output=True, env=my_env, timeout=timeout, check=True, universal_newlines=True)
        else:
            completed_process = subprocess.run(
                command, capture_output=True, env=my_env, check=True, universal_newlines=True)

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


def chunks(lst: List[Any], n: int) -> List[List[Any]]:
    if n == 1:
        return [lst]
    output = []
    for i in range(n):
        output.append(lst[i::n])
    return output


def create_unix_timestamp_folder(parent_path: pathlib.Path) -> pathlib.Path:
    start = str(time.time())
    path = (parent_path / start)
    while path.exists():
        time.sleep(1)
        start = str(time.time())
        path = (parent_path / start)
    path.mkdir(parents=True)
    return path


def merge_files_binary(file_paths: list, output_path: pathlib.Path) -> None:
    with open(output_path, 'wb') as writer:
        for input_file in file_paths:
            with open(input_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)


def parse_input_paths(input_list: Union[None, List[str]], project_name: str,
                      parent_directory: pathlib.Path) -> List[pathlib.Path]:
    if input_list is None:
        input_paths = [pathlib.Path(parent_directory / project_name)]
    else:
        input_paths = []
        for input_path in [pathlib.Path(x) for x in input_list]:
            if input_path.is_absolute():
                input_paths.append(input_path)
            else:
                input_paths.append(parent_directory / input_path)
    return input_paths


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


def query_yes_no(question: str, default: str = "yes") -> str:
    """
    Prompts the user for yes/no input, displaying the specified question text.

    :param str question: The text of the prompt for input.
    :param str default: The default if the user hits <ENTER>. Acceptable values
    are 'yes', 'no', and None.
    :return: 'yes' or 'no'
    """
    valid = {'y': 'yes', 'n': 'no'}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError(f"Invalid default answer: '{default}'")

    choice = default

    while 1:
        user_input = input(question + prompt).lower()
        if not user_input:
            break
        try:
            choice = valid[user_input[0]]
            break
        except (KeyError, IndexError):
            print("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

    return choice


def shutdown(message):
    """
    Terminates program execution with a reason.

    Args:
        message (str): Reason for termination.
    """
    sys.exit(message)
