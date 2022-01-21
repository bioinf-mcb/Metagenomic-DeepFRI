import os
import pathlib
import requests
import shutil
import subprocess
import time

ENV_PATHS = set()


def add_path_to_env(path):
    ENV_PATHS.add(path)


def run_command(command, timeout=-1):
    if type(command) == str:
        command = str.split(command, ' ')

    my_env = os.environ.copy()
    my_env["PATH"] += ":"+str.join(":", ENV_PATHS)

    try:
        if timeout > 0:
            completed_process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env, timeout=timeout)
        else:
            completed_process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    except subprocess.TimeoutExpired:
        raise TimeoutError(f"command {' '.join(command)} timeout")

    if completed_process.stderr != b'':
        error_info = completed_process.stderr.decode()
        raise RuntimeError(f"during execution: {' '.join(command)} exception occurred\n{error_info}")
    else:
        return completed_process.stdout.decode('utf-8')


def download_file(url, path):
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def create_chunks(lst, n):
    return [lst[i::n] for i in range(n)]


def split_dict(d, n):
    keys = list(d.keys())
    for i in range(0, len(keys), n):
        yield {k: d[k] for k in keys[i: i + n]}


def create_unix_time_folder(parent_path):
    parent_path = pathlib.Path(parent_path)
    start = str(time.time())
    path = (parent_path / start)
    while path.exists():
        time.sleep(1)
        start = str(time.time())
        path = (parent_path / start)
    path.mkdir(parents=True)
    return path


def merge_files_binary(file_paths: list, output_path: pathlib.Path):
    with open(output_path, 'wb') as writer:
        for input_file in file_paths:
            with open(input_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)
