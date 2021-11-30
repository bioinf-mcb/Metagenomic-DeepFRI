import os
import subprocess

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


def create_chunks(lst, n):
    return [lst[i::n] for i in range(n)]
