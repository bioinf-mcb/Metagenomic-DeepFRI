import subprocess


def run_command(command, timeout=-1):
    if type(command) == str:
        command = str.split(command, ' ')

    print(str.join(" ", command))

    try:
        if timeout > 0:
            completed_process = subprocess.run(command, timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            completed_process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except subprocess.TimeoutExpired:
        raise TimeoutError(f"command {' '.join(command)} timeout")

    if completed_process.stderr != b'':
        error_info = completed_process.stderr.decode()
        raise RuntimeError(f"during execution: {' '.join(command)} exception occurred\n{error_info}")
    else:
        return completed_process.stdout.decode('utf-8')


def create_chunks(lst, n):
    return [lst[i::n] for i in range(n)]

