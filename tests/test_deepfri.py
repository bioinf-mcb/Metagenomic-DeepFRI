import subprocess

import pytest

TEST_INPUT_DATA = 'tests/data'
TEST_OUTPUT_DATA = 'test_database'


def execute_command(command):
    process = subprocess.run(command.split(), check=False, capture_output=True, text=True)
    if process.returncode != 0:
        raise RuntimeError(f'Command {command} failed with exit code {process.returncode}\n{process.stderr}')
    return process.stdout, process.stderr


@pytest.fixture
def built_database():
    stdout, stderr = execute_command(f'deepfri_db_build -i {TEST_INPUT_DATA} -o {TEST_OUTPUT_DATA}')
    print(stdout)
