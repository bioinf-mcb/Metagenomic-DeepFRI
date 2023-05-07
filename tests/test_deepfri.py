import pathlib
import subprocess

import pytest

INPUT_STRUCTURES = 'tests/data'
OUTPUT_DATABASE = 'test_database'
QUERY_FILE = 'tests/data/small_query.faa'
RESULTS = 'test_results'


def execute_command(command):
    process = subprocess.run(command.split(),
                             check=False,
                             capture_output=True,
                             text=True)
    if process.returncode != 0:
        raise RuntimeError(
            f'Command {command} failed with exit code {process.returncode}\n{process.stderr}'
        )
    return process.stdout, process.stderr


@pytest.fixture
def deepfri_database():
    stdout, stderr = execute_command(
        f'mDeepFRI build-db -i {INPUT_STRUCTURES} -o {OUTPUT_DATABASE}')
    print(stdout)


def test_database(deepfri_database):
    return


@pytest.fixture
def deepri_results(deepfri_database):
    stdout, stderr = execute_command(
        f"mDeepFRI -i {QUERY_FILE} -db {OUTPUT_DATABASE} -o {RESULTS}")
    print(stdout)


def test_results(deepri_results):
    results_path = pathlib.Path(RESULTS)
    assert results_path.exists()
    assert results_path.is_dir()
    assert len(list(results_path.glob('*.tsv'))) == 8


def test_clean(deepri_results):
    execute_command(f'rm -r {OUTPUT_DATABASE} {RESULTS}')
