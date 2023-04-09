import numpy as np
import gzip
import pathlib
import pytest

from meta_deepFRI.structure_files.parsers import parse_pdb, parse_mmcif


@pytest.fixture
def pdb_file():
    test_file = pathlib.Path(__file__).parent / "data" / "structures" / "AF-A0A2Z5TJB0-F1-model_v4.pdb.gz"
    return gzip.open(test_file, "rt")


def test_parse_pdb_sequence(pdb_file):
    sequence, _, _ = parse_pdb(pdb_file)

    assert sequence[:4] == ['MET', 'MET', 'MET', 'MET']


def test_parse_pdb_positions(pdb_file):
    _, positions, _ = parse_pdb(pdb_file)
    expected = np.array([[-47.025, 49.762, -23.074]], dtype=np.float32)
    assert np.array_equal(positions[:1], expected)


def test_parse_pdb_groups(pdb_file):
    _, _, groups = parse_pdb(pdb_file)
    assert groups[:1] == ['A   1']


@pytest.fixture
def mmcif_file():
    test_file = pathlib.Path(__file__).parent / "data" / "structures" / "1gqv-sf.cif.gz"
    return gzip.open(test_file, "rt")
