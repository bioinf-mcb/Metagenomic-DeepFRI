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

    assert (sequence[:4] == ['MET', 'MET', 'MET', 'MET']).all()


def test_parse_pdb_positions(pdb_file):
    _, positions, _ = parse_pdb(pdb_file)
    expected = np.array([[-47.025, 49.762, -23.074]], dtype=np.float32)
    assert np.array_equal(positions[:1], expected)


def test_parse_pdb_groups(pdb_file):
    _, _, groups = parse_pdb(pdb_file)
    assert groups[:1] == ['A   1']


def test_parse_pdb_lengths(pdb_file):
    sequence, positions, groups = parse_pdb(pdb_file)
    assert len(sequence) > 0
    assert len(sequence) == len(positions) == len(groups)


@pytest.fixture
def mmcif_file():
    test_file = pathlib.Path(__file__).parent / "data" / "structures" / "6a0j.cif.gz"
    return gzip.open(test_file, "rt")


def test_parse_mmcif_sequence(mmcif_file):
    sequence, _, _ = parse_mmcif(mmcif_file)
    assert all(sequence[:4] == ['VAL', 'VAL', 'VAL', 'VAL'])


def test_parse_mmcif_positions(mmcif_file):
    _, positions, _ = parse_mmcif(mmcif_file)
    expected = np.array([[-6.123, 31.928, 92.243]], dtype=np.float32)
    assert np.array_equal(positions[:1], expected)


def test_parse_mmcif_groups(mmcif_file):
    _, _, groups = parse_mmcif(mmcif_file)
    assert groups[:1] == ['A22']


def test_parse_mmcif_lengths(mmcif_file):
    sequence, positions, groups = parse_mmcif(mmcif_file)
    assert len(sequence) > 0
    assert len(sequence) == len(positions) == len(groups)
