import gzip
import pathlib

import numpy as np
import pytest

from mDeepFRI.parsers import parse_mmcif, parse_pdb


@pytest.fixture
def pdb_file():
    test_file = pathlib.Path(
        __file__
    ).parent / "data" / "structures" / "AF-A0A2Z5TJB0-F1-model_v4.pdb.gz"
    return gzip.open(test_file, "rb").read().split(b"\n")


def test_parse_pdb_sequence(pdb_file):
    sequence, _, _ = parse_pdb(pdb_file)

    assert sequence[:4] == [b'MET', b'MET', b'MET', b'MET']


def test_parse_pdb_positions(pdb_file):
    _, positions, _ = parse_pdb(pdb_file)
    expected = np.array([-47.025, 49.762, -23.074], dtype=np.float32)
    assert np.array_equal(positions[:3], expected)


def test_parse_pdb_groups(pdb_file):
    _, _, groups = parse_pdb(pdb_file)
    assert groups[:1] == [b'A   1']


def test_parse_pdb_lengths(pdb_file):
    sequence, positions, groups = parse_pdb(pdb_file)
    assert len(sequence) > 0
    assert len(sequence) == len(positions) / 3 == len(groups)


@pytest.fixture
def mmcif_file():
    test_file = pathlib.Path(
        __file__).parent / "data" / "structures" / "6a0j.cif.gz"
    return gzip.open(test_file, "rb").read().split(b"\n")


def test_parse_mmcif_sequence(mmcif_file):
    sequence, _, _ = parse_mmcif(mmcif_file)
    assert sequence[:4] == [b'VAL', b'VAL', b'VAL', b'VAL']


def test_parse_mmcif_positions(mmcif_file):
    _, positions, _ = parse_mmcif(mmcif_file)
    expected = np.array([-6.123, 31.928, 92.243], dtype=np.float32)
    assert np.array_equal(positions[:3], expected)


def test_parse_mmcif_groups(mmcif_file):
    _, _, groups = parse_mmcif(mmcif_file)
    assert groups[:1] == [b'A22']


def test_parse_mmcif_lengths(mmcif_file):
    sequence, positions, groups = parse_mmcif(mmcif_file)
    assert len(sequence) > 0
    assert len(sequence) == len(positions) / 3 == len(groups)
