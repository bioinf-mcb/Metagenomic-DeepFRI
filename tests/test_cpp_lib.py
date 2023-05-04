import filecmp
import pickle
import tempfile
from typing import Callable

import numpy as np
import pytest

from mDeepFRI.CPP_lib import atoms_io
from mDeepFRI.CPP_lib.libAtomDistanceIO import (initialize,
                                                load_aligned_contact_map,
                                                load_contact_map, save_atoms)
from mDeepFRI.CPP_lib.parsers import parse_pdb


@pytest.fixture()
def reference_binary():
    return "tests/data/1S3P-A.bin"


@pytest.fixture()
def expected_contact_map():
    with open("tests/data/1S3P-A.bin.cmap.pkl", "rb") as f:
        expected_contact_map = pickle.load(f)
    return expected_contact_map


def error_wrapper(func: Callable, *args, error_contains: str = None):
    try:
        func(*args)
    except Exception as e:
        if error_contains in str(e):
            return
        else:
            raise e


def test_functions():
    initialize()

    error_wrapper(save_atoms, error_contains="did not match C++ signature")
    error_wrapper(load_contact_map,
                  error_contains="did not match C++ signature")
    error_wrapper(load_aligned_contact_map,
                  error_contains="did not match C++ signature")


def test_load_contact_map(reference_binary, expected_contact_map):
    initialize()
    contact_map_cpp = load_contact_map(reference_binary, 6)
    assert np.array_equal(contact_map_cpp, expected_contact_map)


def test_save_atoms(expected_contact_map):
    with open("tests/data/structures/1S3P-A.pdb") as f:
        _, positions, groups = parse_pdb(f)
    _, group_index = np.unique(groups, return_index=True)
    group_index = np.sort(group_index)
    group_indexes = np.append(group_index,
                              positions.shape[0] / 3).astype(np.int32)

    with tempfile.TemporaryDirectory() as tmpdirname:
        atoms_io.save_atoms_file(positions, group_indexes,
                                 tmpdirname + "/1S3P-A.bin")
        assert filecmp.cmp(tmpdirname + "/1S3P-A.bin", "tests/data/1S3P-A.bin")


def test_load_aligned_contact_map():
    pass
