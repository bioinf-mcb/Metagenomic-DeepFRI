import filecmp
import pickle
import tempfile
from typing import Callable

import numpy as np
import pytest

from mDeepFRI.CPP_lib.atoms_io import load_atoms_file, save_atoms_file
from mDeepFRI.CPP_lib.libAtomDistanceIO import (initialize,
                                                load_aligned_contact_map,
                                                load_contact_map, save_atoms)
from mDeepFRI.CPP_lib.parsers import parse_pdb


@pytest.fixture()
def reference_binary():
    return "tests/data/1S3P-A.bin"


@pytest.fixture()
def struct_params():
    with open("tests/data/structures/1S3P-A.pdb", "rb") as f:
        _, positions, groups = parse_pdb(f.read().split(b"\n"))
    _, group_index = np.unique(groups, return_index=True)
    group_index = np.sort(group_index)
    chain_len = len(positions) / 3
    group_indexes = list(np.append(group_index, chain_len).astype(np.int32))
    return (positions, group_indexes, len(group_indexes))


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


def test_save_atoms(reference_binary, struct_params):
    with tempfile.TemporaryDirectory() as tmpdirname:
        save_atoms_file(struct_params[0], struct_params[1],
                        tmpdirname + "/1S3P-A.bin")
        assert filecmp.cmp(tmpdirname + "/1S3P-A.bin", reference_binary)


def test_load_atoms(reference_binary, struct_params):
    positions_cy, group_indexes_cy, chain_len_cy = load_atoms_file(
        reference_binary)
    assert np.array_equal(positions_cy, np.array(struct_params[0]))
    assert np.array_equal(group_indexes_cy, np.array(struct_params[1]))
    assert chain_len_cy == struct_params[2]


def test_load_aligned_contact_map():
    pass
