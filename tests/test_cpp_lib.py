from typing import Callable

from mDeepFRI.CPP_lib.libAtomDistanceIO import (initialize,
                                                load_aligned_contact_map,
                                                load_contact_map, save_atoms)


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
