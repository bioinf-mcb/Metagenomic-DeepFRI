from meta_deepFRI.CPP_lib.libAtomDistanceIO import initialize, save_atoms, load_contact_map, load_aligned_contact_map


def error_wrapper(func, *args, error):
    try:
        func(*args)
    except Exception as e:
        if error in str(e):
            return
        else:
            raise e


def test_functions():
    initialize()

    error_wrapper(save_atoms, error="did not match C++ signature")
    error_wrapper(load_contact_map, error="did not match C++ signature")
    error_wrapper(load_aligned_contact_map, error="did not match C++ signature")
