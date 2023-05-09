cimport cython
from libc.stdlib cimport free, malloc
from libcpp.string cimport string

import numpy as np

cimport numpy as cnp

cnp.import_array()

cdef extern from "atoms_file_io.h" nogil:
    void SaveAtomsFile(float* positions,
                       size_t atom_count,
                       int* groups,
                       size_t chain_length,
                       string save_path)



@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def save_atoms_file(list positions_list,
                    list groups_list,
                    str save_path):

    # every atom has 3 coordinates (x, y & z)
    cdef:
        int atom_count = int(len(positions_list) / 3)
        int chain_length = int(len(groups_list))
        float * positions_array = <float *>malloc(atom_count * 3 * sizeof(float))
        int * groups_array = <int *>malloc(chain_length * sizeof(int))
        string save_path_bytes = save_path.encode('utf-8')
        int i

    # copy the lists to arrays
    for i in range(atom_count):
        positions_array[i * 3] = positions_list[i * 3]
        positions_array[i * 3 + 1] = positions_list[i * 3 + 1]
        positions_array[i * 3 + 2] = positions_list[i * 3 + 2]

    for i in range(chain_length):
        groups_array[i] = groups_list[i]

    # save the file
    with nogil:
        SaveAtomsFile(positions_array, atom_count,
                      groups_array, chain_length,
                      save_path_bytes)


def load_atom_file(filepath):
    cdef int chain_length = np.fromfile(filepath, count=1, dtype=np.int32)
    cdef cnp.ndarray[cnp.int32_t, ndim=1] groups = np.fromfile(filepath, count=chain_length, dtype=np.int32,
                                                               offset=4, sep='')
    cdef int positions_offset = 4 * (1 + chain_length)

    cdef cnp.ndarray[cnp.float32_t, ndim=1] positions = np.fromfile(filepath, count=groups[-1] * 3,
                                                                    dtype=np.float32, offset=positions_offset, sep='')

    return positions, groups, chain_length
