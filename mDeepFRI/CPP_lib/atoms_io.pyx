cimport cython

import numpy as np

cimport numpy as cnp
from libcpp.string cimport string

import ctypes

cnp.import_array()

cdef extern from "atoms_file_io.h":
    void SaveAtomsFile(float* positions,
                       size_t atom_count,
                       int* groups,
                       size_t chain_length,
                       string save_path)



def save_atoms_file(cnp.ndarray[cnp.float32_t, ndim=2] positions_array,
                    cnp.ndarray[cnp.int32_t, ndim=1] groups_array,
                    str save_path):


    atom_count = positions_array.shape[0]
    chain_length = groups_array.shape[0]

    __write_atoms_file(positions_array,
                       atom_count,
                       groups_array,
                       chain_length,
                       save_path.encode())




@cython.wraparound(False)
@cython.boundscheck(False)
cdef void __write_atoms_file(float[:, ::1] positions_array,
                             int atom_count,
                             int[:] groups_array,
                             int chain_length,
                             string save_path):

    cdef float* position_ptr = &positions_array[0, 0]
    cdef int* groups_ptr = &groups_array[0]

    SaveAtomsFile(position_ptr, atom_count,
                  groups_ptr, chain_length,
                  save_path)
