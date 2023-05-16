cimport cython
from libc.stdlib cimport free, malloc
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp

cnp.import_array()

## TODO: after connecting FoldComp remove bit files load from the function scope

cdef extern from "atoms_file_io.h" nogil:
    void SaveAtomsFile(float* positions,
                       size_t atom_count,
                       int* groups,
                       size_t chain_length,
                       string save_path)

cdef class Coordinate:
    cdef readonly int i, j

    def __init__(self, int i, int j):
        self.i = i
        self.j = j


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


cdef inline load_atoms_file(string filepath):
    cdef int chain_length = np.fromfile(filepath, count=1, dtype=np.int32)

    # offset 4 bytes because we read chain_length already
    cdef int[::1] groups = np.fromfile(filepath, count=chain_length, dtype=np.int32,
                                       offset=4, sep='')

    cdef int positions_offset = 4 * (1 + chain_length)

    cdef float[::1] positions = np.fromfile(filepath, count=groups[-1] * 3,
                                            dtype=np.float32, offset=positions_offset, sep='')

    chain_length = chain_length - 1

    return positions, groups, chain_length


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline float distance(float[::1] array, int i, int j) nogil:
    """
    Calculates the square distance between two atoms.
    Computationally efficient replacement for Euclidian
    distance, as sqrt() is computationally heavy.
    """

    cdef float x = array[i * 3] - array[j * 3]
    cdef float y = array[i * 3 + 1] - array[j * 3 + 1]
    cdef float z = array[i * 3 + 2] - array[j * 3 + 2]
    cdef float distance = x * x + y * y + z * z

    return distance

## TODO: modularize the function
@cython.boundscheck(False)
def load_contact_map(str filepath,
                     float angstrom_contact_threshold,
                     str mode = "matrix"):

    cdef float threshold_sq = angstrom_contact_threshold * angstrom_contact_threshold
    cdef float dist

    cdef int chain_length
    cdef int[::1] groups
    cdef float[::1] positions
    positions, groups, chain_length = load_atoms_file(filepath.encode())

    cdef int i, group_a, group_b, atom_a, atom_b
    cdef bint connected
    cdef cnp.ndarray[ndim=2, dtype=cnp.uint8_t]  contact_map \
        = np.zeros((chain_length, chain_length), dtype=np.bool_)

    cdef cnp.ndarray[ndim=2, dtype=cnp.uint32_t] non_empty

    with nogil:
        # fill diagonal with true
        for group_a in range(chain_length):
            contact_map[group_a, group_a] = True;

            for group_b in range(group_a + 1, chain_length):
                connected = False
                for atom_a in range(groups[group_a], groups[group_a + 1]):
                    for atom_b in range(groups[group_b], groups[group_b + 1]):
                        dist = distance(positions, atom_a, atom_b)
                        if dist <= threshold_sq:
                            connected = True
                            # fill the matrix using symmetry
                            # reduces computational complexity
                            contact_map[group_a, group_b] = True
                            contact_map[group_b, group_a] = True
                    if connected:
                        break

    # returns a binary contact map
    if mode == "matrix":
        return contact_map
    # returns non-empty indices
    if mode == "sparse":
        return np.argwhere(contact_map == True).astype(np.uint32)


def load_aligned_contact_map(str filepath,
                             float angstrom_contact_threshold,
                             str query_alignment,
                             str target_alignment,
                             int generated_contacts):

    cdef cnp.ndarray[ndim=2, dtype=cnp.uint32_t] sparse_target_contacts = load_contact_map(filepath, angstrom_contact_threshold, "sparse")
    # allocate a 2D array of unknown size
    cdef cnp.ndarray[ndim=2, dtype=cnp.uint32_t, mode="c"] sparse_query_contacts = np.empty((0, 2), dtype=np.uint32)
    cdef int query_index = 0
    cdef int target_index = 0
    cdef int i, contact_x, contact_y

    cdef int[::1] target_to_query_indexes = np.zeros(len(query_alignment), dtype=np.int32)
    cdef cnp.uint32_t[::1] row

    for i in range(len(query_alignment)):
        if query_alignment[i] == "-":
            target_to_query_indexes[target_index] = -1
            target_index += 1
        elif target_alignment[i] == "-":
            for j in range(generated_contacts):
                sparse_query_contacts = np.append(sparse_query_contacts,
                                                  np.array([[query_index - j, query_index]], dtype=np.uint32), axis=0)
                sparse_query_contacts = np.append(sparse_query_contacts,
                                                  np.array([[query_index + j, query_index]], dtype=np.uint32), axis=0)
            query_index += 1
        else:
            target_to_query_indexes[target_index] = query_index
            query_index += 1
            target_index += 1

    cdef cnp.ndarray[ndim=2, dtype=cnp.uint8_t]  contact_map \
        = np.zeros((query_index, query_index), dtype=np.bool_)

    sparse_target_contacts = sparse_target_contacts.copy(order="C")
    for row in sparse_target_contacts:
        contact_x = target_to_query_indexes[row[0]]
        contact_y = target_to_query_indexes[row[1]]

        if contact_x != -1 or contact_y != -1:
            continue

        sparse_query_contacts = np.append(sparse_query_contacts,
                                          np.array([[contact_x, contact_y]], dtype=np.uint32), axis=0)


    for i in range(query_index):
        contact_map[i, i] = True

    for row in sparse_query_contacts:
        if row[0] < 0 or row[0] >= query_index:
            continue
        contact_map[row[0], row[1]] = True
        contact_map[row[1], row[0]] = True

    return sparse_target_contacts
