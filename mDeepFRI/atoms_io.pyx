cimport cython
from libc.stdlib cimport free, malloc
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np

cimport numpy as cnp

cnp.import_array()


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

# @cython.boundscheck(False)
# def create_contact_map(np.ndarray
#                        float angstrom_contact_threshold,
#                        str mode = "matrix"):

#     cdef float threshold_sq = angstrom_contact_threshold * angstrom_contact_threshold
#     cdef float dist

#     cdef int chain_length
#     cdef int[::1] groups
#     cdef float[::1] positions
#     positions, groups, chain_length = load_atoms_file(filepath.encode())

#     cdef int i, group_a, group_b, atom_a, atom_b
#     cdef bint connected
#     cdef cnp.ndarray[ndim=2, dtype=cnp.uint8_t]  contact_map \
#         = np.zeros((chain_length, chain_length), dtype=np.bool_)

#     cdef cnp.ndarray[ndim=2, dtype=cnp.uint32_t] non_empty

#     with nogil:
#         # fill diagonal with true
#         for group_a in range(chain_length):
#             contact_map[group_a, group_a] = True;

#             for group_b in range(group_a + 1, chain_length):
#                 connected = False
#                 for atom_a in range(groups[group_a], groups[group_a + 1]):
#                     for atom_b in range(groups[group_b], groups[group_b + 1]):
#                         dist = distance(positions, atom_a, atom_b)
#                         if dist <= threshold_sq:
#                             connected = True
#                             # fill the matrix using symmetry
#                             # reduces computational complexity
#                             contact_map[group_a, group_b] = True
#                             contact_map[group_b, group_a] = True
#                     if connected:
#                         break

#     # returns a binary contact map
#     if mode == "matrix":
#         return contact_map
#     # returns non-empty indices
#     if mode == "sparse":
#         return np.argwhere(contact_map == True).astype(np.uint32)


# def load_aligned_contact_map(str filepath,
#                              float angstrom_contact_threshold,
#                              str query_alignment,
#                              str target_alignment,
#                              int generated_contacts):

#     target_to_query_indexes = np.zeros(len(query_alignment), dtype=np.int32)
#     sparse_target_contacts = load_contact_map(filepath, 6, "sparse")
#     sparse_query_contacts = np.zeros((0, 2), dtype=np.uint32)

#     for i in range(len(query_alignment)):
#         if query_alignment[i] == "-":
#             target_to_query_indexes[target_index] = -1
#             target_index += 1
#         elif target_alignment[i] == "-":
#             for j in range(1, generated_contacts + 1):
#                 # add values to sparse_query_contacts
#                 first_coord = query_index - j
#                 second_coord = query_index + j
#                 if first_coord > 0:
#                     sparse_query_contacts = np.append(sparse_query_contacts,
#                                                         np.array([[first_coord, query_index]], dtype=np.int32), axis=0)

#                 sparse_query_contacts = np.append(sparse_query_contacts,
#                                                     np.array([[query_index + j, query_index]], dtype=np.int32), axis=0)
#             query_index += 1
#         else:
#             target_to_query_indexes[target_index] = query_index
#             query_index += 1
#             target_index += 1


#     for i in range(sparse_target_contacts.shape[0]):
#         contact_x = target_to_query_indexes[sparse_target_contacts[i, 0]]
#         if contact_x < 0:
#             continue
#         contact_y = target_to_query_indexes[sparse_target_contacts[i, 1]]
#         if contact_y < 0:
#             continue

#     sparse_query_contacts = np.vstack([sparse_query_contacts,
#                                       np.array([contact_x, contact_y])])

#     contact_map = np.zeros(query_index * query_index)
#     for i in range(query_index):
#         contact_map[i*query_index + i] = 1

#     for array in sparse_query_contacts:
#         if array[0] < 0:
#             continue
#         if array[0] >= query_index:
#             continue
#         contact_map[array[0]*query_index + array[1]] = 1
#         contact_map[array[1]*query_index + array[0]] = 1

#     contact_map = contact_map.reshape((query_index, query_index))

#     return contact_map

#     cdef cnp.ndarray[ndim=2, dtype=cnp.uint32_t] sparse_target_contacts = load_contact_map(filepath, angstrom_contact_threshold, "sparse")
#     # allocate a 2D array of unknown size
#     cdef vector[pair[int, int]] sparse_query_contacts = []
#     cdef int query_index = 0
#     cdef int target_index = 0
#     cdef int i, contact_x, contact_y

#     cdef int[::1] target_to_query_indexes = np.zeros(len(query_alignment), dtype=np.int32)
#     cdef cnp.uint32_t[::1] row

#     for i in range(len(query_alignment)):
#         if query_alignment[i] == "-":
#             target_to_query_indexes[target_index] = -1
#             target_index += 1
#         elif target_alignment[i] == "-":
#             for j in range(generated_contacts):
#                 sparse_query_contacts = np.append(sparse_query_contacts,
#                                                   np.array([[query_index - j, query_index]], dtype=np.uint32), axis=0)
#                 sparse_query_contacts = np.append(sparse_query_contacts,
#                                                   np.array([[query_index + j, query_index]], dtype=np.uint32), axis=0)
#             query_index += 1
#         else:
#             target_to_query_indexes[target_index] = query_index
#             query_index += 1
#             target_index += 1

#     cdef cnp.ndarray[ndim=2, dtype=cnp.uint8_t]  contact_map \
#         = np.zeros((query_index, query_index), dtype=np.bool_)

#     sparse_target_contacts = sparse_target_contacts.copy(order="C")
#     for row in sparse_target_contacts:
#         contact_x = target_to_query_indexes[row[0]]
#         contact_y = target_to_query_indexes[row[1]]

#         if contact_x != -1 or contact_y != -1:
#             continue

#         sparse_query_contacts = np.append(sparse_query_contacts,
#                                           np.array([[contact_x, contact_y]], dtype=np.uint32), axis=0)


#     for i in range(query_index):
#         contact_map[i, i] = True

#     for row in sparse_query_contacts:
#         if row[0] < 0 or row[0] >= query_index:
#             continue
#         contact_map[row[0], row[1]] = True
#         contact_map[row[1], row[0]] = True

#     return sparse_target_contacts
