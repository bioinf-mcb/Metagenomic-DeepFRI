cimport cython

import numpy as np

cimport numpy as np
from libc.stdlib cimport free, malloc
from libc.string cimport strlen


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef pairwise_sqeuclidean(float[:, ::1] X):

    cdef int n = X.shape[0]
    cdef int m = X.shape[1]
    cdef int i, j, k
    cdef float d, diff
    cdef float[:, ::1] D = np.zeros((n, n), dtype=np.float32)

    with nogil:
        for i in range(n):
            for j in range(i + 1, n):
                d = 0
                for k in range(m):
                    diff = X[i, k] - X[j, k]
                    d += diff * diff
                D[i, j] = d
                D[j, i] = d

    return np.asarray(D)


cpdef align_contact_map(str query_alignment,
                        str target_alignment,
                        np.ndarray[np.int32_t, ndim=2] sparse_target_contact_map,
                        int generated_contacts=2):
    """
    Aligns a contact map based on the alignments of query and target sequences.

    Args:
        query_alignment: The alignment of the query sequence.
        target_alignment: The alignment of the target sequence.
        sparse_target_contact_map: The sparse contact map of the target
                                   sequence represented as a list of tuples (i, j)
                                   indicating contacts between residues iand j.
        generated_contacts: The number of generated contacts to add for gapped
                            regions in the query alignment. Defaults to 2.

    Returns:
        The aligned contact map as a numpy array.

    Algorithm:
    1. Initialize an empty list `sparse_query_contact_map` to store the contacts in the aligned contact map.
    2. Initialize variables `target_index` and `query_index` to track the indices of residues in the target
    and query proteins, respectively.
    3. Initialize an empty dictionary `target_to_query_indices` to map target residues to query residues
    using shift resulting from the alignments.
    4. Iterate over each position in the query alignment:
        - If the query residue is '-', increment the `target_index` and do not add a contact
        to the aligned contact map.
        - If the query residue is not '-', check the target residue:
            - If the target residue is '-', add contacts for the generated region in the query alignment:
                - For each generated contact, add the contact (query_index + j, query_index)
                and (query_index - j, query_index) to the `sparse_query_contact_map` according to generated_contacts.
                - Increment the `query_index`.
            - If the target residue is not '-', map the target residue to the query residue by adding
            an entry in the `target_to_query_indices` dictionary.
                - Increment both the `query_index` and `target_index`.
    5. Translate the target residue indices to query residue indices
    in the `sparse_target_contact_map` by using the `target_to_query_indices` dictionary.
    6. Filter out the contacts that are not present in the query alignment by removing contacts
    with '-1' indices from the `sparse_target_contact_map`.
    7. Add the filtered contacts from the filtered `sparse_target_contact_map` to the `sparse_query_contact_map`.
    8. Build the output contact map with dimensions (query_index, query_index) initialized as all zeros.
    Query index is the number of residues in the query sequence.
    9. Set the diagonal elements of the output contact map to 1.
    10. For each contact (i, j) in the `sparse_query_contact_map`:
        - If i is less than 0 or greater than or equal to `query_index`, skip the contact.
        - Otherwise, set the corresponding elements in the output contact map to 1 symmetrically.
    11. Return the aligned contact map as a numpy array.
    """

    cdef int target_index = 0
    cdef int query_index = 0
    cdef int i, j, k, n
    cdef int *target_to_query_indices = <int *>malloc(len(target_alignment) * sizeof(int))
    cdef int *sparse_query_contact_map = <int *>malloc(len(query_alignment) * len(target_alignment) * 2 * sizeof(int))
    cdef int sparse_map_size = 0
    cdef int output_contact_map_size = query_index * query_index

    # Map target residues to query residues based on the alignments
    for i in range(len(query_alignment)):
        if query_alignment[i] == "-":
            target_to_query_indices[target_index] = -1
            target_index += 1
        else:
            if target_alignment[i] == "-":
                for j in range(1, generated_contacts + 1):
                    sparse_query_contact_map[sparse_map_size] = query_index + j
                    sparse_query_contact_map[sparse_map_size + 1] = query_index
                    sparse_map_size += 2
                    sparse_query_contact_map[sparse_map_size] = query_index - j
                    sparse_query_contact_map[sparse_map_size + 1] = query_index
                    sparse_map_size += 2
                query_index += 1
            else:
                target_to_query_indices[target_index] = query_index
                query_index += 1
                target_index += 1

    # Translate the target residues index to query residues index
    for i in range(sparse_target_contact_map.shape[0]):
        if (target_to_query_indices[sparse_target_contact_map[i, 0]] != -1 and
            target_to_query_indices[sparse_target_contact_map[i, 1]] != -1):
            sparse_query_contact_map[sparse_map_size] = target_to_query_indices[sparse_target_contact_map[i, 0]]
            sparse_query_contact_map[sparse_map_size + 1] = target_to_query_indices[sparse_target_contact_map[i, 1]]
            sparse_map_size += 2

    # Build the output contact map
    cdef np.ndarray[np.int32_t, ndim=2] output_contact_map = np.zeros((query_index, query_index), dtype=np.int32)

    # Fill the diagonal
    for i in range(query_index):
        output_contact_map[i, i] = 1

    # Fill the contacts from the sparse query contact map
    for i in range(0, sparse_map_size, 2):
        j = sparse_query_contact_map[i]
        k = sparse_query_contact_map[i + 1]
        if j < query_index and k < query_index:
            output_contact_map[j, k] = 1
            output_contact_map[k, j] = 1

    free(target_to_query_indices)
    free(sparse_query_contact_map)

    return output_contact_map
