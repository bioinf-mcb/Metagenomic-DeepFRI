

import numpy as np

cimport cython
cimport numpy as cnp
from libcpp.vector cimport vector

# Define types for clarity and speed
ctypedef cnp.int32_t DTYPE_t
from cython.parallel import prange


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef pairwise_sqeuclidean(float[:, ::1] X):

    cdef int n = X.shape[0]
    cdef int m = X.shape[1]
    cdef int i, j, k
    cdef float d, diff

    cdef float[:, ::1] D = np.zeros((n, n), dtype=np.float32)

    with nogil:
        for i in prange(n, schedule='static', chunksize=1):
            for j in range(i + 1, n):
                d = 0.0
                for k in range(m):
                    diff = X[i, k] - X[j, k]
                    d = d + (diff * diff)

                D[i, j] = d
                D[j, i] = d

    return np.asarray(D)


cpdef cnp.ndarray[DTYPE_t, ndim=2] align_contact_map(str query_alignment,
                        str target_alignment,
                        cnp.ndarray[DTYPE_t, ndim=2] sparse_target_contact_map,
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

    # 1. Convert Python strings to C-strings (char*) for O(1) access
    cdef bytes query_bytes = query_alignment.encode('ascii')
    cdef bytes target_bytes = target_alignment.encode('ascii')
    cdef char* q_ptr = query_bytes
    cdef char* t_ptr = target_bytes
    cdef int align_len = len(query_bytes)

    # 2. Use C++ Vectors for dynamic memory management
    # Replaces manual malloc/free and avoids pre-allocation guessing.
    # target_to_query_map stores the query index for every target index.
    cdef vector[int] target_to_query_map
    target_to_query_map.reserve(align_len)

    # flattened_contacts stores pairs (i, j) flatly: [i1, j1, i2, j2...]
    cdef vector[int] sparse_query_contacts
    sparse_query_contacts.reserve(align_len * 2)

    cdef int query_idx = 0
    cdef int target_idx = 0
    cdef int i, j, k, row

    # 3. First Pass: Compute Mapping and Generated Contacts
    for i in range(align_len):
        if q_ptr[i] == 45: # 45 is ASCII for '-'
            # Gap in Query: Target has residue, Query does not.
            # Map current target_idx to -1 (no equivalent in query).
            target_to_query_map.push_back(-1)
            target_idx += 1
        else:
            # Query has residue
            if t_ptr[i] == 45: # Gap in Target
                # Generated contacts logic for insertions in Query
                for j in range(1, generated_contacts + 1):
                    # Add (q + j, q)
                    sparse_query_contacts.push_back(query_idx + j)
                    sparse_query_contacts.push_back(query_idx)
                    # Add (q - j, q)
                    sparse_query_contacts.push_back(query_idx - j)
                    sparse_query_contacts.push_back(query_idx)

                query_idx += 1
            else:
                # Match or Mismatch: Both have residues
                target_to_query_map.push_back(query_idx)
                query_idx += 1
                target_idx += 1

    # 4. Second Pass: Translate Target Contacts
    # Access numpy array via memoryview for raw C speed
    cdef int[:, :] target_contacts_view = sparse_target_contact_map
    cdef int num_target_contacts = sparse_target_contact_map.shape[0]
    cdef int t_res_i, t_res_j
    cdef int q_res_i, q_res_j

    for row in range(num_target_contacts):
        t_res_i = target_contacts_view[row, 0]
        t_res_j = target_contacts_view[row, 1]

        # Bounds check for safety (optional, but good practice with raw vectors)
        if t_res_i < target_to_query_map.size() and t_res_j < target_to_query_map.size():
            q_res_i = target_to_query_map[t_res_i]
            q_res_j = target_to_query_map[t_res_j]

            # If both residues map to valid query positions
            if q_res_i != -1 and q_res_j != -1:
                sparse_query_contacts.push_back(q_res_i)
                sparse_query_contacts.push_back(q_res_j)

    # 5. Build Output Numpy Array
    # We allocate this only once we know the final query_idx size
    cdef cnp.ndarray[DTYPE_t, ndim=2] output_map = np.zeros((query_idx, query_idx), dtype=np.int32)

    # Fill diagonal
    for i in range(query_idx):
        output_map[i, i] = 1

    # Fill contacts from the sparse vector
    cdef int p1, p2
    cdef size_t total_sparse_entries = sparse_query_contacts.size()

    for i in range(0, total_sparse_entries, 2):
        p1 = sparse_query_contacts[i]
        p2 = sparse_query_contacts[i+1]

        # Boundary checks are crucial for generated contacts (which can go out of bounds)
        if 0 <= p1 < query_idx and 0 <= p2 < query_idx:
            output_map[p1, p2] = 1
            output_map[p2, p1] = 1

    return output_map
