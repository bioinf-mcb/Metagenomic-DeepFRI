# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

cimport cython

import numpy as np

cimport numpy as np
from libc.math cimport sqrt
from libc.stdlib cimport free, malloc
from libc.string cimport strlen


@cython.boundscheck(False)
@cython.wraparound(False)
cdef float alignment_sequences_identity(char* query, char* target) nogil:
    cdef int matches = 0
    cdef int i
    cdef float seq_identity
    cdef int query_len = strlen(query)

    # count matches excluding gaps
    for i in range(query_len):
        if query[i] == target[i] and query[i] != b'-' and target[i] != b'-':
            matches += 1

    # calculate identity
    seq_identity = matches / query_len
    return seq_identity


def alignment_identity(str query, str target):
    """
    Calculate sequence identity between two proteins in an alignment.
    Gaps are excluded. The length of sequence is irrelevant, because
    lengths should be identical.

    Args:
        query (str): Query sequence.
        target (str): Target sequence.

    Returns:
        float: Sequence identity between two sequences in an alignment.
    """
    return alignment_sequences_identity(query.encode(), target.encode())


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
