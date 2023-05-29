cimport cython
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
