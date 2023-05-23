from libc.string cimport strlen


cdef alignment_sequences_identity(char* query, char* target) -> float:
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
    cdef int matches = 0
    cdef int i
    cdef float seq_identity
    cdef int query_len = strlen(query)

    # count matches excluding gaps
    for i in range(query_len):
        if query[i] == target[i] and query[i] != b'-' and target[i] != b'-':
            matches += 1

    # calculate identity
    seq_identity = matches / len(query)
    return seq_identity
