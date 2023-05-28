from libc.string cimport strlen


cdef alignment_sequences_identity(char* cigar):
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
    cdef int aln_len = strlen(cigar)

    # count matches excluding gaps
    for i in range(aln_len):
        if cigar[i] == b"M":
            matches += 1

    # calculate identity
    seq_identity = matches / aln_len
    return seq_identity
