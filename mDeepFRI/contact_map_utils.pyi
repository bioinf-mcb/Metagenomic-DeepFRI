"""Type stubs for contact_map_utils module."""

import numpy as np
from numpy.typing import NDArray

def pairwise_sqeuclidean(X: NDArray[np.float32]) -> NDArray[np.float32]:
    """
    Compute pairwise squared Euclidean distances between rows of X.

    Args:
        X: Input array of shape (n, m) where n is the number of points
            and m is the dimensionality.

    Returns:
        Symmetric distance matrix of shape (n, n) containing squared
        Euclidean distances.

    Note:
        This is a highly optimized Cython implementation that uses
        parallel processing with OpenMP for better performance.
    """
    ...

def align_contact_map(
    query_alignment: str,
    target_alignment: str,
    sparse_target_contact_map: NDArray[np.int32],
    generated_contacts: int = 2
) -> NDArray[np.int32]:
    """
    Align a contact map based on sequence alignments.

    Aligns a contact map from a target structure to match a query sequence
    based on their alignment. Handles gaps by generating synthetic contacts
    for inserted regions.

    Args:
        query_alignment: Aligned query sequence (with gaps as '-').
        target_alignment: Aligned target sequence (with gaps as '-').
        sparse_target_contact_map: Sparse contact map of target as array
            of shape (N, 2) where each row [i, j] indicates residues i and j
            are in contact.
        generated_contacts: Number of contacts to generate for gapped regions
            in the query. Defaults to 2.

    Returns:
        Aligned contact map as binary matrix of shape (L, L) where L is the
        query sequence length (excluding gaps). Values are 1 for contacts
        and 0 for non-contacts.

    Algorithm:
        1. Initialize mapping from target to query residue indices
        2. Iterate through alignment positions:
           - Gap in query: skip target residue
           - Gap in target: generate synthetic contacts for query insertion
           - Match/mismatch: record the mapping
        3. Translate target contacts to query coordinates using mapping
        4. Filter out contacts that don't map to valid query positions
        5. Build final contact map matrix with diagonal set to 1
        6. Populate with translated and generated contacts symmetrically

    Note:
        - Input alignments must be the same length
        - Contact indices in sparse_target_contact_map are 0-based
        - Generated contacts help maintain structural continuity in insertions
    """
    ...
