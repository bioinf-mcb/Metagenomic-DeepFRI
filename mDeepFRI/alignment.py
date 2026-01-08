"""
Sequence alignment module using PyOpal for protein alignments.

This module provides utilities for performing pairwise protein sequence alignments
using PyOpal, a Python wrapper for the SIMD-accelerated Opal alignment library.
It handles alignment of query sequences to database hits from MMseqs2 searches
and provides alignment statistics including identity, coverage, and gaps.

The module optimizes alignment by batching sequences and using parallel processing,
making it suitable for aligning large numbers of protein sequences.

Classes:
    AlignmentResult: Container for pairwise alignment results and statistics.

Functions:
    insert_gaps: Insert gaps into sequences based on alignment string.
    _pyopal_align: Internal function for PyOpal alignment.
    align_batch: Align batch of sequences in parallel.
    align_mmseqs_results: Align all MMseqs2 hits to their targets.
"""

import warnings
from functools import partial
from multiprocessing import Pool
from typing import Optional, Tuple

import numpy as np
import pyopal
from scoring_matrices import ScoringMatrix

from mDeepFRI.mmseqs import MMseqsResult
from mDeepFRI.utils import (load_fasta_as_dict, retrieve_fasta_entries_as_dict,
                            stdout_warn)

warnings.showwarning = stdout_warn


def insert_gaps(sequence: str, reference: str,
                alignment_string: str) -> Tuple[str, str]:
    """
    Inserts gaps into query and target sequences.

    Args:
        sequence (str): Query sequence.
        reference (str): Target sequence.
        alignment_string (str): Alignment string.

    Returns:
        gapped_sequence (str): Query sequence with gaps.
        gapped_target (str): Target sequence with gaps.
    """

    sequence_list: list[str] = list(sequence)
    reference_list: list[str] = list(reference)
    alignment_list: list[str] = list(alignment_string)

    for i, a in enumerate(alignment_list):
        if a == "I":
            sequence_list.insert(i, "-")
        elif a == "D":
            reference_list.insert(i, "-")
    return "".join(sequence_list), "".join(reference_list)


class AlignmentResult:
    """
    Container for pairwise protein alignment results and statistics.

    This class stores the results of a protein sequence alignment, including
    the aligned sequences with gaps, alignment statistics (identity, coverage),
    and optional structural information (coordinates, contact maps).

    Attributes:
        query_name (str): Identifier of the query sequence.
        query_sequence (str): Ungapped query protein sequence.
        target_name (str): Identifier of the target (reference) sequence.
        target_sequence (str): Ungapped target protein sequence.
        alignment (str): Alignment string in CIGAR-like format.
            'M' = match/mismatch, 'I' = insertion, 'D' = deletion.
        query_identity (float): Sequence identity as fraction (0.0-1.0) of
            matching residues to alignment length.
        query_coverage (float): Fraction (0.0-1.0) of query sequence covered
            by the alignment.
        target_coverage (float): Fraction (0.0-1.0) of target sequence covered
            by the alignment.
        db_name (str): Name of the database from which target was retrieved.
        gapped_sequence (str): Query sequence with gaps ('-') inserted for alignment.
        gapped_target (str): Target sequence with gaps ('-') inserted for alignment.
        target_coords (np.ndarray, optional): C-alpha atom coordinates from target structure.
        cmap (np.ndarray, optional): Contact map of target structure.
        aligned_cmap (np.ndarray, optional): Contact map aligned to query sequence.

    Example:
        >>> result = AlignmentResult(
        ...     query_name="protein1",
        ...     query_sequence="MSKGEELFT",
        ...     target_name="1GFL_A",
        ...     target_sequence="MSKGEELFTGV",
        ...     alignment="MMMMMMMMMM",
        ...     query_identity=0.90,
        ...     query_coverage=0.82
        ... )
        >>> print(result.gapped_sequence)
        'MSKGEELFT'
    """
    def __init__(self,
                 query_name: str = "",
                 query_sequence: str = "",
                 target_name: str = "",
                 target_sequence: str = "",
                 alignment: str = "",
                 query_identity: Optional[float] = None,
                 query_coverage: Optional[float] = None,
                 target_coverage: Optional[float] = None,
                 db_name: Optional[str] = None,
                 coords: Optional[np.ndarray] = None):

        self.query_name = query_name
        self.query_sequence = query_sequence
        self.target_name = target_name
        self.target_sequence = target_sequence
        self.alignment = alignment
        self.query_identity = query_identity
        self.query_coverage = query_coverage
        self.target_coverage = target_coverage
        self.insert_gaps()
        self.db_name = db_name
        self.coords = coords
        self.target_coords = None
        self.cmap = None
        self.aligned_cmap = None

    def __str__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"query_identity={self.query_identity}, query_coverage={self.query_coverage})"

    def __repr__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"query_identity={self.query_identity}, query_coverage={self.query_coverage})"

    def insert_gaps(self):
        """
        Inserts gaps into query and target sequences.

        Returns:
            AlignmentResult: The object with gapped sequences.
        """

        self.gapped_sequence, self.gapped_target = insert_gaps(
            self.query_sequence, self.target_sequence, self.alignment)


def best_hit_database(query,
                      target_sequences,
                      gap_open: int = 10,
                      gap_extend: int = 1,
                      scoring_matrix: str = "VTML80"):
    """
    Find the best hit in the database and return index.

    Args:
        query (str): The query sequence.
        target_sequences (dict): The target sequences.

    Returns:
        int: The index of the best hit.
        str: The best hit sequence.
    """

    custom_scoring = ScoringMatrix.from_name(scoring_matrix)
    custom_alphabet = custom_scoring.alphabet
    aligner = pyopal.Aligner(scoring_matrix=custom_scoring,
                             gap_open=gap_open,
                             gap_extend=gap_extend)
    target_database = pyopal.Database(target_sequences.values(),
                                      alphabet=custom_alphabet)

    # Retrieve the best hit
    results = aligner.align(query,
                            target_database,
                            mode="score",
                            overflow="buckets",
                            algorithm="nw")
    best_hit = max(results, key=lambda x: x.score)
    best_index = best_hit.target_index
    best_idx = list(target_sequences.keys())[best_index]
    best_seq = target_sequences[best_idx]

    return best_idx, best_seq


def align_pairwise(query,
                   target,
                   gap_open: int = 10,
                   gap_extend: int = 1,
                   scoring_matrix: str = "VTML80"):
    """
    Aligns the query against the target and returns the alignment.

    Args:
        query (str): The query sequence.
        target (str): The target sequence.
        aligner (pyopal.Aligner): The aligner object.

    Returns:
        str: The alignment of the query against the target.
        float: The identity of the alignment.

    """

    custom_scoring = ScoringMatrix.from_name(scoring_matrix)
    custom_alphabet = custom_scoring.alphabet
    aligner = pyopal.Aligner(scoring_matrix=custom_scoring,
                             gap_open=gap_open,
                             gap_extend=gap_extend)
    database = pyopal.Database([target], alphabet=custom_alphabet)
    # Align the sequences
    alignment = aligner.align(query, database, algorithm="nw", mode="full")
    alignment_string = alignment[0].alignment
    identity = alignment[0].identity()
    query_coverage = alignment[0].coverage(reference="query")
    target_coverage = alignment[0].coverage(reference="target")

    return alignment_string, identity, query_coverage, target_coverage


def pairwise_against_database(query_id,
                              query_sequence,
                              target_sequences,
                              gap_open: int = 10,
                              gap_extend: int = 1,
                              scoring_matrix: str = "VTML80"):
    """
    Finds the best alignment of the query against the target.
    """

    best_idx, best_target = best_hit_database(query_sequence, target_sequences,
                                              gap_open, gap_extend,
                                              scoring_matrix)

    # align the query against the best hit
    alignment, identity, query_coverage, target_coverage = align_pairwise(
        query_sequence, best_target, gap_open, gap_extend, scoring_matrix)
    # create an alignment object
    alignment_result = AlignmentResult(query_id,
                                       query_sequence,
                                       best_idx,
                                       best_target,
                                       alignment,
                                       identity,
                                       query_coverage=query_coverage,
                                       target_coverage=target_coverage)
    return alignment_result


def create_partial_database(top_matches, target_sequences):
    """
    Create a partial database from the target sequences.

    Args:

    """

    partial_db = {k: target_sequences[k] for k in top_matches}

    return partial_db


def align_mmseqs_results(best_matches_filepath: str,
                         sequence_db: str,
                         alignment_gap_open: int = 10,
                         alignment_gap_extend: int = 1,
                         threads: int = 1,
                         scoring_matrix: str = "VTML80"):
    """
    Aligns MMseqs2 search results sequence-wise.

    This function takes the best matches from MMseqs2, retrieves the corresponding
    target sequences from the database, and performs pairwise alignment using PyOpal.
    It returns a list of AlignmentResult objects containing alignment statistics
    and gapped sequences.

    Args:
        best_matches_filepath (str): Path to MMseqs2 best matches TSV file.
        sequence_db (str): Path to FASTA database of target sequences.
        alignment_gap_open (int, optional): Gap open penalty. Defaults to 10.
        alignment_gap_extend (int, optional): Gap extension penalty. Defaults to 1.
        threads (int, optional): Number of threads for parallel processing.
            Defaults to 1.
        scoring_matrix (str, optional): Scoring matrix for alignment.
            Defaults to "VTML80".

    Returns:
        List[AlignmentResult]: List of alignment results for each query.
    """

    best_matches = MMseqsResult.from_best_matches(best_matches_filepath)
    # check if there are any best matches
    if best_matches.size == 0:
        return []

    query_dict = load_fasta_as_dict(best_matches.query_fasta)

    # if UniProt header is used, extract second part as sequence ID
    for qid in list(query_dict.keys()):
        if "|" in qid:
            new_qid = qid.split("|")[1]
            query_dict[new_qid] = query_dict.pop(qid)

    unique_queries = {
        query: best_matches.get_query_targets(query)
        for query in best_matches.get_queries()
    }

    target_ids = best_matches.get_targets()
    target_seqs = retrieve_fasta_entries_as_dict(sequence_db, target_ids)

    # create partial databases
    partial_databases = {
        k: create_partial_database(v, target_seqs)
        for k, v in unique_queries.items()
    }

    align_partial = partial(pairwise_against_database,
                            gap_open=alignment_gap_open,
                            gap_extend=alignment_gap_extend,
                            scoring_matrix=scoring_matrix)

    query_ids = list(unique_queries.keys())
    print(query_ids)
    print(query_dict.keys())
    query_sequences = [query_dict[qid] for qid in query_ids]
    target_sequences = [partial_databases[qid] for qid in query_ids]

    # align the query against the best hit
    with Pool(threads) as pool:
        alignments = pool.starmap(
            align_partial, [(query_id, query_sequence, target_sequence)
                            for query_id, query_sequence, target_sequence in
                            zip(query_ids, query_sequences, target_sequences)])

    return alignments
