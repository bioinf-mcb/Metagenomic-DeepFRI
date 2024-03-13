import logging
from collections import defaultdict
from functools import partial
from multiprocessing.pool import Pool

import numpy as np
import pyopal

from mDeepFRI.bio_utils import (AlignmentResult, load_fasta_as_dict,
                                retrieve_fasta_entries_as_dict)
from mDeepFRI.mmseqs import (filter_mmseqs_results, run_mmseqs_search,
                             validate_mmseqs_database)

logger = logging.getLogger(__name__)

##TODO: do not collect large database
##TODO: align 1 query against top-k databases


def best_hit_database(query,
                      target_sequences,
                      gap_open: int = 10,
                      gap_extend: int = 1):
    """
    Find the best hit in the database and return index.

    Args:
        query (str): The query sequence.
        target_sequences (dict): The target sequences.

    Returns:
        int: The index of the best hit in the database.
    """

    aligner = pyopal.Aligner(gap_open=gap_open, gap_extend=gap_extend)
    target_database = pyopal.Database(target_sequences.values())

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


def align_pairwise(query, target, gap_open: int = 10, gap_extend: int = 1):
    """
    Aligns the query against the target and returns the alignment.

    Args:
        query (str): The query sequence.
        target (str): The target sequence.
        aligner (pyopal.Aligner): The aligner object.

    Returns:
        pyopal.Alignment: The alignment of the query against the target.
    """

    aligner = pyopal.Aligner(gap_open=gap_open, gap_extend=gap_extend)
    database = pyopal.Database([target])
    # Align the sequences
    alignment = aligner.align(query, database, algorithm="nw", mode="full")

    return alignment[0].alignment


def pairwise_against_database(query_id,
                              query_sequence,
                              target_sequences,
                              gap_open: int = 10,
                              gap_extend: int = 1):
    """
    Finds the best alignment of the query against the target.
    """

    # find best hits in the database
    best_idx, best_target = best_hit_database(query_sequence, target_sequences,
                                              gap_open, gap_extend)
    # align the query against the best hit
    alignment = align_pairwise(query_sequence, best_target, gap_open,
                               gap_extend)
    # create an alignment object
    alignment_result = AlignmentResult(query_id, query_sequence, best_idx,
                                       best_target, alignment)

    return alignment_result


def create_partial_database(top_matches, target_sequences):
    """
    Create a partial database from the target sequences.

    Args:
        query_id (str): The
    """

    partial_db = {k: target_sequences[k] for k in top_matches}

    return partial_db


def run_alignment(query_file: str,
                  mmseqs_db: str,
                  sequence_db: str,
                  output_path: str,
                  mmseqs_min_bitscore: int = None,
                  mmseqs_max_eval: float = 10e-5,
                  mmseqs_min_identity: float = 0.3,
                  top_k: int = 5,
                  alignment_gap_open: int = 10,
                  alignment_gap_extend: int = 1,
                  threads: int = 1):

    align_partial = partial(pairwise_against_database,
                            gap_open=alignment_gap_open,
                            gap_extend=alignment_gap_extend)
    mmseqs_valid = validate_mmseqs_database(mmseqs_db)
    if not mmseqs_valid:
        raise ValueError(f"MMSeqs2 database not found in {mmseqs_db}.")

    mmseqs_results = run_mmseqs_search(query_file, mmseqs_db, output_path,
                                       threads)
    filtered_mmseqs_results = filter_mmseqs_results(mmseqs_results,
                                                    mmseqs_min_bitscore,
                                                    mmseqs_max_eval,
                                                    mmseqs_min_identity,
                                                    k_best_hits=top_k,
                                                    threads=threads)

    # if MMSeqs2 alignment is empty
    if filtered_mmseqs_results is None:
        alignments = []
    else:
        query_seqs = load_fasta_as_dict(query_file)
        unique_queries = defaultdict(list)
        for result in filtered_mmseqs_results:
            unique_queries[result[0]].append(result[1])

        target_ids = np.unique(filtered_mmseqs_results["target"]).tolist()
        target_seqs = retrieve_fasta_entries_as_dict(sequence_db, target_ids)
        # create partial databases
        partial_databases = {
            k: create_partial_database(v, target_seqs)
            for k, v in unique_queries.items()
        }

        # align the query against the best hit
        with Pool(threads) as pool:
            alignments = pool.starmap(
                align_partial,
                [(query_id, query_seqs[query_id], partial_databases[query_id])
                 for query_id in unique_queries])

    return alignments
