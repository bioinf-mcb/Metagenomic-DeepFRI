import json
from functools import partial
import pandas as pd
import pathlib
import parasail
from multiprocessing.pool import ThreadPool

import logging

from meta_deepFRI.config.names import ALIGNMENTS
from meta_deepFRI.utils.fasta_file_io import SeqFileLoader
from meta_deepFRI.utils.bio_utils import substitution_matrices


def alignment_sequences_identity(query: str, target: str) -> float:
    """
    Calculate sequence identity between two proteins in an alignment.
    Gaps are excluded.

    Args:
        query (str): Query sequence.
        target (str): Target sequence.

    Returns:
        float: Sequence identity between two sequences in an alignment.
    """
    # take the shorter sequence
    if len(query) > len(target):
        query, target = target, query

    # count matches excluding gaps
    matches = 0
    for query_char, target_char in zip(query, target):
        if query_char == "-" or target_char == "-":
            continue
        elif query_char == target_char:
            matches += 1

    # calculate identity
    seq_identity = matches / len(query)
    return seq_identity


# skipped in tests
def align(query: str, target: str, matrix: parasail.bindings_v2.Matrix, gap_open: int,
          gap_extend: int) -> parasail.bindings_v2.Result:
    """
    Align two protein sequences using biopython pairwise2.align.globalms.

    Args:
        query_seq (str): Query protein sequence.
        target_seq (str): Target protein sequence.
        matrix (parasail.bindings_v2.Matrix): Substitution matrix.
        gap_open (int): Penalty score for opening a gap.
        gap_continuation (int): Penalty score for gap continuation.

    Returns:
        parasail.bindings_v2.Result: Outputs Result object with required info.
    """

    alignment = parasail.nw_trace(query, target, gap_open, gap_extend, matrix)

    return alignment


def search_alignments(query_seqs: dict, mmseqs_search_output: pd.DataFrame, target_seqs: SeqFileLoader,
                      output_path: pathlib.Path, matrix: str, alignment_gap_open: float, alignment_gap_extend: float,
                      alignment_min_identity: float, threads: int):

    # format of output JSON file:
    # alignments = dict[query_id]
    #     "target_id": target_id,
    #     "sequence_identity" : alignment_sequence_identity(alignment),
    #     "query_sequence": query protein sequence,
    #     "target_sequence": target protein sequence,
    #     "aln_score": alignment score

    alignment_output_json_path = output_path / ALIGNMENTS

    query_seqs_keys = list(query_seqs.keys())

    filtered = mmseqs_search_output[mmseqs_search_output['query'].isin(query_seqs_keys)]
    logging.info("Filtered %i mmseqs matches", len(mmseqs_search_output) - len(filtered))
    logging.info("Total alignments to check %i", len(filtered))

    queries = list(map(lambda x: query_seqs[x], filtered["query"]))
    targets = list(map(lambda x: target_seqs[x], filtered["target"]))

    # parametrize alignment
    parametrized_align = partial(align,
                                 matrix=substitution_matrices[matrix],
                                 gap_open=alignment_gap_open,
                                 gap_extend=alignment_gap_extend)

    # align with parasail
    with ThreadPool(threads) as pool:
        all_alignments = pool.starmap(parametrized_align, zip(queries, targets))

    alignments_output = dict()
    for i, alignment in enumerate(all_alignments):
        # filter out bad alignments based on alignment sequences identity
        query_sequence = alignment.traceback.query
        target_sequence = alignment.traceback.ref

        sequence_identity = alignment_sequences_identity(query_sequence, target_sequence)

        if sequence_identity > alignment_min_identity:
            query_id = filtered["query"].iloc[i]
            target_id = filtered["target"].iloc[i]

            # if query and it's alignment is not in the list yet - add it to list
            if query_id not in alignments_output:
                alignments_output[query_id] = {
                    "target_id": target_id,
                    "query_sequence": query_sequence,
                    "target_sequence": target_sequence,
                    "sequence_identity": sequence_identity,
                    "aln_score": alignment.score
                }

            # select best alignment based on pairwise2.align.globalms.score
            elif alignment.score > alignments_output[query_id]["aln_score"]:
                alignments_output[query_id] = {
                    "target_id": target_id,
                    "query_sequence": query_sequence,
                    "target_sequence": target_sequence,
                    "sequence_identity": sequence_identity,
                    "aln_score": alignment.score
                }

    json.dump(alignments_output, open(alignment_output_json_path, "w", encoding="utf-8"), indent=4, sort_keys=True)
    return alignments_output
