import json
import pandas as pd
import pathlib
import pathos

import logging

from Bio import pairwise2

from meta_deepFRI.config.names import ALIGNMENTS
from meta_deepFRI.utils.fasta_file_io import SeqFileLoader


# alignment sequence identity return value between 0 and 1
def alignment_sequences_identity(alignment):
    """
    Calculate sequence identity between two sequences in an alignment.

    Args:
        alignment (Bio.pairwise2.Alignment): Alignment object.

    Returns:
        float: Sequence identity between two sequences in an alignment.
    """
    matches = [alignment.seqA[i] == alignment.seqB[i] for i in range(len(alignment.seqA))]
    seq_id = sum(matches) / len(alignment.seqA)
    return seq_id


# skipped in testing as wrapper
def align(query_seq: str, target_seq: str, match: float, missmatch: float, gap_open: float, gap_continuation: float):
    """
    Align two protein sequences using biopython pairwise2.align.globalms.

    Args:
        query_seq (str): Query protein sequence.
        target_seq (str): Target protein sequence.
        match (float): Match score.
        missmatch (float): Missmatch score.
        gap_open (float): Score for opening a gap.
        gap_continuation (float): Score for gap continuation.

    Returns:
        Bio.pairwise2.Alignment: problematic output type. Outputs Alignment object with required info.
    """
    return pairwise2.align.globalms(query_seq,
                                    target_seq,
                                    match,
                                    missmatch,
                                    gap_open,
                                    gap_continuation,
                                    one_alignment_only=True)[0]


def search_alignments(query_seqs: dict, mmseqs_search_output: pd.DataFrame, target_seqs: SeqFileLoader,
                      output_path: pathlib.Path, mmseqs_min_bit_score: float, mmseqs_max_eval: float,
                      mmseqs_min_identity: float, alignment_match: float, alignment_missmatch: float,
                      alignment_gap_open: float, alignment_gap_continuation: float, alignment_min_identity: float,
                      threads):

    # format of output JSON file:
    # alignments = dict[query_id]
    #     "target_id": target_id,
    #     "sequence_identity" : alignment_sequence_identity(alignment)
    #     "alignment": alignment = biopython.alignment
    #         0. seqA = query_sequence
    #         1. seqB = target_sequence
    #         2. score = biopython alignment score
    #         3. start and end of alignment

    alignment_output_json_path = output_path / ALIGNMENTS

    query_seqs_keys = list(query_seqs.keys())
    # MMSeqs2 alginment filters
    if mmseqs_min_identity:
        mmseqs_search_output = mmseqs_search_output.query(f"identity >= {mmseqs_min_identity}")
    if mmseqs_min_bit_score:
        mmseqs_search_output = mmseqs_search_output.query(f"bitscore >= {mmseqs_min_bit_score}")
    if mmseqs_max_eval:
        mmseqs_search_output = mmseqs_search_output.query(f"e_value <= {mmseqs_max_eval}")

    filtered = mmseqs_search_output[mmseqs_search_output['query'].isin(query_seqs_keys)]
    logging.info("Filtered %i mmseqs matches", len(mmseqs_search_output) - len(filtered))
    logging.info("Total alignments to check %i", len(filtered))

    queries = list(map(lambda x: query_seqs[x], filtered["query"]))
    targets = list(map(lambda x: target_seqs[x], filtered["target"]))

    # Couldn't find more elegant solution on how to use repeating values for pathos.multiprocessing
    # Standard multiprocessing.Pool is out of reach due to problematic pairwise2.align.globalms behaviour.
    match = [alignment_match] * len(queries)
    missmatch = [alignment_missmatch] * len(queries)
    gap_open = [alignment_gap_open] * len(queries)
    gap_continuation = [alignment_gap_continuation] * len(queries)

    with pathos.multiprocessing.ProcessingPool(processes=threads) as p:
        all_alignments = p.map(align, queries, targets, match, missmatch, gap_open, gap_continuation)

    alignments_output = dict()
    for i, alignment in enumerate(all_alignments):
        # filter out bad alignments based on alignment sequences identity
        sequence_identity = alignment_sequences_identity(alignment)
        if sequence_identity > alignment_min_identity:
            query_id = filtered["query"].iloc[i]
            target_id = filtered["target"].iloc[i]
            if query_id not in alignments_output:
                alignments_output[query_id] = {
                    "target_id": target_id,
                    "alignment": alignment,
                    "sequence_identity": sequence_identity
                }
            # select best alignment based on pairwise2.align.globalms.score
            elif alignment.score > alignments_output[query_id]["alignment"].score:
                alignments_output[query_id] = {
                    "target_id": target_id,
                    "alignment": alignment,
                    "sequence_identity": sequence_identity
                }

    json.dump(alignments_output, open(alignment_output_json_path, "w"), indent=4, sort_keys=True)
    return alignments_output
