import logging
from dataclasses import dataclass
from functools import partial
from multiprocessing.pool import ThreadPool

import pyopal

from mDeepFRI.alignment_utils import alignment_sequences_identity


@dataclass
class AlignmentResult:
    query_name: str
    query_sequence: str
    target_name: str
    alignment: str
    identity: float


def align_best_score(query_item, database, target_names, gap_open, gap_extend):
    name, sequence = query_item
    results = database.search(sequence,
                              mode="score",
                              algorithm="nw",
                              gap_open=gap_open,
                              gap_extend=gap_extend)
    best = max(results, key=lambda x: x.score)
    best_sequence = database[best.target_index]
    best_name = target_names[best.target_index]
    single_seq = pyopal.Database([sequence])
    alignment = single_seq.search(best_sequence,
                                  mode="full",
                                  algorithm="nw",
                                  gap_open=gap_open,
                                  gap_extend=gap_extend)[0].alignment
    identity = alignment_sequences_identity(alignment.encode())

    if identity < 0.3:
        aln_results = None
    else:
        aln_results = AlignmentResult(name, sequence, best_name, alignment,
                                      identity)

    return aln_results


def align_query(query_seqs: dict, target_seqs: dict, alignment_gap_open: float,
                alignment_gap_extend: float, threads: int):

    logging.info("Pairwise alignment started.")
    logging.info("Aligning %i queries against %i database sequences.",
                 len(query_seqs), len(target_seqs))
    # Default substitution matrix BLOSUM50
    database = pyopal.Database(list(target_seqs.values()))
    align_against_db = partial(align_best_score,
                               database=database,
                               target_names=list(target_seqs.keys()),
                               gap_open=alignment_gap_open,
                               gap_extend=alignment_gap_extend)

    with ThreadPool(threads) as pool:
        all_alignments = pool.map(align_against_db, query_seqs.items())

    all_alignments = list(filter(None, all_alignments))
    logging.info("Found %i alignments.", len(all_alignments))
    logging.info("Pairwise alignment finished.")

    return all_alignments
