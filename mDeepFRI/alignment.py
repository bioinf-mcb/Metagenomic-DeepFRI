import logging
from functools import partial
from multiprocessing.pool import ThreadPool

import pyopal

from mDeepFRI.alignment_utils import alignment_identity


def insert_gaps(sequence, reference, alignment_string):
    sequence = list(sequence)
    reference = list(reference)
    alignment_string = list(alignment_string)

    for i, a in enumerate(alignment_string):
        if a == "I":
            sequence.insert(i, "-")
        elif a == "D":
            reference.insert(i, "-")
    return "".join(sequence), "".join(reference)


class AlignmentResult:
    def __init__(self,
                 query_name,
                 query_sequence,
                 target_name,
                 target_sequence,
                 alignment,
                 gapped_sequence="",
                 gapped_target="",
                 identity=None):
        self.query_name = query_name
        self.query_sequence = query_sequence
        self.target_name = target_name
        self.target_sequence = target_sequence
        self.alignment = alignment
        self.gapped_sequence = gapped_sequence
        self.gapped_target = gapped_target
        self.identity = identity

    def __str__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"identity={self.identity})"

    def __repr__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"identity={self.identity})"

    def insert_gaps(self):
        seq = list(self.query_sequence)
        ref = list(self.target_sequence)
        aln = list(self.alignment)

        for i, a in enumerate(aln):
            if a == "I":
                seq.insert(i, "-")
            elif a == "D":
                ref.insert(i, "-")
        self.gapped_sequence = "".join(seq)
        self.gapped_target = "".join(ref)

        return self

    def calculate_identity(self):
        self.identity = alignment_identity(self.gapped_sequence,
                                           self.gapped_target)

        return self


def align_best_score(query_item: str, database: pyopal.Database,
                     target_names: list, gap_open: int, gap_extend: int,
                     identity_threshold: float):

    query_name, query_sequence = query_item
    results = database.search(query_sequence,
                              mode="score",
                              algorithm="nw",
                              gap_open=gap_open,
                              gap_extend=gap_extend)

    best = max(results, key=lambda x: x.score)
    best_sequence = database[best.target_index]
    best_name = target_names[best.target_index]
    single_seq = pyopal.Database([best_sequence])
    alignment = single_seq.search(query_sequence,
                                  mode="full",
                                  algorithm="nw",
                                  gap_open=gap_open,
                                  gap_extend=gap_extend)[0].alignment

    result = AlignmentResult(query_name, query_sequence, best_name,
                             best_sequence, alignment)
    result.insert_gaps().calculate_identity()

    if result.identity >= identity_threshold:
        return result

    return None


def align_query(query_seqs: dict,
                target_seqs: dict,
                alignment_gap_open: float,
                alignment_gap_extend: float,
                identity_threshold: str = 0.3,
                threads: int = 1):

    logging.info("Pairwise alignment started.")
    logging.info("Aligning %i queries against %i database sequences.",
                 len(query_seqs), len(target_seqs))
    # Default substitution matrix BLOSUM50
    database = pyopal.Database(list(target_seqs.values()))
    align_against_db = partial(align_best_score,
                               database=database,
                               target_names=list(target_seqs.keys()),
                               gap_open=alignment_gap_open,
                               gap_extend=alignment_gap_extend,
                               identity_threshold=identity_threshold)

    with ThreadPool(threads) as pool:
        all_alignments = pool.map(align_against_db, query_seqs.items())

    filtered_alignments = list(filter(None, all_alignments))
    logging.info("Found %i alignments.", len(filtered_alignments))
    logging.info("Pairwise alignment finished.")

    return filtered_alignments
