import logging
from functools import partial
from multiprocessing.pool import ThreadPool

import pyopal

from mDeepFRI.alignment_utils import alignment_identity

logger = logging.getLogger(__name__)


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
    """
    Class for storing pairwise alignment results.

    Attributes:
        query_name (str): Name of the query sequence.
        query_sequence (str): Query sequence.
        target_name (str): Name of the target sequence.
        target_sequence (str): Target sequence.
        alignment (str): Alignment string.
        gapped_sequence (str): Query sequence with gaps.
        gapped_target (str): Target sequence with gaps.
        identity (float): Identity between two sequences.

    Methods:
        insert_gaps: Inserts gaps into query and target sequences.
        calculate_identity: Calculates identity between query and target.
    """
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
        self.db_name = None

    def __str__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"identity={self.identity})"

    def __repr__(self):
        return f"AlignmentResult(query_name={self.query_name}, target_name={self.target_name}, " \
               f"identity={self.identity})"

    def insert_gaps(self):
        """
        Inserts gaps into query and target sequences.

        Args:
            None
        Returns:
            self
        """
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
        """
        Calculates identity between query and target.

        Args:
            None
        Returns:
            self
        """

        self.identity = alignment_identity(self.gapped_sequence,
                                           self.gapped_target)

        return self


def align_best_score(query_item: str, database: pyopal.Database,
                     target_names: list, gap_open: int, gap_extend: int,
                     identity_threshold: float):

    query_name, query_sequence = query_item

    logger.debug("Aligning %s", query_name)
    # current bug in pyOpal
    # https://github.com/althonos/pyopal/issues/3
    try:
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
        logger.debug("Alignment of %s finished.", query_name)
        if result.identity >= identity_threshold:
            pass
        else:
            logger.debug("No alignment for  %s.", query_name)
            result = None
    except RuntimeError:
        logger.debug("Alignment of %s failed.", query_name)
        result = None

    return result


def align_query(query_seqs: dict,
                target_seqs: dict,
                alignment_gap_open: float,
                alignment_gap_extend: float,
                identity_threshold: str = 0.3,
                threads: int = 1):

    logger.info("Pairwise alignment started.")
    logger.info("Aligning %i queries against %i database sequences.",
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
    logger.info("Found %i alignments.", len(filtered_alignments))
    logger.info("Pairwise alignment finished.")

    return filtered_alignments
