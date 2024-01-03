import logging
from functools import partial
from multiprocessing.pool import ThreadPool

import numpy as np
import pyopal

from mDeepFRI.bio_utils import (AlignmentResult, load_query_sequences,
                                retrieve_fasta_entries_as_dict)
from mDeepFRI.mmseqs import (filter_mmseqs_results, run_mmseqs_search,
                             validate_mmseqs_database)

logger = logging.getLogger(__name__)


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


def run_alignment(query_file: str,
                  mmseqs_db: str,
                  sequence_db: str,
                  output_path: str,
                  mmseqs_min_bit_score: int = None,
                  mmseqs_max_eval: float = 10e-5,
                  mmseqs_min_identity: float = 0.3,
                  top_k: int = 30,
                  alignment_gap_open: int = 10,
                  alignment_gap_continuation: int = 1,
                  identity_threshold: float = 0.5,
                  threads: int = 1):

    query_seqs = load_query_sequences(query_file, output_path)
    mmseqs_valid = validate_mmseqs_database(mmseqs_db)
    if not mmseqs_valid:
        raise ValueError(f"MMSeqs2 database not found in {mmseqs_db}.")

    mmseqs_results = run_mmseqs_search(query_file, mmseqs_db, output_path,
                                       threads)
    filtered_mmseqs_results = filter_mmseqs_results(mmseqs_results,
                                                    mmseqs_min_bit_score,
                                                    mmseqs_max_eval,
                                                    mmseqs_min_identity,
                                                    k_best_hits=top_k,
                                                    threads=threads)

    # if MMSeqs2 alignment is empty
    if filtered_mmseqs_results is None:
        pass
    else:
        target_ids = np.unique(filtered_mmseqs_results["target"]).tolist()
        target_seqs = retrieve_fasta_entries_as_dict(sequence_db, target_ids)

        alignments = align_query(query_seqs, target_seqs, alignment_gap_open,
                                 alignment_gap_continuation,
                                 identity_threshold, threads)

    return alignments
