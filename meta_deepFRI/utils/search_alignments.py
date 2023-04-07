import json
import pandas as pd
import pathlib
import pathos

from Bio import pairwise2

from meta_deepFRI.config.names import ALIGNMENTS
from meta_deepFRI.config.job_config import JobConfig
from meta_deepFRI.config import CPU_COUNT
from meta_deepFRI.utils.fasta_file_io import SeqFileLoader


# alignment sequence identity return value between 0 and 1
def alignment_sequences_identity(alignment):
    matches = [alignment.seqA[i] == alignment.seqB[i] for i in range(len(alignment.seqA))]
    seq_id = sum(matches) / len(alignment.seqA)
    return seq_id


# skipped in testing as wrapper
def align(query_seq, target_seq, match, missmatch, gap_open, gap_continuation):
    return pairwise2.align.globalms(query_seq,
                                    target_seq,
                                    match,
                                    missmatch,
                                    gap_open,
                                    gap_continuation,
                                    one_alignment_only=True)[0]


def search_alignments(query_seqs: dict, mmseqs_search_output: pd.DataFrame, target_seqs: SeqFileLoader,
                      task_path: pathlib.Path, job_config: JobConfig):
    """

    :param query_seqs:
    :param mmseqs_search_output:
    :param target_seqs:
    :param task_path:
    :param job_config:
    :return:
    """
    # format of output JSON file:
    # alignments = dict[query_id]
    #     "target_id": target_id,
    #     "sequence_identity" : alignment_sequence_identity(alignment)
    #     "alignment": alignment = biopython.alignment
    #         0. seqA = query_sequence
    #         1. seqB = target_sequence
    #         2. score = biopython alignment score
    #         3. start and end of alignment

    alignment_output_json_path = task_path / ALIGNMENTS
    if alignment_output_json_path.exists():
        return json.load(open(alignment_output_json_path, "r"))

    print(f"MMseqs search output is {len(mmseqs_search_output)} long.")
    query_seqs_keys = list(query_seqs.keys())
    filtered = mmseqs_search_output[mmseqs_search_output['bit_score'] > job_config.MMSEQS_MIN_BIT_SCORE]
    filtered = filtered[filtered['e_value'] < job_config.MMSEQS_MAX_EVAL]
    filtered = filtered[filtered['identity'] > job_config.MMSEQS_MIN_IDENTITY]
    filtered = filtered[filtered['query'].isin(query_seqs_keys)]
    print(f"Filtered {len(mmseqs_search_output) - len(filtered)} mmseqs matches. "
          f"Total alignments to check {len(filtered)}")

    queries = list(map(lambda x: query_seqs[x], filtered["query"]))
    targets = list(map(lambda x: target_seqs[x], filtered["target"]))

    # Couldn't find more elegant solution on how to use repeating values for pathos.multiprocessing
    # Standard multiprocessing.Pool is out of reach due to problematic pairwise2.align.globalms behaviour.
    # todo make some runtime tests, maybe chunkified sequences will perform better
    match = [job_config.PAIRWISE_ALIGNMENT_MATCH] * len(queries)
    missmatch = [job_config.PAIRWISE_ALIGNMENT_MISSMATCH] * len(queries)
    gap_open = [job_config.PAIRWISE_ALIGNMENT_GAP_OPEN] * len(queries)
    gap_continuation = [job_config.PAIRWISE_ALIGNMENT_GAP_CONTINUATION] * len(queries)

    with pathos.multiprocessing.ProcessingPool(processes=CPU_COUNT) as p:
        all_alignments = p.map(align, queries, targets, match, missmatch, gap_open, gap_continuation)

    alignments_output = dict()
    for i, alignment in enumerate(all_alignments):
        # filter out bad alignments based on alignment sequences identity
        sequence_identity = alignment_sequences_identity(alignment)
        if sequence_identity > job_config.ALIGNMENT_MIN_SEQUENCE_IDENTITY:
            query_id = filtered["query"].iloc[i]
            target_id = filtered["target"].iloc[i]
            if query_id not in alignments_output.keys():
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
