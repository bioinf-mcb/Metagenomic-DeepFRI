import json
import pandas as pd
import pathlib
import pathos

from Bio import pairwise2

from CONFIG.RUNTIME_PARAMETERS import *
from utils.seq_file_loader import SeqFileLoader


def alignment_sequences_identity(alignment):
    matches = [alignment.seqA[i] == alignment.seqB[i] for i in range(len(alignment.seqA))]
    seq_id = (100 * sum(matches)) / len(alignment.seqA)
    return seq_id


def align(query_seq, target_seq):
    return pairwise2.align.globalms(query_seq, target_seq,
                                    PAIRWISE_ALIGNMENT_IDENTICAL,
                                    PAIRWISE_ALIGNMENT_MISSMATCH,
                                    PAIRWISE_ALIGNMENT_GAP_OPEN,
                                    PAIRWISE_ALIGNMENT_GAP_EXTENDING,
                                    one_alignment_only=True)[0]


def search_alignments(query_seqs: dict, mmseqs_search_output: pd.DataFrame, target_seqs: SeqFileLoader, work_path: pathlib.Path):
    # format of output JSON file:
    # alignments = dict[query_id]
    #     "target_id": target_id,
    #     "sequence_identity" : alignment_sequence_identity(alignment)
    #     "alignment": alignment = biopython.alignment
    #         seqA = query_sequence
    #         seqB = target_sequence
    #         score = biopython alignment score
    #         start and end of alignment

    # todo test if loading data from json works
    json_file = work_path / "alignments.json"
    if json_file.exists():
        return json.load(open(json_file, "r"))

    filtered_mmseqs_search = mmseqs_search_output[mmseqs_search_output['bit_score'] > MMSEQS_MIN_BIT_SCORE]
    filtered_mmseqs_search = filtered_mmseqs_search[filtered_mmseqs_search['e_value'] < MMSEQS_MAX_EVAL]

    queries = list(map(lambda x: query_seqs[x], filtered_mmseqs_search["query"]))
    targets = list(map(lambda x: target_seqs[x], filtered_mmseqs_search["target"]))

    with pathos.multiprocessing.ProcessingPool(processes=CPU_COUNT) as p:
        alignments = p.map(align, queries, targets)

    query_alignments = dict()
    for i in range(len(queries)):
        alignment = alignments[i]
        query_id = filtered_mmseqs_search["query"].iloc[i]
        target_id = filtered_mmseqs_search["target"].iloc[i]
        sequence_identity = alignment_sequences_identity(alignment)
        if sequence_identity > ALIGNMENT_MIN_SEQUENCE_IDENTITY:
            if query_id not in query_alignments.keys():
                query_alignments[query_id] = {"target_id": target_id, "alignment": alignment,
                                              "sequence_identity": sequence_identity}
                continue

            if alignment.score > query_alignments[query_id]["alignment"].score:
                query_alignments[query_id] = {"target_id": target_id, "alignment": alignment,
                                              "sequence_identity": sequence_identity}

    json.dump(query_alignments, open(json_file, "w"), indent=4, sort_keys=True)
    return query_alignments
