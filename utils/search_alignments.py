import json
import pandas as pd
import pathlib
import pathos

from Bio import pairwise2

from CONFIG.RUNTIME_PARAMETERS import *
from utils.seq_file_loader import SeqFileLoader


def align(query_seq, target_seq):
    best_alignment = None
    best_score = -1
    alignments = pairwise2.align.globalxx(query_seq, target_seq)
    for alignment in alignments:
        if alignment.score > best_score:
            best_score = alignment.score
            best_alignment = alignment
    # save identity score of output
    if best_alignment is not None:
        return best_alignment
    return None


def search_alignments(query_seqs: dict, mmseqs_search_output: pd.DataFrame, target_seqs: SeqFileLoader, work_path: pathlib.Path):
    json_file = work_path / "alignments.json"
    if json_file.exists():
        return json.load(open(json_file, "r"))

    queries = list(map(lambda x: query_seqs[x], mmseqs_search_output["query"]))
    targets = list(map(lambda x: target_seqs[x], mmseqs_search_output["target"]))

    with pathos.multiprocessing.ProcessingPool(processes=CPU_COUNT) as p:
        alignments = p.map(align, queries, targets)

    query_alignments = dict()
    for i in range(len(queries)):
        alignment = alignments[i]
        query_id = mmseqs_search_output["query"].iloc[i]
        target_id = mmseqs_search_output["target"].iloc[i]

        if query_id not in query_alignments.keys():
            query_alignments[query_id] = {"target_id": target_id, "alignment": alignment}
            continue

        if alignment.score > query_alignments[query_id]["alignment"].score:
            query_alignments[query_id] = {"target_id": target_id, "alignment": alignment}

    json.dump(query_alignments, open(json_file, "w"))
    return query_alignments
