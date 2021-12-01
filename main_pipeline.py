import json
import shutil
import time

from Bio import SeqIO
from DeepFRI.deepfrier import Predictor

from CONFIG import *
from CPP_lib.libAtomDistanceIO import load_contact_map
from CPP_lib.libAtomDistanceIO import initialize as initialize_cpp_lib
from utils.run_mmseqs_search import run_mmseqs_search
from utils.search_alignments import search_alignments
# chromwell_process_fasta.py looks like the type of script i need to create


def main_pipeline():
    initialize_cpp_lib()

    work_start = str(time.time())
    work_path = (WORK_PATH / work_start)
    while work_path.exists():
        work_start = str(time.time())
        work_path = (WORK_PATH / work_start)
    work_path.mkdir()

    query_faa_files = list(QUERY_PATH.glob("**/*.faa"))
    with open(work_path / 'merged_query_sequences.faa', 'wb') as writer:
        for seq_file in query_faa_files:
            with open(seq_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)

    mmseqs_search_output = run_mmseqs_search(work_path / 'merged_query_sequences.faa', work_path)

    with open(work_path / 'merged_query_sequences.faa', "r") as f:
        query_seqs = {record.id: record.seq for record in SeqIO.parse(f, "fasta")}
    with open(ATOMS_DATASET_PATH / "merged_sequences.faa", "r") as f:
        target_seqs = {record.id: record.seq for record in SeqIO.parse(f, "fasta")}

    alignments = search_alignments(query_seqs, mmseqs_search_output, target_seqs)
    unaligned_queries = query_seqs.keys() - alignments.keys()

    # todo add this path to config?
    with open("/data/trained_models/model_config.json") as json_file:
        params = json.load(json_file)
    gcn_params = params["gcn"]["models"]
    gcn = Predictor.Predictor(gcn_params["mf"], gcn=True)

    for query_id in alignments.keys():
        alignment = alignments[query_id]
        query_seq = query_seqs[query_id]
        target_id = alignment["target_id"]
        target_contact_map = load_contact_map(str(ATOMS_DATASET_PATH / "positions" / (target_id + ".bin")), ANGSTROM_CONTACT_THRESHOLD)

        # todo align contact map

        gcn.predict_with_cmap(query_seq, target_contact_map, query_id)

    for query_id in unaligned_queries:
        print(f"unable to align {query_id} sequence or contact map. Skipping")

    gcn.export_csv("something.csv", True)
    gcn.save_predictions("something.json")
    shutil.rmtree(work_path)


if __name__ == '__main__':
    main_pipeline()
