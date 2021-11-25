import json
import shutil
import time

from Bio import pairwise2
from Bio import SeqIO

from CONFIG import *
from CPP_lib.libAtomDistanceIO import load_contact_map, initialize
from DeepFRI.deepfrier import Predictor
from utils import add_path_to_env
from run_mmseqs_search import run_mmseqs_search
# chromwell_process_fasta.py looks like the type of script i need to create


def main_pipeline():
    initialize()
    pipeline_start = str(time.time())
    faa_files = list(QUERY_PATH.glob("**/*.faa"))

    work_path = (WORK_PATH / pipeline_start)
    while work_path.exists():
        pipeline_start = str(int(time.time()))
        work_path = (WORK_PATH / pipeline_start)
    work_path.mkdir()

    query_data = dict()
    for query_file in faa_files:
        query_data.update({record.id: record.seq for record in SeqIO.parse(query_file, "fasta")})

    with open(work_path / 'merged_sequences.faa', 'wb') as writer:
        for seq_file in faa_files:
            with open(seq_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)

    output_data = run_mmseqs_search(work_path / 'merged_sequences.faa', work_path)

    # shutil.rmtree(work_path)

    with open("/data/trained_models/model_config.json") as json_file:
        params = json.load(json_file)
    gcn_params = params["gcn"]["models"]
    gcn = Predictor.Predictor(gcn_params["mf"], gcn=True)

    with open(ATOMS_DATASET_PATH / "merged_sequences.faa", "r") as f:
        target_sequences = {record.id: record.seq for record in SeqIO.parse(f, "fasta")}

    for i in range(len(output_data)):
        query_id = output_data["query"].iloc[i]
        query_sequence = query_data[query_id]
        target_id = output_data["target"].iloc[i]
        target_sequence = target_sequences[target_id]

        target_cmap = load_contact_map(str(ATOMS_DATASET_PATH / "positions" / (target_id + ".bin")),
                                       ANGSTROM_CONTACT_THRESHOLD)

        alignment = pairwise2.align.globalxx(query_sequence, target_sequence)[0]
        indexses = []
        if alignment.score == len(query_sequence):
            if len(query_sequence) == len(alignment.seqA):
                aligned_cmap = target_cmap
            else:
                for j in range(len(alignment.seqA)):
                    if alignment.seqA[j] != '-':
                        indexses.append(j)
                aligned_cmap = target_cmap[indexses]
                aligned_cmap = aligned_cmap[:, indexses]
        else:
            # print(i)
            # from Bio.pairwise2 import format_alignment
            # print(format_alignment(*alignment))
            continue

        gcn.predict_with_cmap(query_sequence, aligned_cmap, query_id)

    gcn.export_csv("somethinf.csv", True)
    gcn.save_predictions("somethinf.json")


if __name__ == '__main__':
    main_pipeline()
