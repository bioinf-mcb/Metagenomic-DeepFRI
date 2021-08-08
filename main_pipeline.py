import json
import shutil
import time

from Bio import pairwise2
from Bio import SeqIO
import pandas as pd


from CONFIG import *
from CPP_lib.libAtomDistanceIO import load_contact_map, initialize_CPP_LIB
from DeepFRI.deepfrier import Predictor
from utils import run_command

# chromwell_process_fasta.py looks like the type of script i need to create

def run_mmseqs_search(query_file, work_path):
    target_databases = sorted(list(MMSEQS_DATABASES_PATH.iterdir()))[-1]

    run_command(f"mmseqs createdb {query_file} {work_path / 'queryDB'} --dbtype 1")
    run_command(f"mmseqs search {work_path / 'queryDB'} {target_databases / DEFAULT_MMSEQS_NAME} {work_path / 'resultDB'} {TMP_FOLDER_PATH}")
    run_command(f"mmseqs convertalis {work_path / 'queryDB'} {target_databases / DEFAULT_MMSEQS_NAME} {work_path / 'resultDB'} {work_path / 'resultDB.m8'}")

    column_names = ["query", "target", "identity", "alignment_length", "mismatches", "gap_openings", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score"]
    output_data = pd.read_csv(work_path / 'resultDB.m8', sep="\t", names=column_names)
    shutil.rmtree(work_path)

    # work around pandas Automatic exclusion of “nuisance” columns
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#other-useful-features
    output_data["original_index"] = output_data.index
    grouped = output_data.groupby(["query"], sort=False)
    maximums = grouped.max('identity')
    output_data = output_data.loc[maximums["original_index"]]

    output_data = output_data[output_data.identity > 0.95]
    return output_data


def main_pipeline():
    initialize_CPP_LIB()
    pipeline_start = str(time.time())
    faa_files = list(QUERY_FOLDER_PATH.glob("**/*.faa"))

    query_file = faa_files[0]
    query_data = {record.id: record.seq for record in SeqIO.parse(query_file, "fasta")}

    work_path = (WORK_PATH / pipeline_start)
    while work_path.exists():
        pipeline_start = str(int(time.time()))
        work_path = (WORK_PATH / pipeline_start)
    work_path.mkdir()

    output_data = run_mmseqs_search(query_file, work_path)


    with open("/data/trained_models/model_config.json") as json_file:
        params = json.load(json_file)
    gcn_params = params["gcn"]["models"]
    cnn_params = params["cnn"]["models"]

    gcn = Predictor.Predictor(gcn_params["mf"], gcn=True)

    for i in range(len(output_data)):
        query_id = output_data["query"].iloc[i]
        query_sequence = query_data[query_id]
        target_id = output_data["target"].iloc[i]
        len(query_sequence)
        with open(ATOMS_DATASET_PATH / "seq" / (target_id + ".faa")) as f:
            for record in SeqIO.parse(f, "fasta"):
                target_sequence = record.seq

        target_cmap = load_contact_map(str(ATOMS_DATASET_PATH / "positions" / (target_id + ".bin")), ANGSTROM_CONTACT_THRESHOLD)

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
            print(i)
            from Bio.pairwise2 import format_alignment
            print(format_alignment(*alignment))
            continue

        gcn.predict_with_cmap(query_sequence, aligned_cmap, query_id)

    gcn.export_csv("somethinf.csv",True)
    gcn.save_predictions("somethinf.json")


if __name__ == '__main__':
    main_pipeline()
