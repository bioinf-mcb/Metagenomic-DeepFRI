import json
import shutil
import time

from Bio import SeqIO
from DeepFRI.deepfrier import Predictor

from CONFIG import *
from CPP_lib.libAtomDistanceIO import initialize as initialize_cpp_lib
from CPP_lib.libAtomDistanceIO import load_aligned_contact_map
from utils.run_mmseqs_search import run_mmseqs_search
from utils.search_alignments import search_alignments
from utils.seq_file_loader import SeqFileLoader
# cromwell_process_fasta.py looks like the type of script i need to create


def main_pipeline():
    initialize_cpp_lib()

    work_start = str(time.time())
    work_path = (WORK_PATH / work_start)
    while work_path.exists():
        work_start = str(time.time())
        work_path = (WORK_PATH / work_start)
    work_path.mkdir()

    query_faa_files = list(QUERY_PATH.glob("**/*.faa"))
    if len(query_faa_files) == 0:
        print("No query protein files found terminating")
        return 0
    with open(work_path / 'merged_query_sequences.faa', 'wb') as writer:
        for seq_file in query_faa_files:
            with open(seq_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)
            # todo remove query_file

    query_file = work_path / 'merged_query_sequences.faa'
    mmseqs_search_output = run_mmseqs_search(query_file, work_path)

    target_seqs = SeqFileLoader(SEQ_ATOMS_DATASET_PATH)
    with open(query_file, "r") as f:
        query_seqs = {record.id: record.seq for record in SeqIO.parse(f, "fasta")}

    # alignments[query_id] = {"target_id": target_id, "alignment": alignment}
    alignments = search_alignments(query_seqs, mmseqs_search_output, target_seqs, work_path)
    unaligned_queries = query_seqs.keys() - alignments.keys()

    if not (DATA_ROOT / "trained_models/model_config.json").exists():
        print("Please run post_setup.py script to download model weights")
        return 1
    with open(DATA_ROOT / "trained_models/model_config.json") as json_file:
        models_config = json.loads(json_file.read().replace("./trained_models/", f"{DATA_ROOT}/trained_models/"))

    if len(alignments) > 0:
        gcn_params = models_config["gcn"]["models"]["mf"]
        gcn = Predictor.Predictor(gcn_params, gcn=True)
        print(f"Using GCN for {len(alignments)} proteins")
        for query_id in alignments.keys():
            alignment = alignments[query_id]
            query_seq = query_seqs[query_id]
            target_id = alignment["target_id"]

            query_contact_map = load_aligned_contact_map(str(SEQ_ATOMS_DATASET_PATH / "positions" / (target_id + ".bin")),
                                                         ANGSTROM_CONTACT_THRESHOLD,
                                                         alignment["alignment"].seqA,
                                                         alignment["alignment"].seqB,
                                                         1)
            gcn.predict_with_cmap(query_seq, query_contact_map, query_id)
        gcn.export_csv(work_path / "result_gcn.csv", verbose=False)
        del gcn
    else:
        print("No aligned contact maps found")

    if len(unaligned_queries) > 0:
        cnn_params = models_config["cnn"]["models"]["mf"]
        cnn = Predictor.Predictor(cnn_params, gcn=False)
        print(f"Using CNN for {len(unaligned_queries)} proteins")
        for query_id in unaligned_queries:
            cnn.predict_from_sequence(query_seqs[query_id], query_id)

        cnn.export_csv(work_path / "result_cnn.csv", verbose=False)
        del cnn

    (FINISHED_PATH / work_start).mkdir(parents=True, exist_ok=True)
    shutil.copy(work_path / "merged_query_sequences.faa", FINISHED_PATH / work_start / "sequences.faa")
    shutil.copy(work_path / "result_gcn.csv", FINISHED_PATH / work_start / "result_gcn.csv")
    shutil.copy(work_path / "result_cnn.csv", FINISHED_PATH / work_start / "result_cnn.csv")

    shutil.rmtree(work_path)


if __name__ == '__main__':
    main_pipeline()
