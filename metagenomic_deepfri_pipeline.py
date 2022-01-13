import json

from Bio import SeqIO
from DeepFRI.deepfrier import Predictor

from CONFIG.FOLDER_STRUCTURE import SEQ_ATOMS_DATASET_PATH, DEEPFRI_MODEL_WEIGHTS_JSON_PATH, ATOMS

from CPP_lib.libAtomDistanceIO import initialize as initialize_cpp_lib
from CPP_lib.libAtomDistanceIO import load_aligned_contact_map

from utils.elapsed_time_handler import ElapsedTimeHandler
from utils.run_mmseqs_search import run_mmseqs_search
from utils.search_alignments import search_alignments
from utils.seq_file_loader import SeqFileLoader


def metagenomic_deepfri_pipeline(query_file, target_db, work_path, contact_threshold, generated_contact):
    elapsed_time_handler = ElapsedTimeHandler(work_path / "metadata_runtime.txt")

    mmseqs_search_output = run_mmseqs_search(query_file, target_db, work_path)
    elapsed_time_handler.log("mmseqs2")

    with open(query_file, "r") as f:
        query_seqs = {record.id: record.seq for record in SeqIO.parse(f, "fasta")}
    target_seqs = SeqFileLoader(SEQ_ATOMS_DATASET_PATH)

    # format: alignments[query_id] = {target_id, identity, alignment[seqA = query_seq, seqB = target_seq, score, start, end]}
    elapsed_time_handler.reset()
    alignments = search_alignments(query_seqs, mmseqs_search_output, target_seqs, work_path)
    elapsed_time_handler.log("alignments")
    unaligned_queries = query_seqs.keys() - alignments.keys()

    with open(DEEPFRI_MODEL_WEIGHTS_JSON_PATH) as json_file:
        json_string = json_file.read().replace("./trained_models", f"{DEEPFRI_MODEL_WEIGHTS_JSON_PATH.parent}")
        models_config = json.loads(json_string)

    if len(alignments) > 0:
        print(f"Using GCN for {len(alignments)} proteins")
    if len(unaligned_queries) > 0:
        print(f"Using CNN for {len(unaligned_queries)} proteins")
    initialize_cpp_lib()
    elapsed_time_handler.reset()

    # mf = molecular_function
    # bp = biological_process
    # cc = cellular_component
    # ec = enzyme_commission
    # ['mf', 'bp', 'cc', 'ec']
    for mode in ['mf', 'bp', 'cc', 'ec']:
        print("Processing mode: ", mode)
        if len(alignments) > 0:
            gcn_params = models_config["gcn"]["models"][mode]
            gcn = Predictor.Predictor(gcn_params, gcn=True)
            for query_id in alignments.keys():
                alignment = alignments[query_id]
                query_seq = query_seqs[query_id]
                target_id = alignment["target_id"]

                query_contact_map = load_aligned_contact_map(str(SEQ_ATOMS_DATASET_PATH / ATOMS / (target_id + ".bin")),
                                                             contact_threshold,
                                                             alignment["alignment"].seqA,
                                                             alignment["alignment"].seqB,
                                                             generated_contact)
                gcn.predict_with_cmap(query_seq, query_contact_map, query_id)

            gcn.export_csv(work_path / f"results_gcn_{mode}.csv", verbose=False)
            del gcn
            elapsed_time_handler.log(f"DeepFri_GCN_{mode}")

        if len(unaligned_queries) > 0:
            cnn_params = models_config["cnn"]["models"][mode]
            cnn = Predictor.Predictor(cnn_params, gcn=False)
            for query_id in unaligned_queries:
                cnn.predict_from_sequence(query_seqs[query_id], query_id)

            cnn.export_csv(work_path / f"results_cnn_{mode}.csv", verbose=False)
            del cnn
            elapsed_time_handler.log(f"DeepFri_CNN_{mode}")
