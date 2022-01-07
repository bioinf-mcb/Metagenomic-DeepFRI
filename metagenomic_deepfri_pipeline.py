import json
import shutil

from Bio import SeqIO
from DeepFRI.deepfrier import Predictor

from CONFIG.FOLDER_STRUCTURE import *
from CONFIG.RUNTIME_PARAMETERS import ANGSTROM_CONTACT_THRESHOLD
from CPP_lib.libAtomDistanceIO import initialize as initialize_cpp_lib
from CPP_lib.libAtomDistanceIO import load_aligned_contact_map
from utils.run_mmseqs_search import run_mmseqs_search
from utils.search_alignments import search_alignments
from utils.seq_file_loader import SeqFileLoader


def metagenomic_deepfri_pipeline(query_file, target_db, work_path):
    initialize_cpp_lib()

    mmseqs_search_output = run_mmseqs_search(query_file, target_db, work_path)

    with open(query_file, "r") as f:
        query_seqs = {record.id: record.seq for record in SeqIO.parse(f, "fasta")}
    target_seqs = SeqFileLoader(SEQ_ATOMS_DATASET_PATH)

    # format: alignments[query_id] = {"target_id": target_id, "alignment": alignment}
    alignments = search_alignments(query_seqs, mmseqs_search_output, target_seqs, work_path)
    unaligned_queries = query_seqs.keys() - alignments.keys()

    with open(DEEPFRI_MODEL_WEIGHTS_JSON_PATH) as json_file:
        models_config = json.loads(json_file.read().replace("./trained_models/", f"{DATA_ROOT}/trained_models/"))

    if len(alignments) > 0:
        print(f"Using GCN for {len(alignments)} proteins")
        gcn_params = models_config["gcn"]["models"]["mf"]
        gcn = Predictor.Predictor(gcn_params, gcn=True)
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
        print(f"Using CNN for {len(unaligned_queries)} proteins")
        cnn_params = models_config["cnn"]["models"]["mf"]
        cnn = Predictor.Predictor(cnn_params, gcn=False)
        for query_id in unaligned_queries:
            cnn.predict_from_sequence(query_seqs[query_id], query_id)

        cnn.export_csv(work_path / "result_cnn.csv", verbose=False)
        del cnn

    print("Finished! Saving output files to ", FINISHED_PATH / work_path.name)
    (FINISHED_PATH / work_path.name).mkdir(parents=True, exist_ok=True)
    shutil.copy(query_file, FINISHED_PATH / work_path.name)
    if (work_path / "result_gcn.csv").exists():
        shutil.copy(work_path / "result_gcn.csv", FINISHED_PATH / work_path.name)
    if (work_path / "result_cnn.csv").exists():
        shutil.copy(work_path / "result_cnn.csv", FINISHED_PATH / work_path.name)

    shutil.rmtree(work_path)
