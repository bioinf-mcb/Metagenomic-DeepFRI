import json
import os.path

from DeepFRI.deepfrier import Predictor

from CONFIG.FOLDER_STRUCTURE import SEQ_ATOMS_DATASET_PATH, ATOMS, JOB_CONFIG

from CPP_lib.libAtomDistanceIO import initialize as initialize_cpp_lib
from CPP_lib.libAtomDistanceIO import load_aligned_contact_map

from utils.elapsed_time_logger import ElapsedTimeLogger
from utils.faa_file_io import load_faa_file
from utils.pipeline_utils import find_target_database, load_deepfri_config
from utils.run_mmseqs_search import run_mmseqs_search
from utils.search_alignments import search_alignments
from utils.seq_file_loader import SeqFileLoader

###########################################################################
# in a nutshell:
#
#   load_and_verify_job_data
#   1.  select first .faa file inside job_path
#   2.  filter out proteins that are too long
#   3.  find target database
#
#   metagenomic_deepfri
#   4.  run mmseqs2 search on query and target database
#   5.  find the best alignment for pairs found by mmseqs2 search.
#   6.  If alignment for query exists:
#           DeepFRI GCN for query sequence with aligned target contact map
#       else:
#           DeepFRI CNN for query sequence alone
###########################################################################


def load_and_verify_job_data(job_path, pipeline_config):
    # selects only one .faa file from job_path directory
    query_files = list(job_path.glob("**/*.faa"))
    assert len(query_files) > 0, f"No query .faa files found in {job_path}"
    query_file = query_files[0]
    if len(query_files) > 1:
        print(f"{job_path} contains more than one .faa file. "
              f"Only {query_file} will be processed. {query_files[1:]} will be discarded")

    query_seqs = {record.id: record.seq for record in load_faa_file(query_file)}
    assert len(query_seqs) > 0, f"{query_file} does not contain protein sequences that SeqIO can parse."
    print(f"Found total of {len(query_seqs)} protein sequences in {query_file}")

    # filter out proteins that length is over the CONFIG.RUNTIME_PARAMETERS.MAX_QUERY_CHAIN_LENGTH
    proteins_over_max_length = []
    for query_id in list(query_seqs.keys()):
        if len(query_seqs[query_id]) > pipeline_config["MAX_QUERY_CHAIN_LENGTH"]:
            query_seqs.pop(query_id)
            proteins_over_max_length.append(query_id)

    if len(proteins_over_max_length) > 0:
        print(f"Skipping {len(proteins_over_max_length)} proteins due to sequence length over "
              f"CONFIG.RUNTIME_PARAMETERS.MAX_QUERY_CHAIN_LENGTH. "
              f"\nSkipped protein ids will be saved in metadata_skipped_ids_due_to_max_length.json")
        json.dump(proteins_over_max_length,
                  open(job_path / 'metadata_skipped_ids_due_to_max_length.json', "w"),
                  indent=4,
                  sort_keys=True)
        if len(query_seqs) == 0:
            print(f"All sequences in {query_file} were too long. No sequences will be processed.")

    # select target database
    if os.path.isfile(pipeline_config["target_db"]):
        target_db = pipeline_config["target_db"]
    else:
        # if original database got deleted, try to use the newest database. Its path will be stored in metadata_final_target_db.txt
        print(f"Unable to locate target database {pipeline_config['target_db']}. Looking for the newest one.")
        target_db = find_target_database(pipeline_config["target_db_name"])
        print(f"Found the newest {pipeline_config['target_db_name']} target database. Will be using {target_db}")
        open(job_path / "metadata_final_target_db.txt").write(target_db)
    target_seqs = SeqFileLoader(SEQ_ATOMS_DATASET_PATH / pipeline_config["target_db_name"])

    return query_file, query_seqs, target_db, target_seqs


def metagenomic_deepfri(job_path):
    assert (job_path / JOB_CONFIG).exists(), f"No JOB_CONFIG config file found {job_path / JOB_CONFIG}"
    job_config = json.load(open(job_path / JOB_CONFIG))

    query_file, query_seqs, target_db, target_seqs = load_and_verify_job_data(job_path, job_config)

    if len(query_seqs) == 0:
        print(f"No sequences found. Terminating pipeline.")
        return

    print(f"\nRunning metagenomic_deepfri for {len(query_seqs)} sequences")
    timer = ElapsedTimeLogger(job_path / "metadata_runtime.csv")

    mmseqs_search_output = run_mmseqs_search(query_file, target_db, job_path)
    timer.log("mmseqs2")

    # format: alignments[query_id] = {target_id, identity, alignment[seqA = query_seq, seqB = target_seq, score, start, end]}
    alignments = search_alignments(query_seqs, mmseqs_search_output, target_seqs, job_path, job_config)
    unaligned_queries = query_seqs.keys() - alignments.keys()
    timer.log("alignments")
    if len(alignments) > 0:
        print(f"Using GCN for {len(alignments)} proteins")
    if len(unaligned_queries) > 0:
        print(f"Using CNN for {len(unaligned_queries)} proteins")
    gcn_cnn_count = {"GCN": len(alignments), "CNN": len(unaligned_queries)}
    json.dump(gcn_cnn_count, open(job_path / "metadata_cnn_gcn_counts.json", "w"), indent=4)

    initialize_cpp_lib()
    deepfri_models_config = load_deepfri_config()
    target_db_name = job_config["target_db_name"]

    # DEEPFRI_PROCESSING_MODES = ['mf', 'bp', 'cc', 'ec']
    # mf = molecular_function
    # bp = biological_process
    # cc = cellular_component
    # ec = enzyme_commission
    for mode in job_config["DEEPFRI_PROCESSING_MODES"]:
        timer.reset()
        print("Processing mode: ", mode)
        # GCN for queries with aligned contact map
        if len(alignments) > 0:
            output_file = job_path / f"results_gcn_{mode}.csv"
            if output_file.exists():
                print(f"{output_file.name} already exists.")
            else:
                gcn_params = deepfri_models_config["gcn"]["models"][mode]
                gcn = Predictor.Predictor(gcn_params, gcn=True)
                for query_id in alignments.keys():
                    alignment = alignments[query_id]
                    query_seq = query_seqs[query_id]
                    target_id = alignment["target_id"]

                    generated_query_contact_map = load_aligned_contact_map(
                        str(SEQ_ATOMS_DATASET_PATH / target_db_name / ATOMS / (target_id + ".bin")),
                        job_config["ANGSTROM_CONTACT_THRESHOLD"],
                        alignment["alignment"][0],    # query alignment
                        alignment["alignment"][1],    # target alignment
                        job_config["GENERATE_CONTACTS"])

                    gcn.predict_with_cmap(query_seq, generated_query_contact_map, query_id)

                gcn.export_csv(output_file, verbose=False)
                del gcn
                timer.log(f"deepfri_gcn_{mode}")

        # CNN for queries without satisfying alignments
        if len(unaligned_queries) > 0:
            output_file = job_path / f"results_cnn_{mode}.csv"
            if output_file.exists():
                print(f"{output_file.name} already exists.")
            else:
                cnn_params = deepfri_models_config["cnn"]["models"][mode]
                cnn = Predictor.Predictor(cnn_params, gcn=False)
                for query_id in unaligned_queries:
                    cnn.predict_from_sequence(query_seqs[query_id], query_id)

                cnn.export_csv(output_file, verbose=False)
                del cnn
                timer.log(f"deepfri_cnn_{mode}")

    timer.log_total_time()
