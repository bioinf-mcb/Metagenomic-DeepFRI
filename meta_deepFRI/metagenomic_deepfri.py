import json
import os.path
import pathlib
import logging

from typing import Tuple

from pysam.libcfaidx import FastxFile

from meta_deepFRI.config.names import ATOMS
from meta_deepFRI.DeepFRI.deepfrier import Predictor

from meta_deepFRI.CPP_lib import libAtomDistanceIO
from meta_deepFRI.config.names import SEQ_ATOMS_DATASET_PATH, TARGET_MMSEQS_DB_NAME

from meta_deepFRI.utils.fasta_file_io import SeqFileLoader
from meta_deepFRI.utils.utils import load_deepfri_config
from meta_deepFRI.utils.search_alignments import search_alignments
from meta_deepFRI.utils.mmseqs import run_mmseqs_search

logging.basicConfig(level=logging.DEBUG,
                    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)

###########################################################################
# in a nutshell:
#
#   load_and_verify_job_data
#   1.  select first .faa file inside task_path
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


def check_inputs(query_file: pathlib.Path, database: pathlib.Path,
                 output_path: pathlib.Path) -> Tuple[pathlib.Path, dict, pathlib.Path, SeqFileLoader]:
    """
    Check if input files and directories exist and are valid. Filters out proteins that are too long
    or too short.

    Args:
        query_file (pathlib.Path): Path to a query file with protein sequences.
        database (pathlib.Path): Path to a directory with a pre-built database.
        output_path (pathlib.Path): Path to a directory where results will be saved.

    Raises:
        ValueError: Query file does not contain parsable protein sequences.
        FileNotFoundError: MMSeqs2 database appears to be corrupted.

    Returns:
        Tuple[pathlib.Path, dict, pathlib.Path, SeqFileLoader]: Tuple of query file path, query sequences,
            path to target MMSeqs2 database and target sequences.
    """

    MIN_PROTEIN_LENGTH = 60

    with FastxFile(query_file) as fasta:
        query_seqs = {" ".join([record.name, record.comment]): record.sequence for record in fasta}

    if len(query_seqs) == 0:
        raise ValueError(f"{query_file} does not contain parsable protein sequences.")

    logging.info("Found total of %i protein sequences in %s", len(query_seqs), query_file)

    with open(database / "db_params.json", "r", encoding="utf-8") as f:
        MAX_PROTEIN_LENGTH = json.load(f)["MAX_PROTEIN_LENGTH"]

    prot_len_outliers = {}
    for prot_id, sequence in query_seqs.items():
        prot_len = len(sequence)
        if prot_len > MAX_PROTEIN_LENGTH or prot_len < MIN_PROTEIN_LENGTH:
            prot_len_outliers[prot_id] = prot_len

    for outlier in prot_len_outliers.keys():
        query_seqs.pop(outlier)

    if len(prot_len_outliers) > 0:
        logging.info("Skipping %i proteins due to sequence length outside range %i-%i aa.", len(prot_len_outliers),
                     MIN_PROTEIN_LENGTH, MAX_PROTEIN_LENGTH)
        logging.info("Skipped protein ids will be saved in " \
                     "metadata_skipped_ids_length.json")
        json.dump(prot_len_outliers,
                  open(output_path / 'metadata_skipped_ids_due_to_length.json', "w", encoding="utf-8"),
                  indent=4,
                  sort_keys=True)
        if len(query_seqs) == 0:
            logging.info("All sequences in %s were too long. No sequences will be processed.", query_file)

    # Verify all the files for MMSeqs2 database
    mmseqs2_ext = [
        ".index", ".dbtype", "_h", "_h.index", "_h.dbtype", ".idx", ".idx.index", ".idx.dbtype", ".lookup", ".source"
    ]

    if os.path.isfile(database / TARGET_MMSEQS_DB_NAME):
        target_db = pathlib.Path(database / TARGET_MMSEQS_DB_NAME)
        for ext in mmseqs2_ext:
            assert os.path.isfile(f"{target_db}{ext}")
    else:
        raise FileNotFoundError("MMSeqs2 database appears to be corrupted. Please, rebuild it.")

    target_seqs = SeqFileLoader(database / SEQ_ATOMS_DATASET_PATH)

    return query_file, query_seqs, target_db, target_seqs


## TODO: loading of weights and db as a user-provided parameter
def metagenomic_deepfri(query_file: pathlib.Path, database: pathlib.Path, model_config_json: pathlib.Path,
                        output_path: pathlib.Path, task_path: pathlib.Path):

    query_file, query_seqs, target_db, target_seqs = check_inputs(query_file, database, output_path)

    logging.info("Running metagenomic-DeepFRI for %i sequences", len(query_seqs))
    logging.info("Running MMSeqs2 search for the query against database")
    mmseqs_search_output = run_mmseqs_search(query_file, target_db, output_path)
    mmseqs2_found_proteins = mmseqs_search_output.shape[0]
    logging.info("Found %i proteins in the database", mmseqs2_found_proteins)

    # format: alignments[query_id] = {target_id, identity, alignment[seqA = query_seq, seqB = target_seq, score, start, end]}
    # TODO: remove job_path and job_config from alignments
    alignments = search_alignments(query_seqs, mmseqs_search_output, target_seqs, task_path)
    unaligned_queries = query_seqs.keys() - alignments.keys()
    # timer.log("alignments")
    logging.info(
        "Type of model used depends on whether similar structure was found in the database. If yes, GCN will be used. GCN tends to give better predictions over CNN."
    )
    if len(alignments) > 0:
        logging.info("Using GCN for %i proteins", len(alignments))
    if len(unaligned_queries) > 0:
        logging.info("Using CNN for %i proteins", len(unaligned_queries))
    gcn_cnn_count = {"GCN": len(alignments), "CNN": len(unaligned_queries)}
    json.dump(gcn_cnn_count, open(output_path / "metadata_cnn_gcn_counts.json", "w"), indent=4)

    libAtomDistanceIO.initialize()
    deepfri_models_config = load_deepfri_config(model_config_json)
    target_db_name = job_config.target_db_name

    # DEEPFRI_PROCESSING_MODES = ['mf', 'bp', 'cc', 'ec']
    # mf = molecular_function
    # bp = biological_process
    # cc = cellular_component
    # ec = enzyme_commission
    for mode in job_config.DEEPFRI_PROCESSING_MODES:
        logging.info("Processing mode: %s" % mode)
        # GCN for queries with aligned contact map
        if len(alignments) > 0:
            output_file_name = job_path / f"results_gcn_{mode}"
            if output_file_name.with_suffix('.csv').exists():
                print(f"{output_file_name} results already exists.")
            else:
                gcn_params = deepfri_models_config["gcn"]["models"][mode]
                gcn = Predictor.Predictor(gcn_params, gcn=True)
                for query_id in alignments.keys():
                    alignment = alignments[query_id]
                    query_seq = query_seqs[query_id]
                    target_id = alignment["target_id"]

                    generated_query_contact_map = libAtomDistanceIO.load_aligned_contact_map(
                        str(fsc.SEQ_ATOMS_DATASET_PATH / target_db_name / ATOMS / (target_id + ".bin")),
                        job_config.ANGSTROM_CONTACT_THRESHOLD,
                        alignment["alignment"][0],  # query alignment
                        alignment["alignment"][1],  # target alignment
                        job_config.GENERATE_CONTACTS)

                    gcn.predict_with_cmap(query_seq, generated_query_contact_map, query_id)

                gcn.export_csv(output_file_name.with_suffix('.csv'))
                gcn.export_tsv(output_file_name.with_suffix('.tsv'))
                gcn.export_json(output_file_name.with_suffix('.json'))
                del gcn

        # CNN for queries without satisfying alignments
        if len(unaligned_queries) > 0:
            output_file_name = job_path / f"results_cnn_{mode}"
            if output_file_name.with_suffix('.csv').exists():
                print(f"{output_file_name.name} already exists.")
            else:
                cnn_params = deepfri_models_config["cnn"]["models"][mode]
                cnn = Predictor.Predictor(cnn_params, gcn=False)
                for query_id in unaligned_queries:
                    cnn.predict_from_sequence(query_seqs[query_id], query_id)

                cnn.export_csv(output_file_name.with_suffix('.csv'))
                cnn.export_tsv(output_file_name.with_suffix('.tsv'))
                cnn.export_json(output_file_name.with_suffix('.json'))
                del cnn
