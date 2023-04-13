import json
import os.path
import pathlib
import logging

from typing import Tuple, List

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
        query_seqs = {record.name: record.sequence for record in fasta}

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
def metagenomic_deepfri(query_file: pathlib.Path,
                        database: pathlib.Path,
                        model_config_json: pathlib.Path,
                        output_path: pathlib.Path,
                        output_format: List[str],
                        deepfri_processing_modes: List[str],
                        angstrom_contact_threshold: float = 6,
                        generate_contacts: int = 2,
                        mmseqs_min_bit_score: float = None,
                        mmseqs_max_eval: float = None,
                        mmseqs_min_identity: float = 0.5,
                        alignment_match: float = 2,
                        alignment_missmatch: float = -1,
                        alignment_gap_open: float = -0.5,
                        alignment_gap_continuation: float = -0.1,
                        alignment_min_identity: float = 0.3,
                        threads: int = 1):

    logging.info("Starting metagenomic-DeepFRI.")

    query_file, query_seqs, target_db, target_seqs = check_inputs(query_file, database, output_path)

    logging.info("Running metagenomic-DeepFRI for %i sequences", len(query_seqs))
    logging.info("Running MMSeqs2 search for the query against database")
    mmseqs_search_output = run_mmseqs_search(query_file, target_db, output_path)
    mmseqs2_found_proteins = mmseqs_search_output.shape[0]
    logging.info("Found %i proteins in the database", mmseqs2_found_proteins)

    # format: alignments[query_id] = {target_id, identity, alignment[seqA = query_seq, seqB = target_seq, score, start, end]}
    alignments = search_alignments(query_seqs, mmseqs_search_output, target_seqs, output_path, mmseqs_min_bit_score,
                                   mmseqs_max_eval, mmseqs_min_identity, alignment_match, alignment_missmatch,
                                   alignment_gap_open, alignment_gap_continuation, alignment_min_identity, threads)

    unaligned_queries = query_seqs.keys() - alignments.keys()

    libAtomDistanceIO.initialize()
    deepfri_models_config = load_deepfri_config(model_config_json)

    # deepfri_processing_modes = ['mf', 'bp', 'cc', 'ec']
    # mf = molecular_function
    # bp = biological_process
    # cc = cellular_component
    # ec = enzyme_commission

    for mode in deepfri_processing_modes:
        logging.info("Processing mode: %s", mode)
        # GCN for queries with aligned contact map
        gcn_prots, cnn_prots = len(alignments), len(unaligned_queries)

        if gcn_prots > 0:
            logging.info("Predicting with GCN: %i proteins", gcn_prots)
            output_file_name = output_path / f"results_gcn_{mode}"

            gcn_params = deepfri_models_config["gcn"]["models"][mode]
            gcn = Predictor.Predictor(gcn_params, gcn=True)

            for query_id, alignment in alignments.items():
                logging.info("Predicting %s", query_id)
                query_seq = query_seqs[query_id]
                target_id = alignment["target_id"]

                generated_query_contact_map = libAtomDistanceIO.load_aligned_contact_map(
                    str(database / SEQ_ATOMS_DATASET_PATH / ATOMS / (target_id + ".bin")),
                    angstrom_contact_threshold,
                    alignment["alignment"][0],  # query alignment
                    alignment["alignment"][1],  # target alignment
                    generate_contacts)

                # conversion of cmap to an adequate format
                cmap = generated_query_contact_map.reshape(1, *generated_query_contact_map.shape) * 1

                # running the actual prediction
                gcn.predict_with_cmap(query_seq, cmap, query_id)

                if "tsv" in output_format:
                    gcn.export_tsv(output_file_name.with_suffix('.tsv'))
                if "csv" in output_format:
                    gcn.export_csv(output_file_name.with_suffix('.csv'))
                if "json" in output_format:
                    gcn.export_json(output_file_name.with_suffix('.json'))

            del gcn

        # CNN for queries without satisfying alignments
        if cnn_prots > 0:
            logging.info("Predicting with CNN:  %i proteins", cnn_prots)
            output_file_name = output_path / f"results_cnn_{mode}"

            cnn_params = deepfri_models_config["cnn"]["models"][mode]
            cnn = Predictor.Predictor(cnn_params, gcn=False)
            for query_id in unaligned_queries:
                logging.info("Predicting %s", query_id)
                cnn.predict_from_sequence(query_seqs[query_id], query_id)

            if "tsv" in output_format:
                cnn.export_tsv(output_file_name.with_suffix('.tsv'))
            if "csv" in output_format:
                cnn.export_csv(output_file_name.with_suffix('.csv'))
            if "json" in output_format:
                cnn.export_json(output_file_name.with_suffix('.json'))

            del cnn

    logging.info("meta-DeepFRI finished successfully")
