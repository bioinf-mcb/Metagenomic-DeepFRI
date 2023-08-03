import json
import logging
import os.path
import pathlib
from typing import Dict, List

import numpy as np

from mDeepFRI import MERGED_SEQUENCES, TARGET_MMSEQS_DB_NAME
from mDeepFRI.alignment import align_query
from mDeepFRI.bio_utils import (align_contact_map, calculate_cmap,
                                load_fasta_as_dict,
                                retrieve_fasta_entries_as_dict,
                                retrieve_structure)
from mDeepFRI.mmseqs import filter_mmseqs_results, run_mmseqs_search
from mDeepFRI.predict import Predictor

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def load_query_sequences(query_file, output_path) -> Dict[str, str]:
    """
    Loads query protein sequences from FASTA file. Filters
    out sequences that are too short or too long.

    Args:
        query_file (str): Path to FASTA file with query protein sequences.
        output_path (str): Path to output folder.

    Returns:
        query_seqs (dict): Dictionary with query protein headers to sequences.
    """

    # By DeepFRI design
    MIN_PROTEIN_LENGTH = 60
    MAX_PROTEIN_LENGTH = 1_000

    query_seqs = load_fasta_as_dict(query_file)

    if len(query_seqs) == 0:
        raise ValueError(
            f"{query_file} does not contain parsable protein sequences.")

    logging.info("Found total of %i protein sequences in %s", len(query_seqs),
                 query_file)

    prot_len_outliers = {}
    for prot_id, sequence in query_seqs.items():
        prot_len = len(sequence)
        if prot_len > MAX_PROTEIN_LENGTH or prot_len < MIN_PROTEIN_LENGTH:
            prot_len_outliers[prot_id] = prot_len

    for outlier in prot_len_outliers.keys():
        query_seqs.pop(outlier)

    if len(prot_len_outliers) > 0:
        logging.info(
            "Skipping %i proteins due to sequence length outside range %i-%i aa.",
            len(prot_len_outliers), MIN_PROTEIN_LENGTH, MAX_PROTEIN_LENGTH)
        logging.info("Skipped protein ids will be saved in " \
                     "metadata_skipped_ids_length.json.")
        json.dump(prot_len_outliers,
                  open(output_path / 'metadata_skipped_ids_due_to_length.json',
                       "w",
                       encoding="utf-8"),
                  indent=4,
                  sort_keys=True)

    assert len(query_seqs
               ) > 0, "All proteins were filtered out due to sequence length."

    return query_seqs


def check_mmseqs_database(database: pathlib.Path):
    """
    Check if MMSeqs2 database is intact.

    Args:
        query_file (pathlib.Path): Path to a query file with protein sequences.
        database (pathlib.Path): Path to a directory with a pre-built database.
        output_path (pathlib.Path): Path to a directory where results will be saved.

    Raises:
        FileNotFoundError: MMSeqs2 database appears to be corrupted.

    Returns:
        target_db (pathlib.Path): Path to MMSeqs2 database.
    """

    # Verify all the files for MMSeqs2 database
    mmseqs2_ext = [
        ".index", ".dbtype", "_h", "_h.index", "_h.dbtype", ".idx",
        ".idx.index", ".idx.dbtype", ".lookup", ".source"
    ]

    if os.path.isfile(database / TARGET_MMSEQS_DB_NAME):
        target_db = pathlib.Path(database / TARGET_MMSEQS_DB_NAME)
        for ext in mmseqs2_ext:
            assert os.path.isfile(f"{target_db}{ext}")
    else:
        raise FileNotFoundError(
            "MMSeqs2 database appears to be corrupted. Please, rebuild it.")

    return target_db


def load_deepfri_config(weights: pathlib.Path) -> pathlib.Path:
    """
    Check if DeepFRI weights are valid and load config.
    Args:
        weights :

    Returns:
        pathlib.Path: Path to DeepFRI config.
    """

    assert weights.exists(), f"DeepFRI weights not found at {weights}"
    assert weights.is_dir(
    ), "DeepFRI weights should be a directory, not a file."

    config_path = weights / "model_config.json"
    assert config_path.exists(
    ), "DeepFRI weights are missing model_config.json"

    with open(config_path, "r", encoding="utf-8") as f:
        models_config = json.load(f)

    for net in ["cnn", "gcn"]:
        for model_type, model_path in models_config[net]["models"].items():
            model_name = weights / (pathlib.Path(model_path).name + ".onnx")
            config_name = weights / (pathlib.Path(model_path).name +
                                     "_model_params.json")
            assert model_name.exists(
            ), f"DeepFRI weights are missing {model_type} model at {model_name}"
            assert config_name.exists(
            ), f"DeepFRI weights are missing {model_type} model config at {config_name}"
            # correct config
            models_config[net]["models"][model_type] = str(
                model_name.absolute())
    return models_config


## TODO: structure output folder
def predict_protein_function(
        query_file: str,
        database: str,
        weights: str,
        output_path: str,
        deepfri_processing_modes: List[str] = ["ec", "bp", "mf", "cc"],
        angstrom_contact_threshold: float = 6,
        generate_contacts: int = 2,
        mmseqs_min_bit_score: float = None,
        mmseqs_max_eval: float = 10e-5,
        mmseqs_min_identity: float = 0.3,
        top_k: int = 30,
        alignment_gap_open: float = 10,
        alignment_gap_continuation: float = 1,
        identity_threshold: float = 0.3,
        threads: int = 1):
    """
    Run metagenomic-DeepFRI.
    Args:
        query_file:
        database:
        weights:
        output_path:
        output_format:
        deepfri_processing_modes:
        angstrom_contact_threshold:
        generate_contacts:
        mmseqs_min_bit_score:
        mmseqs_max_eval:
        mmseqs_min_identity:
        alignment_match:
        alignment_missmatch:
        alignment_gap_open:
        alignment_gap_continuation:
        alignment_min_identity:
        threads:

    Returns:

    """
    query_file = pathlib.Path(query_file)
    database = pathlib.Path(database)
    weights = pathlib.Path(weights)
    output_path = pathlib.Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    query_seqs = load_query_sequences(query_file, output_path)
    logging.info("Aligning %i sequences with MMSeqs2.", len(query_seqs))
    target_db = check_mmseqs_database(database)

    mmseqs_results = run_mmseqs_search(query_file, target_db, output_path)
    filtered_mmseqs_results = filter_mmseqs_results(mmseqs_results,
                                                    mmseqs_min_bit_score,
                                                    mmseqs_max_eval,
                                                    mmseqs_min_identity,
                                                    k_best_hits=top_k,
                                                    threads=threads)

    target_ids = np.unique(filtered_mmseqs_results["target"]).tolist()
    target_seqs = retrieve_fasta_entries_as_dict(database / MERGED_SEQUENCES,
                                                 target_ids)

    alignments = align_query(query_seqs, target_seqs, alignment_gap_open,
                             alignment_gap_continuation, identity_threshold,
                             threads)

    aligned_queries = [aln.query_name for aln in alignments]
    unaligned_queries = {
        k: v
        for k, v in query_seqs.items() if k not in aligned_queries
    }
    deepfri_models_config = load_deepfri_config(weights)

    # deepfri_processing_modes = ['mf', 'bp', 'cc', 'ec']
    # mf = molecular_function
    # bp = biological_process
    # cc = cellular_component
    # ec = enzyme_commission

    gcn_prots, cnn_prots = len(aligned_queries), len(unaligned_queries)
    afdb = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/PROPHAGE_REFSEQ_EAN/scratch/databases/afdb_swissprot_v4"
    for mode in deepfri_processing_modes:
        logging.info("Processing mode: %s", mode)
        # GCN for queries with aligned contact map

        if gcn_prots > 0:
            logging.info("Predicting with GCN: %i proteins", gcn_prots)
            output_file_name = output_path / f"results_gcn_{mode}"

            gcn_params = deepfri_models_config["gcn"]["models"][mode]
            gcn = Predictor(gcn_params, threads=threads)

            for aln in alignments:
                query_id = aln.query_name
                logging.info("Predicting %s", query_id)
                string_structure = retrieve_structure(
                    aln.target_name.split(".")[0], afdb)

                sparse_contact_map = calculate_cmap(
                    string_structure,
                    max_seq_len=1000,
                    threshold=angstrom_contact_threshold,
                    mode="sparse")

                aligned_cmap = align_contact_map(
                    aln.gapped_sequence,
                    aln.gapped_target,
                    sparse_contact_map,
                    generated_contacts=generate_contacts)

                # running the actual prediction
                gcn.predict_function(seqres=aln.query_sequence,
                                     cmap=aligned_cmap,
                                     chain=query_id)

            gcn.export_tsv(str(output_file_name.with_suffix('.tsv')))
            del gcn

        # CNN for queries without satisfying alignments
        if cnn_prots > 0:
            logging.info("Predicting with CNN: %i proteins", cnn_prots)
            output_file_name = output_path / f"results_cnn_{mode}"

            cnn_params = deepfri_models_config["cnn"]["models"][mode]
            cnn = Predictor(cnn_params, threads=threads)
            for query_id in unaligned_queries:
                logging.info("Predicting %s", query_id)
                cnn.predict_function(seqres=query_seqs[query_id],
                                     chain=query_id)

            cnn.export_tsv(str(output_file_name.with_suffix('.tsv')))
            del cnn

    logging.info("meta-DeepFRI finished successfully")
