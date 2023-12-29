import csv
import logging
import pathlib
from functools import partial
from multiprocessing import Pool
from typing import List, Tuple

import numpy as np

from mDeepFRI.alignment import align_query
from mDeepFRI.bio_utils import (load_query_sequences,
                                retrieve_align_contact_map,
                                retrieve_fasta_entries_as_dict)
from mDeepFRI.database import build_database
from mDeepFRI.mmseqs import (filter_mmseqs_results, run_mmseqs_search,
                             validate_mmseqs_database)
from mDeepFRI.predict import Predictor
from mDeepFRI.utils import load_deepfri_config, remove_intermediate_files

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


# TODO: input an MMSeqs DB as query
def predict_protein_function(
        query_file: str,
        databases: Tuple[str],
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
        remove_intermediate=False,
        overwrite=False,
        threads: int = 1):

    MAX_SEQ_LEN = 1000
    query_file = pathlib.Path(query_file)
    weights = pathlib.Path(weights)
    output_path = pathlib.Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    deepfri_models_config = load_deepfri_config(weights)

    # version 1.1 drops support for ec
    if deepfri_models_config["version"] == "1.1":
        # remove "ec" from processing modes
        deepfri_processing_modes = [
            mode for mode in deepfri_processing_modes if mode != "ec"
        ]
        logging.info("EC number prediction is not supported in version 1.1.")

    assert len(
        deepfri_processing_modes) > 0, "No valid processing modes selected."

    # design solution
    # database is built in the same directory
    # where the structure database is stored

    deepfri_dbs = []
    for database in databases:
        database = pathlib.Path(database)
        db = build_database(
            input_path=database,
            output_path=database.parent,
            overwrite=overwrite,
            threads=threads,
        )
        deepfri_dbs.append(db)

    # SEQUENCE ALIGNMENT
    query_seqs = load_query_sequences(query_file, output_path)
    logging.info("Aligning %i sequences with MMSeqs2.", len(query_seqs))
    aligned_queries = {}

    for db in deepfri_dbs:
        mmseqs_valid = validate_mmseqs_database(db.mmseqs_db)
        if not mmseqs_valid:
            raise ValueError(
                f"MMSeqs2 database not found in {database.parent}.")

        mmseqs_results = run_mmseqs_search(query_file, db.mmseqs_db,
                                           output_path, threads)

        filtered_mmseqs_results = filter_mmseqs_results(mmseqs_results,
                                                        mmseqs_min_bit_score,
                                                        mmseqs_max_eval,
                                                        mmseqs_min_identity,
                                                        k_best_hits=top_k,
                                                        threads=threads)

        # if MMSeqs2 alignment is empty
        if filtered_mmseqs_results is None:
            pass
        else:
            target_ids = np.unique(filtered_mmseqs_results["target"]).tolist()
            target_seqs = retrieve_fasta_entries_as_dict(
                db.sequence_db, target_ids)

            alignments = align_query(query_seqs, target_seqs,
                                     alignment_gap_open,
                                     alignment_gap_continuation,
                                     identity_threshold, threads)

            # set a db name for alignments
            for aln in alignments:
                aln.db_name = db.name

            aligned_queries.update({aln.query_name: aln for aln in alignments})

    unaligned_queries = {
        k: v
        for k, v in query_seqs.items() if k not in aligned_queries.keys()
    }

    # deepfri_processing_modes = ['mf', 'bp', 'cc', 'ec']
    # mf = molecular_function
    # bp = biological_process
    # cc = cellular_component
    # ec = enzyme_commission

    gcn_prots, cnn_prots = len(aligned_queries), len(unaligned_queries)

    # CONTACT MAP ALIGNMENT
    aligned_cmaps = []
    if gcn_prots > 0:
        for db in deepfri_dbs:
            # select only alignments from db
            db_alignments = [
                v for v in aligned_queries.values() if v.db_name == db.name
            ]

            partial_align = partial(retrieve_align_contact_map,
                                    database=db.foldcomp_db,
                                    max_seq_len=MAX_SEQ_LEN,
                                    threshold=angstrom_contact_threshold,
                                    generated_contacts=generate_contacts)
            logging.info("Aligning contact maps for %i proteins", gcn_prots)

            with Pool(threads) as p:
                partial_cmaps = p.map(partial_align, db_alignments)

            aligned_cmaps.extend(partial_cmaps)

    logging.info("Aligned %i contact maps", len(aligned_cmaps))

    output_file_name = output_path / "results.tsv"
    output_buffer = open(output_file_name, "w", encoding="utf-8")
    csv_writer = csv.writer(output_buffer, delimiter="\t")
    csv_writer.writerow([
        'Protein',
        'GO_term/EC_number',
        'Score',
        'Annotation',
        'Neural_net',
        'DeepFRI_mode',
        'DB_hit',
        'DB_name',
        'Identity',
    ])

    # FUNCTION PREDICTION
    for i, mode in enumerate(deepfri_processing_modes):
        # GCN
        if gcn_prots > 0:
            net_type = "gcn"
            logging.info("Processing mode: %s; %i/%i", mode, i + 1,
                         len(deepfri_processing_modes))
            # GCN for queries with aligned contact map
            logging.info("Predicting with GCN: %i proteins", gcn_prots)
            gcn_path = deepfri_models_config[net_type][mode]

            gcn = Predictor(gcn_path, threads=threads)

            for i, (aln, aligned_cmap) in enumerate(aligned_cmaps):
                logging.info("Predicting %s; %i/%i", aln.query_name, i + 1,
                             gcn_prots)
                # running the actual prediction
                prediction_rows = gcn.predict_function(
                    seqres=aln.query_sequence,
                    cmap=aligned_cmap,
                    chain=aln.query_name)
                # writing the results to the output file
                for row in prediction_rows:
                    row.extend([net_type, mode])
                    # corrected name for FoldComp inconsistency
                    row.extend([
                        aln.target_name.rsplit(".", 1)[0], aln.db_name,
                        aln.identity
                    ])
                    csv_writer.writerow(row)

            del gcn

        # CNN for queries without satisfying alignments
        if cnn_prots > 0:
            net_type = "cnn"
            logging.info("Predicting with CNN: %i proteins", cnn_prots)
            cnn_path = deepfri_models_config[net_type][mode]
            cnn = Predictor(cnn_path, threads=threads)
            for i, query_id in enumerate(unaligned_queries):
                logging.info("Predicting %s; %i/%i", query_id, i + 1,
                             cnn_prots)
                prediction_rows = cnn.predict_function(
                    seqres=query_seqs[query_id], chain=query_id)
                for row in prediction_rows:
                    row.extend([net_type, mode])
                    row.extend([np.nan, np.nan, np.nan])
                    csv_writer.writerow(row)

            del cnn

    output_buffer.close()

    # sort predictions by score
    with open(output_file_name, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = sorted(reader, key=lambda row: float(row[2]), reverse=True)

    with open(output_file_name, "w", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)

    if remove_intermediate:
        for db in deepfri_dbs:
            remove_intermediate_files([db.sequence_db, db.mmseqs_db])

    logging.info("meta-DeepFRI finished successfully.")
