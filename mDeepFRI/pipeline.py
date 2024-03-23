import csv
import logging
import pathlib
from functools import partial
from multiprocessing import Pool
from typing import List, Tuple

import numpy as np

from mDeepFRI.alignment import run_alignment
from mDeepFRI.bio_utils import load_fasta_as_dict, retrieve_align_contact_map
from mDeepFRI.database import build_database
from mDeepFRI.pdb import create_pdb_mmseqs
from mDeepFRI.predict import Predictor
from mDeepFRI.utils import load_deepfri_config, remove_intermediate_files

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def predict_protein_function(
        query_file: str,
        databases: Tuple[str],
        weights: str,
        output_path: str,
        deepfri_processing_modes: List[str] = ["ec", "bp", "mf", "cc"],
        angstrom_contact_threshold: float = 6,
        generate_contacts: int = 2,
        mmseqs_min_bitscore: float = None,
        mmseqs_max_eval: float = 10e-5,
        mmseqs_min_identity: float = 0.5,
        top_k: int = 5,
        alignment_gap_open: float = 10,
        alignment_gap_continuation: float = 1,
        identity_threshold: float = 0.5,
        remove_intermediate=False,
        overwrite=False,
        threads: int = 1,
        skip_pdb: bool = False,
        min_length: int = 60,
        max_length: int = 1000):

    MIN_SEQ_LEN = min_length
    MAX_SEQ_LEN = max_length
    logger.info("DeepFRI protein sequence limit: %i-%i aa.", MIN_SEQ_LEN,
                MAX_SEQ_LEN)

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
        logger.info("EC number prediction is not supported in version 1.1.")

    assert len(
        deepfri_processing_modes) > 0, "No valid processing modes selected."

    deepfri_dbs = []
    # PDB100 database
    if not skip_pdb:
        logger.info(
            "Creating PDB100 database. This may take a bit during a first run."
        )
        pdb100 = create_pdb_mmseqs()
        deepfri_dbs.append(pdb100)
        logger.info("PDB100 database created.")

    # design solution
    # database is built in the same directory
    # where the structure database is stored
    for database in databases:
        database = pathlib.Path(database)
        db = build_database(
            input_path=database,
            output_path=database.parent,
            overwrite=overwrite,
            threads=threads,
        )
        deepfri_dbs.append(db)

    query_seqs = load_fasta_as_dict(query_file)
    aligned_cmaps = []

    for db in deepfri_dbs:
        # SEQUENCE ALIGNMENT
        # calculate already aligned sequences
        aligned = len(aligned_cmaps)
        logger.info("Aligning %s sequences against %s.", len(query_seqs),
                    db.name)

        # align new sequences agains db
        alignments = run_alignment(query_file, db.mmseqs_db, db.sequence_db,
                                   output_path, mmseqs_min_bitscore,
                                   mmseqs_max_eval, mmseqs_min_identity, top_k,
                                   alignment_gap_open,
                                   alignment_gap_continuation, threads)

        # if anything aligned
        if not alignments:
            logger.info("No alignments found for %s.", db.name)
            continue
        # filter alignments by identity
        alignments = [
            aln for aln in alignments if aln.identity > identity_threshold
        ]

        if not alignments:
            logger.info("All alignments below identity threshold for %s.",
                        db.name)
            continue

        # set a db name for alignments
        for aln in alignments:
            aln.db_name = db.name

        aligned_queries = [aln[0].query_name for aln in aligned_cmaps]
        new_alignments = {
            aln.query_name: aln
            for aln in alignments if aln.query_name not in aligned_queries
        }

        # CONTACT MAP ALIGNMENT
        # initially designed as a separate step
        # some protein structures in PDB are not formatted correctly
        # so contact map alignment fails for them
        # for this cases we replace closest experimental structure with
        # closest predicted structure if available
        # if no alignments were found - report
        logger.info("Aligning contact maps for %i proteins against %s.",
                    len(new_alignments), db.name)

        partial_align = partial(retrieve_align_contact_map,
                                database=db.foldcomp_db,
                                threshold=angstrom_contact_threshold,
                                generated_contacts=generate_contacts)

        with Pool(threads) as p:
            partial_cmaps = p.map(partial_align, new_alignments.values())

        # filter errored contact maps
        # returned as Tuple[AlignmentResult, None] from `retrieve_align_contact_map`
        partial_cmaps = [cmap for cmap in partial_cmaps if cmap[1] is not None]
        aligned_cmaps.extend(partial_cmaps)

    # report alignments for each database
    logger.info("Aligned %i sequences.", len(aligned_queries) - aligned)

    aligned_queries = [aln.query_name for aln in alignments]
    unaligned_queries = {
        k: v
        for k, v in query_seqs.items() if k not in aligned_queries
    }

    # deepfri_processing_modes = ['mf', 'bp', 'cc', 'ec']
    # mf = molecular_function
    # bp = biological_process
    # cc = cellular_component
    # ec = enzyme_commission

    # sort cmaps by length of query sequence
    aligned_cmaps = sorted(aligned_cmaps,
                           key=lambda x: len(x[0].query_sequence))

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
        gcn_prots = len(aligned_cmaps)
        if gcn_prots > 0:
            net_type = "gcn"
            logger.info("Processing mode: %s; %i/%i", mode, i + 1,
                        len(deepfri_processing_modes))
            # GCN for queries with aligned contact map
            logger.info("Predicting with GCN: %i proteins", gcn_prots)
            gcn_path = deepfri_models_config[net_type][mode]

            gcn = Predictor(gcn_path, threads=threads)

            for i, (aln, aligned_cmap) in enumerate(aligned_cmaps):

                ### PROTEIN LENGTH CHECKS
                if len(aln.query_sequence) < MIN_SEQ_LEN:
                    logger.info("Skipping %s; sequence too short %i aa",
                                aln.query_name, len(aln.query_sequence))
                    continue

                elif len(aln.query_sequence) > MAX_SEQ_LEN:
                    logger.info("Skipping %s; sequence too long - %i aa",
                                aln.query_name, len(aln.query_sequence))
                    continue

                logger.info("Predicting %s; %i/%i", aln.query_name, i + 1,
                            gcn_prots)
                # running the actual prediction
                prediction_rows = gcn.predict_function(
                    seqres=aln.query_sequence,
                    cmap=aligned_cmap,
                    chain=str(aln.query_name))
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
        cnn_prots = len(unaligned_queries)
        if cnn_prots > 0:
            net_type = "cnn"
            logger.info("Predicting with CNN: %i proteins", cnn_prots)
            cnn_path = deepfri_models_config[net_type][mode]
            cnn = Predictor(cnn_path, threads=threads)
            for i, query_id in enumerate(unaligned_queries):
                seq = query_seqs[query_id]
                if len(seq) > MAX_SEQ_LEN:
                    logger.info("Skipping %s; sequence too long %i", query_id,
                                len(seq))
                    continue

                elif len(seq) < MIN_SEQ_LEN:
                    logger.info("Skipping %s; sequence too short %i", query_id,
                                len(seq))
                    continue

                logger.info("Predicting %s; %i/%i", query_id, i + 1, cnn_prots)
                prediction_rows = cnn.predict_function(
                    seqres=query_seqs[query_id], chain=str(query_id))
                for row in prediction_rows:
                    row.extend([net_type, mode])
                    row.extend([np.nan, np.nan, np.nan])
                    csv_writer.writerow(row)

            del cnn

    output_buffer.close()

    # sort predictions by protein name and score
    with open(output_file_name, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        # row[0] - protein name
        # row[2] - DeepFRI score
        rows = sorted(reader, key=lambda row: (str(row[0]), -float(row[2])))

    with open(output_file_name, "w", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)

    if remove_intermediate:
        for db in deepfri_dbs:
            remove_intermediate_files([db.sequence_db, db.mmseqs_db])

    logger.info("meta-DeepFRI finished successfully.")
