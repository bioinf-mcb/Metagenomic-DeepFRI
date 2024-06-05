import csv
import logging
import pathlib
import sys
from functools import partial
from multiprocessing import Pool
from typing import Iterable, List, Tuple

import numpy as np
from tqdm import tqdm

from mDeepFRI import BAR_FORMAT, DEEPFRI_MODES, OUTPUT_HEADER
from mDeepFRI.alignment import align_mmseqs_results
from mDeepFRI.bio_utils import build_align_contact_map
from mDeepFRI.database import build_database
from mDeepFRI.mmseqs import MMseqsResult, QueryFile
from mDeepFRI.pdb import create_pdb_mmseqs, extract_calpha_coords
from mDeepFRI.predict import Predictor
from mDeepFRI.utils import load_deepfri_config, remove_intermediate_files

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    '[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def hierarchical_database_search(query_file: str,
                                 output_path: str,
                                 databases: Iterable[str] = [],
                                 sensitivity: float = 5.7,
                                 min_seq_len: int = None,
                                 max_seq_len: int = None,
                                 min_bits: float = 0,
                                 max_eval: float = 1e-5,
                                 min_ident: float = 0.5,
                                 min_coverage: float = 0.9,
                                 top_k: int = 5,
                                 skip_pdb: bool = False,
                                 overwrite: bool = False,
                                 threads: int = 1):

    output_path = pathlib.Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    # load initial sequences
    query_file = QueryFile(filepath=query_file)
    query_file.load_sequences()
    # logging variable
    sequence_num_start = len(query_file.sequences)
    # filter out sequences
    if min_seq_len or max_seq_len:
        query_file.filter_sequences(
            lambda x: min_seq_len <= len(x) <= max_seq_len)

    for idx, seq in query_file.filtered_out.items():
        logger.info(f"Skipping {idx}; sequence length {len(seq)} aa.")

    dbs = []
    # PDB100 database
    if not skip_pdb:
        logger.info("Creating PDB100 database.")
        pdb100 = create_pdb_mmseqs(threads=threads)
        dbs.append(pdb100)
        logger.info("PDB100 database created.")

    for database in databases:
        database = pathlib.Path(database)
        db = build_database(
            input_path=database,
            output_path=database.parent,
            overwrite=overwrite,
            threads=threads,
        )
        dbs.append(db)

    aligned_total = 0

    seqs_to_search = query_file.sequences.copy()
    for db in dbs:

        results = query_file.search(db.mmseqs_db,
                                    sensitivity=sensitivity,
                                    eval=max_eval,
                                    threads=threads)

        filtered = results.apply_filters(min_cov=min_coverage,
                                         min_bits=min_bits,
                                         min_ident=min_ident)

        try:
            best_matches = filtered.find_best_matches(top_k, threads=threads)
        except ValueError:
            best_matches = MMseqsResult([], results.query_fasta,
                                        db.sequence_db)

        mmseqs_results_path = output_path / f"{db.name}_results.tsv"
        # save intermediate results
        best_matches.save(mmseqs_results_path)
        # store the location of the result for the next step
        db.mmseqs_result = mmseqs_results_path

        # catch error if no matches to database
        # a case from phage proteins
        try:
            unique_hits = np.unique(best_matches["query"])
        except IndexError:
            unique_hits = np.array([])

        if "pdb100" in db.name:
            pdb_hits = unique_hits
        else:
            unique_hits = [hit for hit in unique_hits if hit not in pdb_hits]
        aligned_db = len(unique_hits)
        aligned_total += aligned_db
        aligned_perc = round(aligned_db / sequence_num_start * 100, 2)
        total_perc = round(aligned_total / sequence_num_start * 100, 2)
        logger.info(f"Aligned {aligned_db}/{sequence_num_start} "
                    f"({aligned_perc:.2f}%) proteins against {db.name}.")
        logger.info(
            f"Aligned {aligned_total}/{sequence_num_start} ({total_perc:.2f}%) proteins in total."
        )

        # this mechanism decreases the amount of sequences
        # on each iteration. Drastically improves execution times
        # for large datasets.
        # PDB100 hits are aligned second time to experimental
        # structures in order to save failed contact map alignemnts.
        if 'pdb100' not in db.name:
            query_file.remove_sequences(unique_hits)

    # reassign initial sequences for further processing
    query_file.sequences = seqs_to_search

    return query_file, dbs


def predict_protein_function(
        query_file: str,
        databases: Tuple[str],
        weights: str,
        output_path: str,
        deepfri_processing_modes: List[str] = ["ec", "bp", "mf", "cc"],
        angstrom_contact_threshold: float = 6,
        generate_contacts: int = 2,
        mmseqs_sensitivity: float = 5.7,
        mmseqs_min_bitscore: float = 0,
        mmseqs_max_eval: float = 1e-5,
        mmseqs_min_identity: float = 0.5,
        mmseqs_min_coverage: float = 0.9,
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

    # load DeepFRI model
    deepfri_models_config = load_deepfri_config(weights)
    # version 1.1 drops support for ec
    if deepfri_models_config["version"] == "1.1":
        # remove "ec" from processing modes
        deepfri_processing_modes = [
            mode for mode in deepfri_processing_modes if mode != "ec"
        ]
        logger.info("EC number prediction is not supported in version 1.1.")

    if len(deepfri_processing_modes) == 0:
        raise ValueError("No processing modes selected.")

    weights = pathlib.Path(weights)
    output_path = pathlib.Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    query_file, deepfri_dbs = hierarchical_database_search(
        query_file=query_file,
        output_path=output_path / "database_search",
        databases=databases,
        sensitivity=mmseqs_sensitivity,
        min_seq_len=min_length,
        max_seq_len=max_length,
        min_bits=mmseqs_min_bitscore,
        max_eval=mmseqs_max_eval,
        min_ident=mmseqs_min_identity,
        min_coverage=mmseqs_min_coverage,
        top_k=top_k,
        skip_pdb=skip_pdb,
        overwrite=overwrite,
        threads=threads)

    aligned_cmaps = []
    for db in deepfri_dbs:
        # SEQUENCE ALIGNMENT
        # calculate already aligned sequences
        alignments = align_mmseqs_results(
            best_matches_filepath=db.mmseqs_result,
            sequence_db=db.sequence_db,
            alignment_gap_open=alignment_gap_open,
            alignment_gap_extend=alignment_gap_continuation,
            threads=threads)

        # filter alignments by identity
        alignments = [
            aln for aln in alignments
            if aln.query_identity > identity_threshold
        ]

        try:
            # set a db name for alignments
            for aln in alignments:
                aln.db_name = db.name

            aligned_queries = [aln[0].query_name for aln in aligned_cmaps]
            new_alignments = {
                aln.query_name: aln
                for aln in alignments if aln.query_name not in aligned_queries
                and aln.query_name in query_file.sequences
            }

            # CONTACT MAP ALIGNMENT
            # initially designed as a separate step
            # some protein structures in PDB are not formatted correctly
            # so contact map alignment fails for them
            # for this cases we replace closest experimental structure with
            # closest predicted structure if available
            # if no alignments were found - report

            query_ids = [aln.query_name for aln in new_alignments.values()]
            target_ids = [
                aln.target_name.rsplit(".", 1)[0]
                for aln in new_alignments.values()
            ]

            # extract structural information
            # in form of C-alpha coordinates
            coords = extract_calpha_coords(db, target_ids, query_ids, threads)

            for aln, coord in zip(new_alignments.values(), coords):
                aln.coords = coord

        # troubleshoot cases where alignments are empty
        except IndexError:
            logger.info("No alignments found for %s.", db.name)
            new_alignments = {}
            continue

        # if new alignments are empty - result is empty as well
        partial_map_align = partial(build_align_contact_map,
                                    threshold=angstrom_contact_threshold,
                                    generated_contacts=generate_contacts)

        with Pool(threads) as p:
            cmaps = list(p.map(partial_map_align, new_alignments.values()))

        # filter errored contact maps
        # returned as Tuple[AlignmentResult, None] from `retrieve_align_contact_map`
        partial_cmaps = [cmap for cmap in cmaps if cmap[1] is not None]
        aligned_cmaps.extend(partial_cmaps)
        aligned_database = round(
            len(partial_cmaps) / len(query_file.sequences) * 100, 2)
        aligned_total = round(
            len(aligned_cmaps) / len(query_file.sequences) * 100, 2)
        logger.info(
            f"Aligned {len(partial_cmaps)}/{len(query_file.sequences)} ({aligned_database}%) "
            f"proteins against {db.name} [without length ivalid].")
        logger.info(
            f"Aligned {len(aligned_cmaps)}/{len(query_file.sequences)} ({aligned_total}%) "
            "proteins in total [without length invalid].")

    # FUNCTION PREDICTION
    aligned_queries = [aln[0].query_name for aln in aligned_cmaps]
    unaligned_queries = {
        query_id: seq
        for query_id, seq in query_file.sequences.items()
        if query_id not in aligned_queries
    }

    # sort cmaps by length of query sequence
    aligned_cmaps = sorted(aligned_cmaps,
                           key=lambda x: len(x[0].query_sequence))
    # sort unaligned queries by length
    unaligned_queries = dict(
        sorted(unaligned_queries.items(), key=lambda x: len(x[1])))

    output_file_name = output_path / "results.tsv"
    output_buffer = open(output_file_name, "w", encoding="utf-8")
    csv_writer = csv.writer(output_buffer, delimiter="\t")
    csv_writer.writerow(OUTPUT_HEADER)

    for i, mode in enumerate(deepfri_processing_modes):
        logger.info("Processing mode: %s; %i/%i", DEEPFRI_MODES[mode], i + 1,
                    len(deepfri_processing_modes))
        # GCN
        gcn_prots = len(aligned_cmaps)
        if gcn_prots > 0:
            net_type = "gcn"

            # GCN for queries with aligned contact map
            gcn_path = deepfri_models_config[net_type][mode]

            gcn = Predictor(gcn_path, threads=threads)

            for i, (aln, aligned_cmap) in tqdm(
                    enumerate(aligned_cmaps),
                    total=gcn_prots,
                    miniters=len(aligned_cmaps) // 10,
                    desc=f"Predicting with GCN ({DEEPFRI_MODES[mode]})",
                    bar_format=BAR_FORMAT):
                # writing the results to the output file

                prediction_rows = gcn.predict_function(
                    seqres=aln.query_sequence,
                    cmap=aligned_cmap,
                    chain=str(aln.query_name))

                for row in prediction_rows:
                    deepfri_info = [net_type, mode]
                    row.extend(deepfri_info)

                    # additional alignment info
                    # corrected name for FoldComp inconsistency

                    row.extend([
                        aln.target_name.rsplit(".", 1)[0], aln.db_name,
                        aln.query_identity, aln.query_coverage
                    ])
                    csv_writer.writerow(row)

            del gcn

        # CNN for queries without satisfying alignments
        cnn_prots = len(unaligned_queries)
        if cnn_prots > 0:
            net_type = "cnn"
            cnn_path = deepfri_models_config[net_type][mode]
            cnn = Predictor(cnn_path, threads=threads)
            for i, query_id in tqdm(
                    enumerate(unaligned_queries),
                    total=cnn_prots,
                    miniters=len(unaligned_queries) // 10,
                    desc=f"Predicting with CNN ({DEEPFRI_MODES[mode]})",
                    bar_format=BAR_FORMAT):

                prediction_rows = cnn.predict_function(
                    seqres=unaligned_queries[query_id], chain=str(query_id))
                for row in prediction_rows:
                    row.extend([net_type, mode])
                    row.extend([np.nan, np.nan, np.nan])
                    csv_writer.writerow(row)

            del cnn

    output_buffer.close()

    if remove_intermediate:
        for db in deepfri_dbs:
            remove_intermediate_files([db.sequence_db, db.mmseqs_db])

    logger.info("meta-DeepFRI finished successfully.")
