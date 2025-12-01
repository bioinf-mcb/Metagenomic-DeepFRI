import csv
import logging
import pathlib
import pickle
import sys
from functools import partial
from multiprocessing import Pool
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np
from tqdm import tqdm

from mDeepFRI import DEEPFRI_MODES
from mDeepFRI.alignment import align_mmseqs_results
from mDeepFRI.bio_utils import build_align_contact_map
from mDeepFRI.database import Database, build_database
from mDeepFRI.mmseqs import MMseqsResult, QueryFile
from mDeepFRI.pdb import create_pdb_mmseqs, extract_calpha_coords
from mDeepFRI.predict import Predictor
from mDeepFRI.utils import (get_json_values, load_deepfri_config,
                            remove_intermediate_files)

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
logger.propagate = False
formatter = logging.Formatter(
    '[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

ALIGNMENT_HEADER = [
    "query_id", "aligned", "target_id", "db_name", "query_identity",
    "query_coverage", "target_coverage"
]
FINAL_OUTPUT_HEADER = [
    "protein", "network_type", "prediction_mode", "go_term", "score",
    "go_name", "aligned", "target_id", "db_name", "query_identity",
    "query_coverage", "target_coverage"
]

NAN_ALIGNMENT_INFO = [np.nan] * 6


# load sequences in query filtered by length
def load_query_file(query_file: str, min_length: int = None, max_length=None):
    query_file = QueryFile(filepath=query_file)
    query_file.load_sequences()
    # filter out sequences
    if min_length or max_length:
        query_file.filter_sequences(
            lambda x: min_length <= len(x) <= max_length)

    return query_file


def hierarchical_database_search(query_file: QueryFile,
                                 output_path: str,
                                 databases: Iterable[str] = [],
                                 mmseqs_sensitivity: float = 5.7,
                                 min_bits: float = 0,
                                 max_eval: float = 1e-5,
                                 min_ident: float = 0.5,
                                 min_coverage: float = 0.9,
                                 top_k: int = 5,
                                 skip_pdb: bool = False,
                                 overwrite: bool = False,
                                 tmpdir: str = None,
                                 threads: int = 1):

    output_path = pathlib.Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # logging variable
    sequence_num_start = len(query_file.sequences)

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
    pdb_hits = set()

    for db in dbs:
        results = query_file.search(db.mmseqs_db,
                                    mmseqs_sensitivity=mmseqs_sensitivity,
                                    eval=max_eval,
                                    threads=threads,
                                    tmpdir=tmpdir)

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
            all_hits = np.unique(best_matches["query"])
        except IndexError:
            all_hits = np.array([])
        # cover skip_pdb case
        unique_hits = all_hits

        if "pdb100" in db.name:
            pdb_hits.update(all_hits)
        elif not skip_pdb:
            unique_hits = [hit for hit in all_hits if hit not in pdb_hits]

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
            query_file.remove_sequences(all_hits)

    return dbs


def align_pairwise():
    pass


def _initialize_processing_modes(modes: List[str],
                                 config: Dict[str, Any]) -> List[str]:
    """
    Filters processing modes based on the model config version.
    """
    filtered_modes = list(modes)
    # version 1.1 drops support for ec
    if config.get("version") == "1.1":
        if "ec" in filtered_modes:
            filtered_modes.remove("ec")
            logger.info(
                "EC number prediction is not supported in version 1.1.")

    if len(filtered_modes) == 0:
        raise ValueError("No processing modes selected.")
    return filtered_modes


def _run_prediction_loop(predictor, data_iterable: iter, data_len: int,
                         net_type: str, tsv_writer: csv.writer,
                         description: str):
    """
    A helper function to run a prediction loop for either GCN or CNN.
    """
    BAR_FORMAT = "{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}], {rate_fmt}{postfix}"
    # Assuming BAR_FORMAT and sys.stdout are available in scope
    # or passed as arguments if this is in a different module.
    for item in tqdm(data_iterable,
                     total=data_len,
                     desc=description,
                     bar_format=BAR_FORMAT,
                     file=sys.stdout):
        if net_type == "gcn":
            # item is (aln, aligned_cmap)
            aln, aligned_cmap = item
            query_id = aln.query_name
            pred_vector = predictor.forward_pass(seqres=aln.query_sequence,
                                                 cmap=aligned_cmap)
        else:  # net_type == "cnn"
            # item is (query_id, sequence)
            query_id, sequence = item
            pred_vector = predictor.forward_pass(seqres=sequence)

        out_row = [query_id, net_type] + pred_vector.tolist()
        tsv_writer.writerow(out_row)


def predict_protein_function(
        query_file: QueryFile,
        databases: Tuple[Database],
        weights: str,
        output_path: str,
        deepfri_processing_modes: List[str] = ["ec", "bp", "mf", "cc"],
        angstrom_contact_threshold: float = 6,
        generate_contacts: int = 2,
        alignment_gap_open: float = 10,
        alignment_gap_continuation: float = 1,
        alignment_min_identity: float = 0.5,
        alignment_min_coverage: float = 0.9,
        remove_intermediate=False,
        threads: int = 1,
        save_structures: bool = False,
        save_cmaps: bool = False):

    # load DeepFRI model
    deepfri_models_config = load_deepfri_config(weights)
    deepfri_processing_modes = _initialize_processing_modes(
        deepfri_processing_modes, deepfri_models_config)

    weights = pathlib.Path(weights)
    output_path = pathlib.Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    aligned_cmaps = []
    for db in databases:
        # SEQUENCE ALIGNMENT
        # calculate already aligned sequences
        alignments = align_mmseqs_results(
            best_matches_filepath=db.mmseqs_result,
            sequence_db=db.sequence_db,
            alignment_gap_open=alignment_gap_open,
            alignment_gap_extend=alignment_gap_continuation,
            threads=threads)

        # filter alignments by identity and coverage
        alignments = [
            aln for aln in alignments
            if aln.query_identity > alignment_min_identity
            and aln.query_coverage > alignment_min_coverage
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

            # remove broken structures
            if db.name == "highquality_clust30":
                data_path = pathlib.Path(__file__).parent / "assets"
                # convert to abspath
                data_path = data_path.resolve()
                with open(data_path / "highquality_clust30_error_ids.pkl",
                          "rb") as f:
                    error_ids = pickle.load(f)
                # filter out broken structures
                new_alignments = {
                    query_name: aln
                    for query_name, aln in new_alignments.items()
                    if aln.target_name not in error_ids
                }

            query_ids = [aln.query_name for aln in new_alignments.values()]
            target_ids = [
                aln.target_name.rsplit(".", 1)[0]
                for aln in new_alignments.values()
            ]

            # extract structural information
            # in form of C-alpha coordinates
            if save_structures:
                save_dir = output_path / "structures" / db.name
                save_dir.mkdir(parents=True, exist_ok=True)
            else:
                save_dir = None

            coords = extract_calpha_coords(db,
                                           target_ids,
                                           query_ids,
                                           save_directory=save_dir,
                                           threads=threads)

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

    if save_cmaps:
        cmap_dir = output_path / "contact_maps"
        cmap_dir.mkdir(parents=True, exist_ok=True)
        for i, (aln, cmap) in enumerate(aligned_cmaps):
            cmap_file = cmap_dir / f"{aln.query_name}.npy"
            np.save(cmap_file, cmap)

    aligned_queries = [aln[0].query_name for aln in aligned_cmaps]
    unaligned_queries = {
        query_id: seq
        for query_id, seq in query_file.sequences.items()
        if query_id not in aligned_queries
    }

    # WRITE ALIGNMENT RESULTS
    alignment_results_file = output_path / "alignment_summary.tsv"

    with open(alignment_results_file, "w", encoding="utf-8") as aln_output:
        tsv_writer = csv.writer(aln_output, delimiter="\t")
        tsv_writer.writerow(ALIGNMENT_HEADER)
        for aln, _ in aligned_cmaps:
            tsv_writer.writerow([
                aln.query_name, True, aln.target_name, aln.db_name,
                aln.query_identity, aln.query_coverage, aln.target_coverage
            ])
        for query_id in unaligned_queries:
            tsv_writer.writerow(
                [query_id, False, np.nan, np.nan, np.nan, np.nan, np.nan])

    ### FUNCTION PREDICTION ###
    # sort cmaps by length of query sequence
    aligned_cmaps = sorted(aligned_cmaps,
                           key=lambda x: len(x[0].query_sequence))
    # sort unaligned queries by length
    unaligned_queries = dict(
        sorted(unaligned_queries.items(), key=lambda x: len(x[1])))

    # output_file_name = output_path / "results.tsv"
    # output_buffer = open(output_file_name, "w", encoding="utf-8")
    # csv_writer = csv.writer(output_buffer, delimiter="\t")
    # csv_writer.writerow(OUTPUT_HEADER)

    matrices = {}
    json_configs = {}
    for i, mode in enumerate(deepfri_processing_modes):
        # load model go terms
        model_path = deepfri_models_config["gcn"][mode]
        config_path = model_path.rsplit(".", 1)[0] + "_model_params.json"
        json_configs[mode] = config_path
        GOTERMS = get_json_values(config_path, "goterms")

        # create output file for each mode
        output_matrix = output_path / f"prediction_matrix_{mode}.tsv"
        matrices[mode] = output_matrix
        with open(output_matrix, "w", encoding="utf-8") as output_buffer:
            tsv_writer = csv.writer(output_buffer, delimiter="\t")
            tsv_writer.writerow(["protein", "network_type"] + GOTERMS)

            logger.info("Processing mode: %s; %i/%i", DEEPFRI_MODES[mode],
                        i + 1, len(deepfri_processing_modes))

            # GCN prediction
            gcn_prots = len(aligned_cmaps)
            if gcn_prots > 0:
                net_type = "gcn"

                # GCN for queries with aligned contact map
                gcn_path = deepfri_models_config[net_type][mode]
                gcn = Predictor(gcn_path, threads=threads)
                _run_prediction_loop(
                    predictor=gcn,
                    data_iterable=aligned_cmaps,
                    data_len=len(aligned_cmaps),
                    net_type=net_type,
                    tsv_writer=tsv_writer,
                    description=f"Predicting with GCN ({DEEPFRI_MODES[mode]})")
                del gcn  # Explicitly free memory

            # CNN for queries without satisfying alignments
            cnn_prots = len(unaligned_queries)
            if cnn_prots > 0:
                net_type = "cnn"
                cnn_path = deepfri_models_config[net_type][mode]
                cnn = Predictor(cnn_path, threads=threads)
                _run_prediction_loop(
                    predictor=cnn,
                    data_iterable=unaligned_queries.items(),
                    data_len=len(unaligned_queries),
                    net_type=net_type,
                    tsv_writer=tsv_writer,
                    description=f"Predicting with CNN ({DEEPFRI_MODES[mode]})")
                del cnn  # Explicitly free memory

    ### FORMAT AND CREATE FINAL OUTPUT FILES ###
    # combine mode-specific matrices into a single file
    # open and load alignment data
    with open(alignment_results_file, "r", encoding="utf-8") as aln_input:
        tsv_reader = csv.reader(aln_input, delimiter="\t")
        next(tsv_reader)  # skip header
        alignment_data = {row[0]: row[1:] for row in tsv_reader}

    final_output = output_path / "results.tsv"
    with open(final_output, "w", encoding="utf-8") as fout:
        fout.write("\t".join(FINAL_OUTPUT_HEADER) + "\n")
        for mode, matrix_file in matrices.items():
            json_path = json_configs[mode]
            GONAMES = get_json_values(json_path, "gonames")
            with open(matrix_file, "r", encoding="utf-8") as matrix_input:
                # convert wide to narrow format for GO terms with score > 0.1
                tsv_reader = csv.reader(matrix_input, delimiter="\t")
                # get term names from header
                header = next(tsv_reader)
                terms = header[2:]  # skip first two columns (Protein and Type)
                term_to_name = {
                    term: name
                    for term, name in zip(terms, GONAMES)
                }
                # get ids with scores > 0.1
                for row in tsv_reader:
                    query_id = row[0]
                    net_type = row[1]
                    scores = row[2:]
                    term_score = {
                        terms[i]: float(scores[i])
                        for i in range(len(terms)) if float(scores[i]) >= 0.1
                    }
                    sorted_term_score = dict(
                        sorted(term_score.items(),
                               key=lambda item: item[1],
                               reverse=True))
                    # print results and add go names
                    for term, score in sorted_term_score.items():
                        go_name = term_to_name.get(term, "Unknown")
                        aln_info = alignment_data.get(query_id, [np.nan] * 6)
                        fout.write(
                            f"{query_id}\t{net_type}\t{DEEPFRI_MODES[mode]}\t{term}\t{score:.4f}\t{go_name}\t"
                            f"\t{aln_info[1]}\t{aln_info[2]}\t{aln_info[3]}\t{aln_info[4]}\n"
                        )

    if remove_intermediate:
        for db in databases:
            remove_intermediate_files([db.sequence_db, db.mmseqs_db])

    logger.info("meta-DeepFRI finished successfully.")
