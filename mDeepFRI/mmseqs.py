import logging
import os
import tempfile
from functools import partial
from multiprocessing.pool import ThreadPool
from pathlib import Path

import numpy as np
from pysam import tabix_compress

import mDeepFRI
from mDeepFRI.utils import run_command

MMSEQS_COLUMN_NAMES = [
    "query", "target", "identity", "alignment_length", "mismatches",
    "gap_openings", "query_start", "query_end", "target_start", "target_end",
    "e_value", "bit_score"
]

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def createdb(sequences_file, db_path):
    """
    Converts FASTA file to a DB format needed for MMseqs2.
    This should generate five files,
    e.g. queryDB, queryDB_h and its corresponding index file queryDB.index,
    queryDB_h.index and queryDB.lookup from the FASTA QUERY.fasta input sequences.

    sequence_file (str): path to FASTA file.
    db_path (str): path to output db file.

    Returns:
        None
    """
    run_command(f"mmseqs createdb {sequences_file} {db_path} --dbtype 1")


def createindex(db_path):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(f"mmseqs createindex {db_path} {tmp_path}")


def search(query_db: str, target_db: str, result_db: str, threads: int = 1):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(
            f"mmseqs search --threads {threads} {query_db} {target_db} {result_db} {tmp_path}"
        )


def convertalis(query_db: str, target_db: str, result_db: str,
                output_file: str):
    run_command(
        f"mmseqs convertalis {query_db} {target_db} {result_db} {output_file}")


def extract_fasta_foldcomp(foldcomp_db: str,
                           output_file: str,
                           threads: int = 1):
    """
    Extracts FASTA from database
    """
    foldcomp_bin = Path(mDeepFRI.__path__[0]).parent / "foldcomp_bin"

    # run command
    run_command(
        f"{foldcomp_bin} extract --fasta -t {threads} {foldcomp_db} {output_file}"
    )
    # gzip fasta file
    tabix_compress(output_file, str(output_file) + ".gz", force=True)
    # remove unzipped file
    os.remove(output_file)
    return Path(str(output_file) + ".gz")


def create_target_database(foldcomp_fasta_path: str,
                           mmseqs_db_path: str) -> None:
    """
    Extracts sequences from compressed FoldComp database.

    Args:
        foldcomp_db_path (str): Path to FoldComp database.
        new_db_path (str): Path to new MMSeqs database.

    Returns:
        None
    """
    createdb(foldcomp_fasta_path, mmseqs_db_path)
    logger.info("Indexing new target mmseqs2 database %s", mmseqs_db_path)
    createindex(mmseqs_db_path)


def validate_mmseqs_database(database: str):
    """
    Check if MMSeqs2 database is intact.

    Args:
        database (str): Path to MMSeqs2 database.

    Returns:
        is_valid (bool): True if database is intact.

    """

    # Verify all the files for MMSeqs2 database
    mmseqs2_ext = [
        ".index", ".dbtype", "_h", "_h.index", "_h.dbtype", ".idx",
        ".idx.index", ".idx.dbtype", ".lookup", ".source"
    ]

    target_db = Path(database)
    is_valid = True
    for ext in mmseqs2_ext:
        if not os.path.isfile(f"{target_db}{ext}"):
            logging.debug(f"{target_db}{ext} is missing.")
            is_valid = False
            break

    return is_valid


def run_mmseqs_search(query_file: str,
                      target_db: str,
                      output_path: str,
                      threads: int = 1) -> Path:
    """Creates a database from query sequences and runs mmseqs2 search against database.

    Args:
        query_file (str): Path to query FASTA file.
        target_db (str): Path to target MMSeqs2 database.
        output_path (str): Path to output folder.
        threads (int): Number of threads to use.

    Returns:
        output_file (pathlib.Path): Path to MMSeqs2 search results.
    """
    query_file = Path(query_file)
    target_db = Path(target_db)
    output_path = Path(output_path)

    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / Path(target_db.stem + ".search_results.tsv")
    query_db = str(output_path / 'query.mmseqsDB')
    createdb(query_file, query_db)

    with tempfile.TemporaryDirectory() as tmp_path:

        result_db = str(Path(tmp_path) / 'search_resultDB')
        search(query_db, target_db, result_db, threads=threads)

        # Convert results to tabular format
        convertalis(query_db, target_db, result_db, output_file)

    return output_file


def filter_mmseqs_results(results_file: str,
                          min_bit_score: float = None,
                          max_evalue: float = None,
                          min_identity: float = None,
                          k_best_hits: int = 5,
                          threads: int = 1) -> np.recarray:
    """
    Filters MMSeqs results retrieving only k best hits based on identity
    above specified thresholds. Allows number of paiwise alignments
    in the next step of pipeline.

    Args:
        results_file (pathlib.Path): Path to MMSeqs2 search results.
        min_bit_score (float): Minimum bit score.
        max_evalue (float): Maximum e-value.
        min_identity (float): Minimum identity.
        k_best_hits (int): Number of best hits to keep.
        threads (int): Number of threads to use.

    Returns:
        output (numpy.recarray): Filtered results.
    """
    def select_top_k(query, db, k=30):
        return db[db["query"] == query][:k]

    output = np.recfromcsv(results_file,
                           delimiter="\t",
                           encoding="utf-8",
                           names=MMSEQS_COLUMN_NAMES)

    # check if output is not empty
    if output.size == 0:
        final_database = None

    else:
        # MMSeqs2 alginment filters
        if min_identity:
            output = output[output['identity'] >= min_identity]
        if min_bit_score:
            output = output[output['bit_score'] >= min_bit_score]
        if max_evalue:
            output = output[output['e_value'] <= max_evalue]

        # Get k best hits
        output.sort(order=["query", "identity", "e_value"], kind="quciksort")
        top_k_db = partial(select_top_k, db=output[::-1], k=k_best_hits)
        with ThreadPool(threads) as pool:
            top_k_chunks = pool.map(top_k_db, np.unique(output["query"]))

        final_database = np.concatenate(top_k_chunks)

        logger.info("%i pairs after filtering with k=%i best hits.",
                    final_database.shape[0], k_best_hits)

    return final_database
