# Create logger
import logging
from pathlib import Path

from mDeepFRI import TARGET_MMSEQS_DB_NAME
from mDeepFRI.mmseqs import create_target_database, extract_fasta_foldcomp

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def build_database(
    input_path: str,
    output_path: str,
    threads: int,
) -> None:
    """
    Extracts FASTA file from FoldComp database. Creates MMSeqs2 database and index.

    Args:
        input_path (str): path to a FoldComp database with compressed structures.
        output_path (str): path to folder where the database for DeepFRI will be created.
        threads (int): number of threads to use.

    Returns:
        None
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    output_sequences = output_path / (input_path.stem + ".fasta")
    extract_fasta_foldcomp(input_path, output_sequences, threads)
    create_target_database(output_sequences,
                           output_path / TARGET_MMSEQS_DB_NAME)
