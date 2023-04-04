import time
import pytest
import tempfile
import os.path

from meta_deepFRI.utils import (bio_utils, elapsed_time_logger)


def test_protein_letters():
    """
    Test protein letter encoding dictionaries

    Returns:
        None
    """

    expected = {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E',
        'PHE': 'F',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LYS': 'K',
        'LEU': 'L',
        'MET': 'M',
        'ASN': 'N',
        'PRO': 'P',
        'GLN': 'Q',
        'ARG': 'R',
        'SER': 'S',
        'THR': 'T',
        'VAL': 'V',
        'TRP': 'W',
        'TYR': 'Y',
        'ASX': 'B',
        'XAA': 'X',
        'GLX': 'Z',
        'XLE': 'J',
        'SEC': 'U',
        'PYL': 'O',
        'UNK': 'X'
    }
    assert bio_utils.PROTEIN_LETTERS == expected


def test_elapsed_time_logger():
    """
    Test ElapsedTimeLogger

    Returns:
        None
    """

    # Create a temporary file for logging
    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        tmp_file_path = tmp_file.name

    # Create an ElapsedTimeLogger instance
    logger = elapsed_time_logger.ElapsedTimeLogger(path=tmp_file_path)

    # Log some times
    logger.log("step1")
    time.sleep(1)
    logger.log("step2")

    # Log total time
    logger.log_total_time()

    # Read the log file and check its contents
    with open(tmp_file_path, "r", encoding="utf-8") as f:
        log_contents = f.read()

    # extract numbers from log_contents
    log_contents = log_contents.replace("step1,", "")
    log_contents = log_contents.replace("step2,", "")
    log_contents = log_contents.replace("total_time,", "")
    log_contents = log_contents.replace("\n", ",")
    log_contents = log_contents[:-1]    # remove last comma
    # turn into floats
    log_contents = [float(x) for x in log_contents.split(",")]

    # Check that the first two numbers are close to 1 second
    assert pytest.approx(0, abs=0.01) == log_contents[0]
    assert pytest.approx(1, abs=0.01) == log_contents[1]
    # Check that the last number is close to 1 second1
    assert pytest.approx(1, abs=0.01) == log_contents[2]

    # Delete the temporary file
    os.remove(tmp_file_path)
