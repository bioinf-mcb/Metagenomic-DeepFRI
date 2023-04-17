import tempfile
import os
import pytest

from functools import partial
from Bio import pairwise2

default_pair_align = partial(
    pairwise2.align.globalms, match=2, mismatch=-1, open=-0.5, extend=-0.1, one_alignment_only=True)

from meta_deepFRI.utils import (bio_utils, search_alignments, utils)
# hash_sequence_id, encode_faa_ids, load_fasta_file, write_fasta_file


def test_protein_letters():

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


@pytest.mark.parametrize("alignment, expected_output", [
    (default_pair_align('ATCG', 'ATCG')[0], 1.0),
    (default_pair_align('ATCG', 'AGCG')[0], 0.6),
    (default_pair_align('ATCG', 'CGTA')[0], 0.3),
    (default_pair_align('ATCG', 'ATCGG')[0], 0.8),
])
def test_alignment_sequences_identity(alignment, expected_output):
    # Call the function with the test input
    result = round(search_alignments.alignment_sequences_identity(alignment), 1)

    # Check whether the result matches the expected output
    assert result == expected_output


def test_load_deepfri_config(tmp_path):
    """Test load_deepfri_config() function."""
    tmp_json_file = tmp_path / "model_config.json"
    tmp_json_file.write_text('{"model1": "./trained_models/model1.h5", "model2": "./trained_models/model2.h5"}')

    config = utils.load_deepfri_config(str(tmp_json_file))

    assert config["model1"] == str(tmp_path / "model1.h5")
    assert config["model2"] == str(tmp_path / "model2.h5")


def test_load_deepfri_config_file_not_found():
    """Test FileNotFoundError is raised if config file is not found."""
    with pytest.raises(FileNotFoundError):
        utils.load_deepfri_config("non_existing_file.json")
