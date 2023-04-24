import pytest

from meta_deepFRI.utils import bio_utils, utils


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


def test_run_command():
    command = "echo 'Hello World!'"
    result = utils.run_command(command).strip()
    assert result == "Hello World!"
