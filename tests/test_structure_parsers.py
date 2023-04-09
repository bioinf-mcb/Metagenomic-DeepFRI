import numpy as np
import gzip
import pathlib
import pytest

from meta_deepFRI.structure_files.parse_structure_file import SeqAtoms, read_structure_file


@pytest.mark.parametrize(
    "file_name, excepted_id, excepted_amino_group, excepted_position, excepted_group, excepted_length",
    [("6a0j.cif.gz", "6a0j", ['VAL', 'VAL', 'VAL', 'VAL'], [-6.123, 31.928, 92.243], ['A22'], 3443),
     ("AF-A0A2Z5TJB0-F1-model_v4.pdb.gz", "AF-A0A2Z5TJB0-F1-model_v4", ['MET', 'MET', 'MET', 'MET'], [-47.025, 49.762, -23.074], ['A   1'], 3956)])  # yapf: disable
def test_proper_structure_files(file_name, excepted_id, excepted_amino_group, excepted_position, excepted_group,
                                excepted_length):
    file_path = pathlib.Path("tests") / "data" / "structures" / file_name
    output: SeqAtoms = read_structure_file(file_path)

    assert output.protein_id == excepted_id
    assert (output.atom_amino_group[:4] == excepted_amino_group).all()
    assert (output.positions[:1] == np.array([excepted_position], dtype=np.float32)).all()
    assert output.groups[:1] == excepted_group
    assert len(output.atom_amino_group) == len(output.positions) == len(output.groups) == excepted_length
