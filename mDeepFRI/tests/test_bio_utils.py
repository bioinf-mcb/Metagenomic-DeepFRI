import os
import unittest

import requests
from biotite.structure.io.pdb import PDBFile

from mDeepFRI.bio_utils import get_residues_coordinates


class TestGetResiduesCoordinates(unittest.TestCase):
    def setUp(self) -> None:
        afdb_link = "https://alphafold.ebi.ac.uk/files/AF-F4HVG8-F1-model_v6.pdb"
        afdb_path = "./AF-A7YWM6-F1-model_v4.pdb"
        with open(afdb_path, "wb") as f:
            f.write(requests.get(afdb_link).content)
        self.afdb_structure = PDBFile.read(afdb_path).get_structure()[0]

    def test_predicted_default_chain(self):
        sequence, coordinates = get_residues_coordinates(self.afdb_structure)
        self.assertEqual(sequence[:10], "MLLSAIASQT")
        self.assertAlmostEqual(coordinates[:10].sum(), -25.926003, places=5)

    def test_predicted_invalid_chain(self):
        chain = "B"
        with self.assertRaises(ValueError):
            sequence, coordinates = get_residues_coordinates(
                self.afdb_structure, chain)

    def tearDown(self) -> None:
        os.remove("./AF-A7YWM6-F1-model_v4.pdb")
