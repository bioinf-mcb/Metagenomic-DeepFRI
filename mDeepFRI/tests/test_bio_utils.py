import os
import unittest

import numpy as np
import requests
from biotite.structure.io.pdb import PDBFile

from mDeepFRI.bio_utils import get_residues_coordinates
from mDeepFRI.contact_map_utils import pairwise_sqeuclidean


class TestPairwiseSqeuclidean(unittest.TestCase):
    def test_pairwise_sqeuclidean(self):
        np.random.seed(42)
        expected = np.array(
            [[0, 1.01354558, 0.12442072], [1.01354558, 0, 0.99467713],
             [0.12442072, 0.99467713, 0]],
            dtype=np.float32)

        matrix = np.random.rand(3, 3).astype(np.float32)
        result = pairwise_sqeuclidean(matrix)
        np.allclose(result, expected)


class TestGetResiduesCoordinates(unittest.TestCase):
    def setUp(self) -> None:
        afdb_link = "https://alphafold.ebi.ac.uk/files/AF-A7YWM6-F1-model_v4.pdb"
        afdb_path = "./AF-A7YWM6-F1-model_v4.pdb"
        with open(afdb_path, "wb") as f:
            f.write(requests.get(afdb_link).content)
        self.afdb_structure = PDBFile.read(afdb_path).get_structure()[0]

    def test_predicted_default_chain(self):
        sequence, coordinates = get_residues_coordinates(self.afdb_structure)
        self.assertEqual(sequence[:10], "MAPLVAQLLF")
        self.assertAlmostEqual(coordinates[:10].sum(), 199.20398, places=5)

    def test_predicted_invalid_chain(self):
        chain = "B"
        with self.assertRaises(ValueError):
            sequence, coordinates = get_residues_coordinates(
                self.afdb_structure, chain)

    def tearDown(self) -> None:
        os.remove("./AF-A7YWM6-F1-model_v4.pdb")
