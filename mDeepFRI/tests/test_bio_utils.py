import unittest

import numpy as np

from mDeepFRI.alignment_utils import alignment_identity, pairwise_sqeuclidean
from mDeepFRI.bio_utils import insert_gaps


class TestBioUtils(unittest.TestCase):
    def test_insert_gaps(self):
        self.assertEqual(insert_gaps('AACT', 'AAT', 'MMDM'), ('AACT', 'AA-T'))
        self.assertEqual(insert_gaps('AAT', 'AATC', 'MMMI'), ('AAT-', 'AATC'))
        self.assertEqual(insert_gaps('AAT', 'FGTC', 'XXMI'), ('AAT-', 'FGTC'))

    def test_identity(self):
        self.assertEqual(round(alignment_identity('AASDS', 'AASDS'), 2), 1)
        self.assertEqual(round(alignment_identity('AASDS', 'ADSDS'), 2), 0.8)
        self.assertEqual(round(alignment_identity('AASSS', 'A-S-S'), 2), 0.6)
        self.assertEqual(round(alignment_identity('A--SS', 'A-S--'), 2), 0.2)

    def test_pairwise_sqeuclidean(self):
        np.random.seed(42)
        expected = np.array(
            [[0, 1.01354558, 0.12442072], [1.01354558, 0, 0.99467713],
             [0.12442072, 0.99467713, 0]],
            dtype=np.float32)

        matrix = np.random.rand(3, 3).astype(np.float32)
        result = pairwise_sqeuclidean(matrix)
        np.allclose(result, expected)
