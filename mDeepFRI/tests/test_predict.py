import unittest

import numpy as np

from mDeepFRI.predict import seq2onehot


class TestSeq2OneHot(unittest.TestCase):
    def test_empty_sequence(self):
        seq = ""
        result = seq2onehot(seq)
        self.assertEqual(result.shape, (0, 26))
        self.assertEqual(result.dtype, np.float32)
        self.assertTrue(np.all(result == np.array([]).reshape(0, 26)))

    def test_single_character_sequence(self):
        seq = "D"
        result = seq2onehot(seq)
        self.assertEqual(result.shape, (1, 26))
        self.assertTrue(np.all(result == np.array([[0, 1] + [0] * 24])))

    def test_sequence(self):
        seq = "-DGU"
        result = seq2onehot(seq)
        self.assertEqual(result.shape, (4, 26))
        expected = np.array([[1, 0, 0, 0] + [0] * 22, [0, 1, 0, 0] + [0] * 22,
                             [0, 0, 1, 0] + [0] * 22, [0, 0, 0, 1] + [0] * 22])
        self.assertTrue(np.all(result == expected))

    def test_invalid_character(self):
        seq = "J"
        with self.assertRaises(ValueError):
            seq2onehot(seq)
