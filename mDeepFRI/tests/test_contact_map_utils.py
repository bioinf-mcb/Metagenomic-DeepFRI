import unittest

import numpy as np

try:
    from mDeepFRI.contact_map_utils import (align_contact_map,
                                            pairwise_sqeuclidean)
except ImportError:
    raise ImportError(
        "Failed to import align_contact_map from mDeepFRI.contact_map_utils."
        "Run python setup.py build_ext --inplace to compile Cython extensions."
    )


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


class SetContactMapUtilsTestCase(unittest.TestCase):
    def test_identity_alignment(self):
        """Test simple 1:1 alignment with no gaps."""
        q_align = "AB"
        t_align = "AB"
        # Target has contact between 0 and 1
        target_contacts = np.array([[0, 1]], dtype=np.int32)

        result = align_contact_map(q_align, t_align, target_contacts)

        expected = np.array([[1, 1], [1, 1]], dtype=np.int32)

        np.testing.assert_array_equal(result, expected)

    def test_gap_in_query_deletion(self):
        """
        Query:  A-C (Deletion of residue B)
        Target: ABC

        Target Contacts:
        (0,1) -> A-B -> Should be removed (B is gap in query)
        (1,2) -> B-C -> Should be removed (B is gap in query)
        (0,2) -> A-C -> Should be kept (mapped to 0-1 in query)
        """
        q_align = "A-C"
        t_align = "ABC"
        target_contacts = np.array([[0, 1], [1, 2], [0, 2]], dtype=np.int32)

        result = align_contact_map(q_align, t_align, target_contacts)

        # Result should be 2x2 (A, C)
        # We expect a contact between A(0) and C(1)
        expected = np.array([[1, 1], [1, 1]], dtype=np.int32)

        np.testing.assert_array_equal(result, expected)

    def test_gap_in_target_insertion(self):
        """
        Query:  ABC (Insertion of B)
        Target: A-C

        Indices:
        Target: A=0, C=1
        Query:  A=0, B=1, C=2

        Target Contacts: (0, 1) -> A-C

        Generated contacts (gen=1) at B(1):
        - Contact (1+1, 1) -> (2, 1) -> C-B
        - Contact (1-1, 1) -> (0, 1) -> A-B
        """
        q_align = "ABC"
        t_align = "A-C"
        target_contacts = np.array([[0, 1]], dtype=np.int32)

        # generated_contacts=1 (neighbors only)
        result = align_contact_map(q_align,
                                   t_align,
                                   target_contacts,
                                   generated_contacts=1)

        expected = np.array(
            [
                [1, 1, 1],  # 0-0, 0-1(gen), 0-2(mapped A-C)
                [1, 1, 1],  # 1-0(gen), 1-1, 1-2(gen)
                [1, 1, 1]  # 2-0(mapped), 2-1(gen), 2-2
            ],
            dtype=np.int32)

        np.testing.assert_array_equal(result, expected)

    def test_large_input_stress(self):
        """Simple stress test to check for segmentation faults or memory errors."""
        N = 100
        q_align = "A" * N
        t_align = "A" * N
        # Create a chain of contacts
        target_contacts = np.array([[i, i + 1] for i in range(N - 1)],
                                   dtype=np.int32)

        result = align_contact_map(q_align, t_align, target_contacts)
        self.assertEqual(result.shape, (N, N))
        self.assertEqual(result[0, 1], 1)
