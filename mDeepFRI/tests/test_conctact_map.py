import unittest

import numpy as np

from mDeepFRI.contact_map import CAlphaCoordinates


class TestCAlphaCoords(unittest.TestCase):
    def test_init_valid_coords(self):
        coords = np.array([[1, 2, 3], [4, 5, 6]])
        alpha_coords = CAlphaCoordinates(coords)
        assert np.array_equal(alpha_coords.coords, coords)

    def test_init_invalid_dims(self):
        coords = np.array([[1, 2], [3, 4]])
        with self.assertRaisesRegex(ValueError, "Coordinates are not 3D."):
            CAlphaCoordinates(coords)

    def test_init_negative_coords(self):
        coords = np.array([[-1, 2, 3], [4, 5, 6]])
        with self.assertRaisesRegex(ValueError,
                                    "Coordinates contain negative values."):
            CAlphaCoordinates(coords)

    def test_calculate_distance_map(self):
        coords = np.array([[0, 0, 0], [1, 1, 1]])
        alpha_coords = CAlphaCoordinates(coords)
        distance_map = alpha_coords.calculate_distance_map(
            distance="sqeuclidean")
        expected_distances = np.array([[0, 3], [3, 0]], dtype=np.float32)
        assert np.allclose(np.sqrt(distance_map.distance_map),
                           expected_distances)

    def test_calculate_distance_map_invalid_function(self):
        coords = np.array([[0, 0, 0], [1, 1, 1]])
        alpha_coords = CAlphaCoordinates(coords)
        with self.assertRaisesRegex(NotImplementedError,
                                    "Distance metric not implemented."):
            alpha_coords.calculate_distance_map(distance="euclidean")

    def test_calculate_contact_map(self):
        coords = np.array([[0, 0, 0], [5, 0, 0], [10, 0, 0]])
        alpha_coords = CAlphaCoordinates(coords)
        contact_map = alpha_coords.calculate_contact_map(threshold=6.0)
        expected_contact_map = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]])
        assert np.array_equal(contact_map.cmap, expected_contact_map)
