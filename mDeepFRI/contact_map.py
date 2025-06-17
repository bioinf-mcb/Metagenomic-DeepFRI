import numpy as np

from mDeepFRI.contact_map_utils import pairwise_sqeuclidean


class CAlphaCoordinates:
    def __init__(self, structure_id: str, coords: np.ndarray):

        self.structure_id = structure_id
        self.coords = coords
        # check that coords are 3D
        if coords.shape[1] != 3:
            raise ValueError("Coordinates are not 3D.")

    def calculate_distance_map(self, distance="sqeuclidean"):
        """
        Calculate distance map from C-alpha coordinates.

        Args:
            distance (str): Distance metric.
        Returns:
            np.ndarray: Distance map.
        """
        if distance == "sqeuclidean":
            distances = pairwise_sqeuclidean(self.coords.astype(np.float32))
        else:
            raise NotImplementedError("Distance metric not implemented.")

        return DistanceMap(distances)

    def calculate_contact_map(self, threshold=6.0):
        """
        Calculate contact map from PDB string.

        Args:
            pdb_string (str): PDB file read into string.
            max_seq_len (int): Maximum sequence length.
            threshold (float): Distance threshold for contact map.
            mode (str): Output mode. Either "matrix" or "sparse".

        Returns:
            np.ndarray: Contact map.
        """

        distance_map = self.calculate_distance_map()
        cmap = distance_map.calculate_contacts(threshold**2)

        return cmap


class DistanceMap:
    def __init__(self, distance_map):
        self.distance_map = distance_map
        # check if values are positive
        if not np.all(distance_map >= 0):
            raise ValueError("Distance map contains negative values.")
        # check if diagonal is zero
        if not np.all(np.diag(distance_map) == 0):
            raise ValueError("Distance map diagonal is not zero.")
        # check if values are symmetric
        if not np.allclose(distance_map, distance_map.T):
            raise ValueError("Distance map is not symmetric.")

    def calculate_contacts(self, threshold: int):
        """
        Calculate contacts from distance map.

        Args:
            threshold (int): Distance threshold for contact.

        Returns:
            np.ndarray: Contact map.
        """
        cmap = (self.distance_map < threshold).astype(np.int32)
        return ContactMap(cmap)


class ContactMap:
    def __init__(self, cmap):
        self.cmap = cmap
        # check if cmap is symmetric
        if not np.allclose(cmap, cmap.T):
            raise ValueError("Contact map is not symmetric.")
        # check if values in range [0, 1]
        if not np.all(np.isin(cmap, [0, 1])):
            raise ValueError("Contact map values not in range [0, 1].")

    def sparsify(self):
        """
        Convert contact map to sparse format.

        Returns:
            np.ndarray: Sparse contact map.
        """
        return np.argwhere(self.cmap == 1).astype(np.int32)
