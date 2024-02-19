import unittest

from mDeepFRI.bio_utils import insert_gaps


class TestBioUtils(unittest.TestCase):
    def test_insert_gaps(self):
        self.assertEqual(insert_gaps("ACGT", "ACGT", "MMMM")[0], "ACGT")
        self.assertEqual(insert_gaps("ACT", "ACCT", "MMIM")[0], "AC-T")
        self.assertEqual(insert_gaps("ACGT", "ACT", "MMDM")[1], "AC-T")
