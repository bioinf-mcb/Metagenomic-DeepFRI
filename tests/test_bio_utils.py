import unittest

from mDeepFRI.alignment_utils import alignment_identity
from mDeepFRI.bio_utils import insert_gaps


class TestBioUtils(unittest.TestCase):
    def testInsertGaps(self):
        self.assertEqual(insert_gaps('AACT', 'AAT', 'MMDM'), ('AACT', 'AA-T'))
        self.assertEqual(insert_gaps('AAT', 'AATC', 'MMMI'), ('AAT-', 'AATC'))
        self.assertEqual(insert_gaps('AAT', 'FGTC', 'XXMI'), ('AAT-', 'FGTC'))

    def testIdentity(self):
        self.assertEqual(round(alignment_identity('AASDS', 'AASDS'), 2), 1)
        self.assertEqual(round(alignment_identity('AASDS', 'ADSDS'), 2), 0.8)
        self.assertEqual(round(alignment_identity('AASSS', 'A-S-S'), 2), 0.6)
        self.assertEqual(round(alignment_identity('A--SS', 'A-S--'), 2), 0.2)
