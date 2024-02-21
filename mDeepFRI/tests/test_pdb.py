import pathlib
import unittest

from mDeepFRI.pdb import create_pdb_mmseqs


class TestPDB(unittest.TestCase):
    def test_create_pdb_mmseqs(self):
        pdb100 = create_pdb_mmseqs()

        self.assertTrue(pathlib.Path(pdb100.sequence_db).exists())
        self.assertTrue(pathlib.Path(pdb100.mmseqs_db).exists())
