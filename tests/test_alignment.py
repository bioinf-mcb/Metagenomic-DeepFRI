import unittest

import pyopal

from mDeepFRI.alignment import best_hit_database


class TestAlignment(unittest.TestCase):
    def test_best_hit_database(self):
        queries = {
            "query_seq":
            "MAGFLKVVQLLAKYGSKAVQWAWANKGKILDWLNAGQAIDWVVSKIKQILGIK"
        }
        targets = dict(
            seq1="MESILDLQELETSEEESALMAASTVSNNC",
            seq2="MKKAVIVENKGCATCSIGAACLVDGPIPDFEIAGATGLFGLWG",
            seq3="MAGFLKVVQILAKYGSKAVQWAWANKGKILDWINAGQAIDWVVEKIKQILGIK",
            seq4="MAGFLKVVQILAKYGSKAVQWAWANKGKILDWINAGQAIDWVVEKIKQILGIK",
        )
        aligner = pyopal.Aligner()
        database = pyopal.Database(targets.values())

        best_index = best_hit_database(queries["query_seq"], database, aligner)
        self.assertEqual(best_index, 2)
