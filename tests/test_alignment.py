import unittest

import pyopal

from mDeepFRI.alignment import align_pairwise, best_hit_database


class TestAlignment(unittest.TestCase):
    def setUp(self) -> None:
        self.queries = dict(
            query_seq="MAGFLKVVQLLAKYGSKAVQWAWANKGKILDWLNAGQAIDWVVS")
        self.targets = dict(
            seq1="MESILDLQELETSEEESALMAASTVSNNC",
            seq2="MKKAVIVENKGCATCSIGAACLVDGPIPDFEIAGATGLFGLWG",
            seq3="MAGFLKVVQILAKYGSKAVQWAWANKGKILDWINAGQAIDWVVE",
            seq4="MAGFLKVVQILAKYGSKAVQWAWANKGKILDWINAGQAIDWVVE",
        )
        self.aligner = pyopal.Aligner()
        return

    def test_best_hit_database(self):
        database = pyopal.Database(self.targets.values())

        best_index = best_hit_database(self.queries["query_seq"], database,
                                       self.aligner)
        self.assertEqual(best_index, 2)

    def test_align_pairwise(self):
        alignment = align_pairwise(self.queries["query_seq"],
                                   self.targets["seq3"], self.aligner)
        self.assertEqual(alignment,
                         "MMMMMMMMMXMMMMMMMMMMMMMMMMMMMMMMXMMMMMMMMMMX")
