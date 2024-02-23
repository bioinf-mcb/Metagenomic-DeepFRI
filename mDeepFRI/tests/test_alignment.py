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

        return super().setUp()

    def test_best_hit_database(self):

        best_hit, _ = best_hit_database(self.queries["query_seq"],
                                        self.targets)
        self.assertEqual(best_hit, "seq3")

    def test_align_pairwise(self):
        alignment = align_pairwise(self.queries["query_seq"],
                                   self.targets["seq3"])
        self.assertEqual(alignment,
                         "MMMMMMMMMXMMMMMMMMMMMMMMMMMMMMMMXMMMMMMMMMMX")
