import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from mDeepFRI.mmseqs import MMseqsResult, QueryFile
from mDeepFRI.pdb import create_pdb_mmseqs


class TestQueryFile(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.fasta_file = os.path.join(self.temp_dir.name, "test.fasta")
        self.database = create_pdb_mmseqs()
        with open(self.fasta_file, "w") as f:
            f.write(">seq1\nATGC\n>seq2\nNNNNNCATGC\n>seq3\nATGCATGCATGC\n" \
                    ">seq4\nAAAAAAAAAACT")

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_init(self):
        query_file = QueryFile(self.fasta_file)
        self.assertEqual(query_file.filepath, self.fasta_file)
        self.assertEqual(query_file.sequences, {})

    def test_load_ids(self):
        query_file = QueryFile(self.fasta_file)
        query_file.load_ids(["seq1", "seq2"])
        self.assertEqual(len(query_file.sequences), 2)
        self.assertIn("seq1", query_file.sequences)
        self.assertIn("seq2", query_file.sequences)

    def test_load_sequences(self):
        query_file = QueryFile(self.fasta_file)
        query_file.load_sequences()
        self.assertEqual(len(query_file.sequences), 4)
        self.assertIn("seq1", query_file.sequences)
        self.assertIn("seq2", query_file.sequences)
        self.assertIn("seq3", query_file.sequences)
        self.assertIn("seq4", query_file.sequences)

    def test_remove_sequences(self):
        query_file = QueryFile(self.fasta_file)
        query_file.load_sequences()
        query_file.remove_sequences(["seq1", "seq2"])
        self.assertEqual(len(query_file.sequences), 2)
        self.assertNotIn("seq1", query_file.sequences)
        self.assertNotIn("seq2", query_file.sequences)
        self.assertIn("seq3", query_file.sequences)
        self.assertIn("seq4", query_file.sequences)

    def test_filter_sequences(self):
        query_file = QueryFile(self.fasta_file)
        query_file.load_sequences()
        query_file.filter_sequences(lambda seq: len(seq) <= 10)
        self.assertEqual(len(query_file.sequences), 2)
        self.assertIn("seq1", query_file.sequences)
        self.assertIn("seq2", query_file.sequences)

        query_file.filter_sequences(lambda seq: 'ATG' in seq)
        self.assertEqual(len(query_file.sequences), 2)
        self.assertIn("seq1", query_file.sequences)
        self.assertIn("seq2", query_file.sequences)

        query_file.filter_sequences(lambda seq: seq.count('N') < 5)
        self.assertEqual(len(query_file.sequences), 1)
        self.assertIn("seq1", query_file.sequences)
        self.assertNotIn("seq3", query_file.sequences)

    def test_search(self):
        query_file = QueryFile(self.fasta_file)
        query_file.load_sequences()
        result = query_file.search(self.database.mmseqs_db)
        self.assertIsInstance(result, MMseqsResult)
        self.assertIsNotNone(result.query_fasta)
        self.assertIsNotNone(result.database)


class TestMMseqsResult(unittest.TestCase):
    def setUp(self):
        self.data = np.rec.array([
            ("seq1", "target1", 10, 20, 0.3, 40, 1e-3),
            ("seq1", "target2", 20, 30, 0.4, 50, 1e-4),
            ("seq2", "target3", 30, 40, 0.5, 60, 1e-5),
            ("seq2", "target4", 40, 50, 0.6, 70, 1e-6),
            ("seq3", "target5", 50, 60, 0.7, 80, 1e-7),
        ],
                                 dtype=[("query", "U10"), ("target", "U10"),
                                        ("qcov", int), ("tcov", int),
                                        ("fident", float), ("bits", float),
                                        ("evalue", float)])
        self.query_fasta = Path("path/to/query.fasta").resolve()
        self.database = Path("path/to/database").resolve()
        self.result = MMseqsResult(self.data, self.query_fasta, self.database)

    def test_init(self):
        self.assertIsInstance(self.result, MMseqsResult)
        self.assertEqual(self.result.query_fasta,
                         Path(self.query_fasta).resolve())
        self.assertEqual(self.result.database, Path(self.database).resolve())

    def test_columns(self):
        np.array_equal(
            self.result.columns,
            np.array([
                "query", "target", "qcov", "tcov", "fident", "bits", "evalue"
            ]))

    def test_save_tsv(self):
        filepath = "test.tsv"
        self.result.save(filepath)
        with open(filepath, "r") as f:
            lines = f.readlines()
            self.assertEqual(
                lines[0],
                "query\ttarget\tqcov\ttcov\tfident\tbits\tevalue\tquery_file\tdatabase_file\n"
            )
            self.assertEqual(
                lines[1],
                f"seq1\ttarget1\t10\t20\t0.3\t40.0\t{float(1e-03)}\t{self.query_fasta}\t{self.database}\n"
            )

    def test_save_npz(self):
        filepath = "test.npz"
        self.result.save(filepath, filetype="npz")
        loaded_data = np.load(filepath)
        self.assertTrue(
            np.array_equal(loaded_data["arr_0"], self.result.result_arr))

    def test_apply_filters(self):
        filtered = self.result.apply_filters(min_cov=30,
                                             min_ident=0.4,
                                             min_bits=50)
        self.assertEqual(len(filtered), 3)

    def test_find_best_matches(self):
        best_matches = self.result.find_best_matches(k=2)
        self.assertEqual(len(best_matches), 5)

    def test_from_mmseqs_result(self):
        filepath = "test.tsv"
        with open(filepath, "w", newline="") as f:
            f.write("query\ttarget\tqcov\ttcov\tfident\tbits\tevalue\n")
            f.write(f"seq1\ttarget1\t10\t20\t0.3\t40.0\t{float(1e-03)}\n")
            f.write(f"seq1\ttarget2\t20\t30\t0.4\t50.0\t{float(1e-04)}\n")
        result = MMseqsResult.from_mmseqs_result(filepath, self.query_fasta,
                                                 self.database)
        self.assertIsInstance(result, MMseqsResult)
        self.assertEqual(result.query_fasta, Path(self.query_fasta).resolve())
        self.assertEqual(result.database, Path(self.database).resolve())


if __name__ == "__main__":
    unittest.main()
