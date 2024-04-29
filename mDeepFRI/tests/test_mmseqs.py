import gzip
import unittest
from pathlib import Path
from unittest import mock

from mDeepFRI.mmseqs import QueryFile


class TestQueryFile(unittest.TestCase):
    def setUp(self):
        self.filepath = "mDeepFRI/tests/data/test.fa"
        # create empty file
        Path(self.filepath).touch()
        self.query_file = QueryFile(self.filepath)

    def test_load_sequences(self):
        # Mock the FastxFile
        with mock.patch("mDeepFRI.mmseqs.FastxFile") as mock_fastx_file:
            mock_entry1 = mock.MagicMock()
            mock_entry1.name = "seq1"
            mock_entry1.sequence = "ATCG"
            mock_entry2 = mock.MagicMock()
            mock_entry2.name = "seq2"
            mock_entry2.sequence = "CGTA"
            mock_fastx_file.return_value.__enter__.return_value = [
                mock_entry1, mock_entry2
            ]

            self.query_file.load_sequences()
            expected_sequences = {"seq1": "ATCG", "seq2": "CGTA"}
            self.assertEqual(self.query_file.sequences, expected_sequences)

    def test_load_ids(self):
        # Mock the FastaFile
        with mock.patch("mDeepFRI.mmseqs.FastaFile") as mock_fasta_file:
            mock_fasta_file.return_value.fetch.side_effect = ["ATCG", "CGTA"]
            self.query_file.load_ids(["seq1", "seq2"])
            expected_sequences = {"seq1": "ATCG", "seq2": "CGTA"}
            self.assertEqual(self.query_file.sequences, expected_sequences)

    def test_load_ids_gzip(self):
        # gzip a filepath
        gzip_filepath = self.filepath + ".gz"
        with open(self.filepath, "rb") as f_in:
            with gzip.open(gzip_filepath, "wb") as f_out:
                f_out.writelines(f_in)
        self.query_file.filepath = gzip_filepath

        # Mock the FastaFile
        with mock.patch("mDeepFRI.mmseqs.FastaFile") as mock_fasta_file:
            mock_fasta_file.return_value.fetch.side_effect = ["ATCG", "CGTA"]
            self.query_file.load_ids(["seq1", "seq2"])
            expected_sequences = {"seq1": "ATCG", "seq2": "CGTA"}
            self.assertEqual(self.query_file.sequences, expected_sequences)

        # reset filepath
        self.query_file.filepath = self.filepath
        # delete gzip file
        Path(gzip_filepath).unlink()

    def test_load_ids_not_found(self):
        # Test for FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            self.query_file.filepath = "path/to/non-existent.fasta"
            self.query_file.load_ids(["seq1", "seq2"])

    def test_load_ids_invalid_id(self):
        # Test for ValueError
        with mock.patch("mDeepFRI.mmseqs.FastaFile") as mock_fasta_file:
            mock_fasta_file.return_value.fetch.side_effect = KeyError
            with self.assertRaises(ValueError):
                self.query_file.load_ids(["seq1", "seq2"])

    def test_filter_sequences_min_length(self):
        self.query_file.sequences = {"seq1": "ATCG", "seq2": "CGTAY"}

        # Test min_length
        self.query_file.filter_sequences(min_length=5)
        expected_sequences = {"seq2": "CGTAY"}
        self.assertEqual(self.query_file.sequences, expected_sequences)
        self.assertEqual(self.query_file.too_long, ["seq1"])

    def test_filter_sequences_max_length(self):
        self.query_file.sequences = {"seq1": "ATCG", "seq2": "CGTAY"}
        # Test max_length
        self.query_file.filter_sequences(max_length=4)
        expected_sequences = {"seq1": "ATCG"}
        self.assertEqual(self.query_file.sequences, expected_sequences)
        self.assertEqual(self.query_file.too_short, ["seq2"])

    def test_filter_sequences_empty(self):
        # Test empty sequences
        self.query_file.sequences = {}
        with self.assertRaises(ValueError):
            self.query_file.filter_sequences(min_length=5)

    # @mock.patch("mDeepFRI.mmseqs._createdb")
    # @mock.patch("mDeepFRI.mmseqs._createindex")
    # @mock.patch("mDeepFRI.mmseqs._search")
    # @mock.patch("mDeepFRI.mmseqs._convertalis")
    # def test_search(self, mock_convertalis, mock_search, mock_createindex, mock_createdb):
    #     self.query_file.sequences = {"seq1": "ATCG", "seq2": "CGTA"}

    #     # Test with valid sensitivity value
    #     self.query_file.search("path/to/database", sensitivity=5.0)

    #     # Test with invalid sensitivity value
    #     with self.assertRaises(ValueError):
    #         self.query_file.search("path/to/database", sensitivity=8.0)

    #     # Test with empty sequences
    #     self.query_file.sequences = {}
    #     with mock.mock_open(self.filepath, "r"):
    #         self.query_file.search("path/to/database")

    def tearDown(self):
        Path(self.filepath).unlink()


if __name__ == "__main__":
    unittest.main()
