import unittest
from tempfile import NamedTemporaryFile

from mDeepFRI.utils import (load_fasta_as_dict, retrieve_fasta_entries_as_dict,
                            run_command)


class TestCommands(unittest.TestCase):
    def test_run_command(self):
        command = "echo 'Hello World!'"
        result = run_command(command).strip()
        self.assertEqual(result, "Hello World!")


class TestFastaFunctions(unittest.TestCase):
    def setUp(self):
        self.fasta_content = ">seq1\nATGC\n>seq2\nATGCATGC\n>seq3\nATGCATGCATGC\n"
        self.fasta_file = NamedTemporaryFile(suffix=".fasta", mode="w")
        self.fasta_file.write(self.fasta_content)
        self.fasta_file.flush()

        self.gzip_fasta_content = self.fasta_content.encode()
        self.gzip_fasta_file = NamedTemporaryFile(suffix=".fasta.gz",
                                                  mode="wb")
        self.gzip_fasta_file.write(self.gzip_fasta_content)
        self.gzip_fasta_file.flush()

    def test_load_fasta_as_dict(self):
        fasta_dict = load_fasta_as_dict(self.fasta_file.name)
        self.assertEqual(len(fasta_dict), 3)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertIn("seq3", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")
        self.assertEqual(fasta_dict["seq3"], "ATGCATGCATGC")

        # test with gzip file
        fasta_dict = load_fasta_as_dict(self.gzip_fasta_file.name)
        self.assertEqual(len(fasta_dict), 3)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertIn("seq3", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")
        self.assertEqual(fasta_dict["seq3"], "ATGCATGCATGC")

    def test_retrieve_fasta_entries_as_dict(self):
        entries = ["seq1", "seq2"]
        fasta_dict = retrieve_fasta_entries_as_dict(self.fasta_file.name,
                                                    entries)
        self.assertEqual(len(fasta_dict), 2)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")

        # test with gzip file
        fasta_dict = retrieve_fasta_entries_as_dict(self.gzip_fasta_file.name,
                                                    entries)
        self.assertEqual(len(fasta_dict), 2)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")

    def tearDown(self):
        self.fasta_file.close()
        self.gzip_fasta_file.close()
