import unittest

from mDeepFRI import utils


class TestUtils(unittest.TestCase):
    def test_run_command(self):
        command = "echo 'Hello World!'"
        result = utils.run_command(command).strip()
        self.assertEqual(result, "Hello World!")

    def test_mmseqs(self):
        command = "mmseqs"
        result = utils.run_command(command).strip()
        self.assertEqual(result.split()[0], "MMseqs2")
