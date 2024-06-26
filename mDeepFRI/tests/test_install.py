import subprocess
import unittest

from mDeepFRI.mmseqs import FOLDCOMP_PATH, MMSEQS_PATH


class TestSetup(unittest.TestCase):
    def test_foldcomp(self):
        """
        Test if foldcomp binary was downloaded properly
        and can be invoked from Python.
        """
        stdout = subprocess.run([FOLDCOMP_PATH],
                                capture_output=True,
                                check=True,
                                text=True).stdout
        assert "foldcomp compress" in stdout

    def test_mmseqs(self):
        """
        Test if mmseqs binary was downloaded properly
        and can be invoked from Python.
        """

        stdout = subprocess.run([MMSEQS_PATH],
                                capture_output=True,
                                check=True,
                                text=True).stdout
        assert "mmseqs" in stdout
