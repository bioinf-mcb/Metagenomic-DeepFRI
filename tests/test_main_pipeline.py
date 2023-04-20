import subprocess

import pytest

from meta_deepFRI.utils import utils


def test_deepfri_errors():
    with pytest.raises(RuntimeError):
        result = utils.run_command("deepfri").strip()
