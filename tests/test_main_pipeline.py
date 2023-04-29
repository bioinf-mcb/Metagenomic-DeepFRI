import pytest

from mDeepFRI.utils import utils


def test_deepfri_errors():
    with pytest.raises(RuntimeError):
        utils.run_command("mDeepFri").strip()
