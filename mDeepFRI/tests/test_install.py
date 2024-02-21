import subprocess
from pathlib import Path

import mDeepFRI


# is a tricky test, will fail if ran from a wrong directory
def test_foldcomp_download():
    """
    Test if foldcomp binary was downloaded properly
    and can be invoked from Python.
    """
    foldcomp_bin = Path(mDeepFRI.__path__[0]).parent / "foldcomp"
    stdout = subprocess.run([foldcomp_bin],
                            capture_output=True,
                            check=True,
                            text=True).stdout
    assert "foldcomp compress" in stdout
