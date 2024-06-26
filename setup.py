import configparser
import os
import re
import shutil
import tarfile
from pathlib import Path

import requests
from setuptools import Command, Extension, find_packages, setup
from setuptools.command.build import build as _build
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Utils -----------------------------------------------------------------

os.environ["CC"] = "gcc"
os.environ["CXX"] = "g++"


def _detect_target_machine(platform):
    if platform == "win32":
        return "x86"
    return platform.rsplit("-", 1)[-1]


def _detect_target_cpu(platform):
    machine = _detect_target_machine(platform)
    if re.match("^mips", machine):
        return "mips"
    elif re.match("^(aarch64|arm64)$", machine):
        return "aarch64"
    elif re.match("^arm", machine):
        return "arm"
    elif re.match("(x86_64)|(x86)|(AMD64|amd64)|(^i.86$)", machine):
        return "x86"
    elif re.match("^(powerpc|ppc)", machine):
        return "ppc"
    return None


def _detect_target_system(platform):
    if platform.startswith("win"):
        return "windows"
    elif platform.startswith("macos"):
        return "macos"
    elif platform.startswith("linux"):
        return "linux_or_android"
    elif platform.startswith("freebsd"):
        return "freebsd"
    return None


def _detect_cpu_features():
    import archspec.cpu
    host = archspec.cpu.host()
    features = host.to_dict()["features"]
    if "avx2" in features:
        return "avx2"
    elif "sse4_1" in features:
        return "sse41"
    else:
        return "sse2"


FOLDCOMP_BINARIES = {
    "linux_or_android":
    "https://mmseqs.com/foldcomp/foldcomp-linux-x86_64.tar.gz",
    "aarch64": "https://mmseqs.com/foldcomp/foldcomp-linux-arm64.tar.gz",
    "macos": "https://mmseqs.com/foldcomp/foldcomp-macos-universal.tar.gz"
}

MMSEQS_BINARIES = {
    "avx2": "https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz",
    "sse41": "https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz",
    "sse2": "https://mmseqs.com/latest/mmseqs-linux-sse2.tar.gz",
    "aarch64": "https://mmseqs.com/latest/mmseqs-linux-arm64.tar.gz",
    "ppc": "https://mmseqs.com/latest/mmseqs-linux-ppc64le-power8.tar.gz"
}


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def _download_file(url, path):
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def _download_foldcomp(system, cpu, output_path):
    # select appropriate binary
    if cpu == "aarch64" and system == "linux":
        build = "aarch64"
    else:
        build = system

    url = FOLDCOMP_BINARIES[build]

    output_path = Path(output_path)
    foldcomp_tar = output_path / "foldcomp.tar.gz"

    _download_file(url, foldcomp_tar)
    # untar file
    with tarfile.open(foldcomp_tar, "r:gz") as archive:
        archive.extract("foldcomp", output_path)
    # rename binary after extraction to foldcomp_bin
    (output_path / "foldcomp").rename(output_path / "foldcomp_bin")

    # remove tar file
    foldcomp_tar.unlink()


def _download_mmseqs(system, cpu, features, output_path):
    # select appropriate binary
    if cpu == "aarch64":
        build = "aarch64"
    elif cpu == "ppc" and system == "linux":
        build = "ppc"
    else:
        build = features

    url = MMSEQS_BINARIES[build]

    output_path = Path(output_path)
    mmseqs_tar = output_path / "mmseqs.tar.gz"
    # download file
    _download_file(url, mmseqs_tar)
    # untar file
    with tarfile.open(mmseqs_tar, "r:gz") as archive:
        archive.extractall(path=output_path)

    # remove tar file
    mmseqs_tar.unlink()


# --- Commands ------------------------------------------------------------------


class build_binaries(Command):
    """
    A custom command to download foldcomp and mmseqs binaries.
    """

    description = "Download foldcomp and mmseqs binaries."
    user_options = [
        ("force", "f",
         "Force download even if the binaries are already present."),
        ("inplace", "i", "Download the binaries in the source directory."),
    ]

    def initialize_options(self):
        self.force = False
        self.inplace = False

    def finalize_options(self):
        _build_py = self.get_finalized_command("build_py")
        self.build_lib = _build_py.build_lib

    def run(self):
        if self.inplace:
            output_path = Path(".")
        else:
            output_path = Path(self.build_lib)

        if not self.force:
            if (output_path / "foldcomp_bin").exists() and (
                    output_path / "mmseqs_bin").exists():
                return

        _build_ext = self.get_finalized_command("build_ext")
        _build_ext.run()

        self.announce("Downloading foldcomp", level=5)
        _download_foldcomp(_build_ext.target_system, _build_ext.target_cpu,
                           output_path)

        self.announce("Downloading mmseqs", level=5)
        _download_mmseqs(_build_ext.target_system, _build_ext.target_cpu,
                         _build_ext.target_features, output_path)


class build_ext(_build_ext):
    def initialize_options(self) -> None:
        _build_ext.initialize_options(self)
        self.target_machine = None
        self.target_cpu = None
        self.target_system = None
        self.target_features = None

    def finalize_options(self) -> None:
        _build_ext.finalize_options(self)

        import builtins
        builtins.__NUMPY_SETUP__ = False

        import numpy
        self.include_dirs.append(numpy.get_include())

    def run(self):
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Failed to import Cython") from cythonize

        self.extensions = cythonize(self.extensions,
                                    compiler_directives={
                                        "linetrace": True,
                                        "language_level": 3
                                    })
        for ext in self.extensions:  # this fixes a bug with setuptools
            ext._needs_stub = False
        _build_ext.run(self)
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_features = _detect_cpu_features()

        if self.target_system == "windows":
            raise RuntimeError("Windows is not supported.")


class build(_build):
    def run(self):
        _build.run(self)
        _build_binaries = self.get_finalized_command("build_binaries")
        _build_binaries.run()
        _build.run(self)


class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """
    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires",
              str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


SRC_DIR = "mDeepFRI"
PACKAGES = [SRC_DIR]

install_requires = ["cython", "numpy", "requests"]
setup_requires = ["cython", "requests", "numpy"]
extra_compile_args = ["-O3"]

EXTENSIONS = [
    Extension("mDeepFRI.predict",
              sources=[SRC_DIR + "/predict.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=extra_compile_args,
              define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")
                             ]),
    Extension("mDeepFRI.contact_map_utils",
              sources=[SRC_DIR + "/contact_map_utils.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=extra_compile_args,
              define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")
                             ]),
]

extras = {}
extras["dev"] = ["pre-commit"]

setup(
    include_package_data=True,
    ext_modules=EXTENSIONS,
    extras_require=extras,
    install_requires=install_requires,
    cmdclass={
        'build_ext': build_ext,
        'build_binaries': build_binaries,
        'build': build,
        'sdist': sdist
    },
    packages=find_packages(),
    package_data={
        "foldcomp": ["/foldcomp_bin"],
        "mmseqs": ["/mmseqs"]
    },
)
