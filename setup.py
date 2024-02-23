import configparser
import os
import re
import shutil
import tarfile
from distutils.util import convert_path
from pathlib import Path

import numpy as np
import requests
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Utils -----------------------------------------------------------------


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


FOLDCOMP_BINARIES = {
    "linux_or_android":
    "https://mmseqs.com/foldcomp/foldcomp-linux-x86_64.tar.gz",
    "aarch64": "https://mmseqs.com/foldcomp/foldcomp-linux-arm64.tar.gz",
    "macos": "https://mmseqs.com/foldcomp/foldcomp-macos-universal.tar.gz",
    "windows": "https://mmseqs.com/foldcomp/foldcomp-windows-x64.zip"
}

main_ns = {}
ver_path = convert_path('mDeepFRI/__init__.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)


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


class build_ext(_build_ext):
    def initialize_options(self) -> None:
        _build_ext.initialize_options(self)
        self.target_machine = None
        self.target_cpu = None
        self.target_system = None

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
        _download_foldcomp(self.target_system, self.target_cpu, self.build_lib)


# --- Commands ------------------------------------------------------------------


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
extra_compile_args = ["-std=c++14", "-O3"]

EXTENSIONS = [
    Extension("mDeepFRI.predict",
              sources=[SRC_DIR + "/predict.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=extra_compile_args),
    Extension("mDeepFRI.alignment_utils",
              sources=[SRC_DIR + "/alignment_utils.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=extra_compile_args),
]

extras = {}
extras["dev"] = ["pre-commit"]

setup(
    name="mDeepFRI",
    version=main_ns['__version__'],
    description=
    "Pipeline for searching and aligning contact maps for proteins, then running DeepFri's GCN.",
    long_description=read("README.md"),
    long_description_content_type='text/markdown',
    keywords="protein function metagenomics deep neural network",
    author=" Valentyn Bezshapkin, Piotr Kucharski",
    author_email=
    "valentyn.bezshapkin@micro.biol.ethz.ch, soliareofastorauj@gmail.com",
    url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    download_url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "mDeepFRI = mDeepFRI.cli:main",
        ],
    },
    ext_modules=EXTENSIONS,
    include_dirs=[np.get_include()],
    extras_require=extras,
    install_requires=install_requires,
    cmdclass={'build_ext': build_ext},
    license="BSD-3-Clause",
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
)
