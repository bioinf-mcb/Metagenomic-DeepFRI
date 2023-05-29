import os
import platform
import shutil
import tarfile
from distutils.util import convert_path
from pathlib import Path

import numpy as np
import requests
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

FOLDCOMP_BINARIES = {
    "LINUX": "https://mmseqs.com/foldcomp/foldcomp-linux-x86_64.tar.gz",
    "LINUX_ARM64 ": "https://mmseqs.com/foldcomp/foldcomp-linux-arm64.tar.gz",
    "MAC": "https://mmseqs.com/foldcomp/foldcomp-macos-universal.tar.gz",
    "WIN": "https://mmseqs.com/foldcomp/foldcomp-windows-x64.zip"
}

main_ns = {}
ver_path = convert_path('mDeepFRI/__init__.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def download_file(url, path):
    with requests.get(url, stream=True) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def download_foldcomp(output_path):
    output_path = Path(output_path)
    foldcomp_tar = output_path / "foldcomp.tar.gz"
    foldcomp_bin = output_path / "foldcomp"
    # check if foldcomp is downloaded
    if not foldcomp_bin.exists():
        # get OS
        info = platform.uname()
        system = info[0]
        arch = info[4]

        if system == "Linux":
            if arch == "ARM64":
                url = FOLDCOMP_BINARIES["LINUX_ARM64"]
            else:
                url = FOLDCOMP_BINARIES["LINUX"]
        elif system == "Windows":
            url = FOLDCOMP_BINARIES["WINDOWS"]
        elif system == "Darwin":
            url = FOLDCOMP_BINARIES["MACOS"]
        download_file(url, foldcomp_tar)
        # untar file
        with tarfile.open(foldcomp_tar, "r:gz") as archive:
            archive.extract("foldcomp", output_path)
        # remove tar file
        foldcomp_tar.unlink()


class build_ext(_build_ext):
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
        download_foldcomp(self.build_lib)


SRC_DIR = "mDeepFRI"
PACKAGES = [SRC_DIR]

install_requires = ["cython", "numpy"]
setup_requires = ["cython"]

EXTENSIONS = [
    Extension("mDeepFRI.atoms_io",
              sources=[SRC_DIR + "/atoms_io.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=["-std=c++17", "-O3"],
              extra_link_args=["-std=c++17"]),
    Extension("mDeepFRI.parsers",
              sources=[SRC_DIR + "/parsers.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=["-std=c++17", "-O3"]),
    Extension("mDeepFRI.predict",
              sources=[SRC_DIR + "/predict.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=["-std=c++17", "-O3"]),
    Extension("mDeepFRI.alignment_utils",
              sources=[SRC_DIR + "/alignment_utils.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=["-std=c++17", "-O3"]),
]

setup(
    name="mDeepFRI",
    version=main_ns['__version__'],
    description=
    "Pipeline for searching and aligning contact maps for proteins, then running DeepFri's GCN.",
    long_description=read("README.md"),
    long_description_content_type='text/markdown',
    keywords="protein function metagenomics deep neural network",
    author="Piotr Kucharski, Valentyn Bezshapkin",
    author_email=
    "soliareofastorauj@gmail.com, valentyn.bezshapkin@micro.biol.ethz.ch",
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
    install_requires=install_requires,
    cmdclass={'build_ext': build_ext},
    license="GNU GPLv3",
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
