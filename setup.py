import os
from distutils.util import convert_path

import numpy as np
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

main_ns = {}
ver_path = convert_path('mDeepFRI/__init__.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


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


SRC_DIR = "mDeepFRI"
PACKAGES = [SRC_DIR]

install_requires = ["cython"]

EXTENSIONS = [
    Extension("mDeepFRI.CPP_lib.atoms_io",
              sources=[
                  SRC_DIR + "/CPP_lib/atoms_io.pyx",
                  SRC_DIR + "/CPP_lib/atoms_file_io.cpp",
                  SRC_DIR + "/CPP_lib/load_contact_maps.cpp"
              ],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=["-std=c++17", "-O3"],
              extra_link_args=["-std=c++17"]),
    Extension("mDeepFRI.CPP_lib.parsers",
              sources=[SRC_DIR + "/CPP_lib/parsers.pyx"],
              language="c++",
              libraries=["stdc++"],
              extra_compile_args=["-std=c++17", "-O3"])
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
    cmdclass={
        'build_ext': build_ext,
    },
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
