import os

import numpy as np
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext

from mDeepFRI import __version__

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


# class CMakeExtension(Extension):

#     def __init__(self, name):
#         # don't invoke the original build_ext for this special extension
#         super().__init__(name, sources=[])

# class build_ext(_build_ext):

#     def run(self):
#         for ext in self.extensions:
#             self.build_cmake(ext)
#         super().run()

#     def build_cmake(self, ext):
#         cwd = pathlib.Path().absolute()

#         # these dirs will be created in build_py, so if you don't have
#         # any python sources to bundle, the dirs will be missing
#         build_temp = pathlib.Path(self.build_temp)
#         build_temp.mkdir(parents=True, exist_ok=True)
#         extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
#         extdir.mkdir(parents=True, exist_ok=True)

#         # example of cmake args
#         python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
#         python_include_path = str(pathlib.Path(
#             sys.executable).parent.parent.absolute()) + f'/include/python{python_version}/'

#         config = 'Debug' if self.debug else 'Release'
#         cmake_args = [
#             '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
#             '-DPY_INCLUDE_PATH=' + python_include_path,
#             '-DCMAKE_BUILD_TYPE=' + config,
#         ]

#         # example of build args
#         build_args = ['--config', config, '--', '-j4']

#         os.chdir(str(build_temp))
#         self.spawn(['cmake', str(cwd)] + cmake_args)
#         if not self.dry_run:
#             self.spawn(['cmake', '--build', '.'] + build_args)
#         # Troubleshooting: if fail on line above then delete all possible
#         # temporary CMake files including "CMakeCache.txt" in top level dir.
#         os.chdir(str(cwd))


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

setup_requires = ["cython"]
install_requires = []

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
]

setup(
    name="mDeepFRI",
    version=__version__,
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
