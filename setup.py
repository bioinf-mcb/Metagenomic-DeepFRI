import os
import pathlib
import sys

from setuptools import find_namespace_packages, setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig

# TODO: Add __init__.py file to meta_deepFRI/DeepFRI/


class CMakeExtension(Extension):

    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        # example of cmake args
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        python_include_path = str(pathlib.Path(
            sys.executable).parent.parent.absolute()) + f'/include/python{python_version}/'

        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DPY_INCLUDE_PATH=' + python_include_path,
            '-DCMAKE_BUILD_TYPE=' + config,
        ]

        # example of build args
        build_args = ['--config', config, '--', '-j4']

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        # Troubleshooting: if fail on line above then delete all possible
        # temporary CMake files including "CMakeCache.txt" in top level dir.
        os.chdir(str(cwd))


setup(
    name="meta_deepFRI",
    version="0.2.0",
    description="Pipeline for searching and aligning contact maps for proteins, then running DeepFri's GCN.",
    author="Piotr Kucharski, Valentyn Bezshapkin",
    author_email="soliareofastorauj@gmail.com, valentyn.bezshapkin@micro.biol.ethz.ch",
    url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    download_url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "deepfri = scripts.main_pipeline:main",
            "deepfri_db_build = scripts.build_database:main",
        ],
    },
    ext_modules=[CMakeExtension("meta_deepFRI/CPP_lib/library_definition")],
    cmdclass={
        'build_ext': build_ext,
    },
    license="GNU GPLv3",
    packages=find_namespace_packages(),
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
