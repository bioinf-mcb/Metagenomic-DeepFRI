from setuptools import find_packages
from setuptools import setup

setup(
    name="meta_deepFRI",
    version="0.2.0",
    description="Pipeline for searching and aligning contact maps for proteins, then running DeepFri's GCN.",
    author="Piotr Kucharski",
    author_email="soliareofastorauj@gmail.com",
    url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    download_url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    license="GNU GPLv3",
    packages=find_packages(),
)
