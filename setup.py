from setuptools import find_packages
from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="meta_deepFRI",
    version="0.2.0",
    description="Pipeline for searching and aligning contact maps for proteins, then running DeepFri's GCN.",
    author="Piotr Kucharski",
    author_email="soliareofastorauj@gmail.com",
    url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    download_url="https://github.com/bioinf-mcb/Metagenomic-DeepFRI",
    license="GNU GPLv3",
    install_requires=required,
    packages=find_packages(),
)

print("If installing manually, please follow additional steps")
print("sudo apt-get install libboost-numpy1.71-dev libboost-python1.71-dev")
print("sudo apt install mmseqs2")
print("python post_setup.py")
