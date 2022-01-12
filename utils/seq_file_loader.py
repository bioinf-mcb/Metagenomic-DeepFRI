import pathlib
from Bio import SeqIO
from CONFIG.FOLDER_STRUCTURE import SEQUENCES


class SeqFileLoader:
    def __init__(self, path):
        self.path = pathlib.Path(path)

    def __getitem__(self, protein_id):
        with open(self.path / SEQUENCES / (protein_id + ".faa"), "r") as f:
            sequence = SeqIO.parse(f, "fasta").seq
        return sequence
