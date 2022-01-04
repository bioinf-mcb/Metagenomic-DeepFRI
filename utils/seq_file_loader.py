import pathlib

from Bio import SeqIO


class SeqFileLoader:
    def __init__(self, path):
        self.path = pathlib.Path(path)

    def __getitem__(self, protein_id):
        with open(self.path / "seq" / protein_id, "r") as f:
            sequence = SeqIO.parse(f, "fasta").seq
        return sequence
