import pathlib

from Bio import SeqIO

from CONFIG.FOLDER_STRUCTURE import SEQUENCES, ATOMS


class SeqFileLoader:
    def __init__(self, path):
        self.path = pathlib.Path(path)

    def __getitem__(self, target_id):
        assert (self.path / SEQUENCES / (target_id + ".faa")).exists(), f"Target database {self.path.name} contains ID for SEQUENCE.faa that is not in the {self.path / SEQUENCES / (target_id + '.faa')}. " \
                                                                        f"It should not happen. Please create new target database."
        assert (self.path / ATOMS / (target_id + ".bin")).exists(), f"Target database {self.path.name} contains ID to ATOM.bin that is not in the {self.path / ATOMS / (target_id + '.bin')}. " \
                                                                    f"It should not happen. Please create new target database."

        with open(self.path / SEQUENCES / (target_id + ".faa"), "r") as f:
            sequence = SeqIO.read(f, "fasta").seq
        return sequence
