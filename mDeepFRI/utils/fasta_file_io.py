import pathlib

from pysam.libcfaidx import FastxFile

from mDeepFRI import ATOMS, SEQUENCES


class SeqFileLoader:
    def __init__(self, path):
        self.path = pathlib.Path(path)

    def __getitem__(self, target_id):
        protein_file = (self.path / SEQUENCES / (target_id + ".faa"))
        structure_file = (self.path / ATOMS / (target_id + ".bin"))
        assert protein_file.exists(), f"Target database {self.path.name} does not contain sequence" \
                                      f"{self.path / SEQUENCES / (target_id + '.faa')}. " \
                                      "It should not happen. Please create new target database."
        assert structure_file.exists(), f"Target database {self.path.name} does not contain structure" \
                                        f"{self.path / ATOMS / (target_id + '.bin')}. " \
                                        "It should not happen. Please create new target database."

        with FastxFile(self.path / SEQUENCES / (target_id + ".faa"), "r") as f:
            sequence = [record.sequence for record in f][0]

        return sequence
