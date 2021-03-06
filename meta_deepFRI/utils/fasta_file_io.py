import dataclasses
import pathlib

from Bio import SeqIO

from meta_deepFRI.config.names import SEQUENCES, ATOMS


@dataclasses.dataclass
class SeqRecord:
    id: str
    seq: str


def load_fasta_file(file):
    seq_records = []
    with open(file, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].replace("\n", "")
                seq_records.append(SeqRecord(id=seq_id, seq=""))
            else:
                seq_records[-1].seq += line.replace("\n", "")
    return seq_records


def write_fasta_file(seq_records, path):
    with open(path, "w") as f:
        for seq_record in seq_records:
            f.write(f">{seq_record.id}\n{seq_record.seq}\n")


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
