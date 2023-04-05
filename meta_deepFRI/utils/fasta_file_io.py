import hashlib
import dataclasses
import pathlib

from Bio import SeqIO

from meta_deepFRI.config.names import SEQUENCES, ATOMS


@dataclasses.dataclass
class SeqRecord:
    id: str
    seq: str


def load_fasta_file(file):
    """
    Loads FASTA file.

    Args:
        file (str): path to a FASTA file.

    Returns:
        List[SeqRecord]: list of records from FASTA file.
    """
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
    """Writes a FASTA file

    Args:
        seq_records (List[SeqRecords]): list of FASTA records.
        path (str): path to the output file.

    Returns:
        None
    """
    try:
        with open(path, "w", encoding="utf-8") as f:
            for seq_record in seq_records:
                f.write(f">{seq_record.id}\n{seq_record.seq}\n")
    except IsADirectoryError as e:
        raise IsADirectoryError(f"Path {path} is a directory. Please provide a path to a file.") from e
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Path {path} does not exist. Please provide a valid path.") from e


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


def hash_sequence_id(sequence: str):
    """Return the SHA256 encoding of protein sequence

    Args:
        sequence (str): Aminoacid sequence of the protein.

    Returns:
        str: SHA256 encoding of the protein sequence.
    """
    return hashlib.sha256(bytes(sequence, encoding='utf-8')).hexdigest()


def encode_faa_ids(path):
    seq_records = load_fasta_file(path)
    hash_lookup_dict = {}

    for i in range(len(seq_records)):
        seq_id_hash = hash_sequence_id(seq_records[i].id)
        while seq_id_hash in hash_lookup_dict.keys():
            seq_id_hash = hash_sequence_id(seq_records[i].id + "1")
        hash_lookup_dict[seq_id_hash] = seq_records[i].id
        seq_records[i].id = seq_id_hash

    faa_hashed_ids_path = str(path) + ".hashed_ids"
    write_fasta_file(seq_records, faa_hashed_ids_path)
    return faa_hashed_ids_path, hash_lookup_dict