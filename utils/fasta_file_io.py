import dataclasses


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
