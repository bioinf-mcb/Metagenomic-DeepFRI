import numpy as np
from pysam import FastaFile, FastxFile


def load_fasta_as_dict(fasta_file):
    """Load FASTA file as dict"""
    with FastxFile(fasta_file) as fasta:
        fasta_dict = {entry.name: entry.sequence for entry in fasta}
    return fasta_dict


def retrieve_fasta_entries_as_dict(fasta_file, entries):
    """Retrieve selected FASTA entries as dict"""
    fasta_dict = dict()
    with FastaFile(fasta_file) as fasta:
        for name in entries:
            fasta_dict[name] = fasta.fetch(name)
    return fasta_dict


# copied from Biopython to remove dependency
protein_letters_1to3 = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}
protein_letters_1to3_extended = {
    **protein_letters_1to3,
    **{
        "B": "Asx",
        "X": "Xaa",
        "Z": "Glx",
        "J": "Xle",
        "U": "Sec",
        "O": "Pyl"
    },
}

protein_letters_3to1_extended = {
    value: key
    for key, value in protein_letters_1to3_extended.items()
}

PROTEIN_LETTERS = dict()
for k, v in protein_letters_3to1_extended.items():
    PROTEIN_LETTERS[str.upper(k)] = v
PROTEIN_LETTERS["UNK"] = "X"


def seq2onehot(seq):
    """Create 26-dim embedding"""
    chars = [
        '-', 'D', 'G', 'U', 'L', 'N', 'T', 'K', 'H', 'Y', 'W', 'C', 'P', 'V',
        'S', 'O', 'I', 'E', 'F', 'X', 'Q', 'A', 'B', 'Z', 'R', 'M'
    ]
    vocab_size = len(chars)
    vocab_embed = dict(zip(chars, range(vocab_size)))

    # Convert vocab to one-hot
    vocab_one_hot = np.zeros((vocab_size, vocab_size), int)
    for _, val in vocab_embed.items():
        vocab_one_hot[val, val] = 1

    embed_x = [vocab_embed[v] for v in seq]
    seqs_x = np.array([vocab_one_hot[j, :] for j in embed_x])

    return seqs_x
