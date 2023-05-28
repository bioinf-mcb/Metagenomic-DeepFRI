import numpy as np
import parasail
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


# parasail martices
substitution_matrices = {
    "blosum30": parasail.blosum30,
    "blosum35": parasail.blosum35,
    "blosum40": parasail.blosum40,
    "blosum45": parasail.blosum45,
    "blosum50": parasail.blosum50,
    "blosum55": parasail.blosum55,
    "blosum60": parasail.blosum60,
    "blosum62": parasail.blosum62,
    "blosum65": parasail.blosum65,
    "blosum70": parasail.blosum70,
    "blosum75": parasail.blosum75,
    "blosum80": parasail.blosum80,
    "blosum85": parasail.blosum85,
    "blosum90": parasail.blosum90,
    "blosum100": parasail.blosum100,
    "pam10": parasail.pam10,
    "pam20": parasail.pam20,
    "pam30": parasail.pam30,
    "pam40": parasail.pam40,
    "pam50": parasail.pam50,
    "pam60": parasail.pam60,
    "pam70": parasail.pam70,
    "pam80": parasail.pam80,
    "pam90": parasail.pam90,
    "pam100": parasail.pam100,
    "pam110": parasail.pam110,
    "pam120": parasail.pam120,
    "pam130": parasail.pam130,
    "pam140": parasail.pam140,
    "pam150": parasail.pam150,
    "pam160": parasail.pam160,
    "pam170": parasail.pam170,
    "pam180": parasail.pam180,
    "pam190": parasail.pam190,
    "pam200": parasail.pam200,
    "pam210": parasail.pam210,
    "pam220": parasail.pam220,
    "pam230": parasail.pam230,
    "pam240": parasail.pam240,
    "pam250": parasail.pam250,
    "pam260": parasail.pam260,
    "pam270": parasail.pam270,
    "pam280": parasail.pam280,
    "pam290": parasail.pam290,
    "pam300": parasail.pam300,
    "pam310": parasail.pam310,
    "pam320": parasail.pam320,
    "pam330": parasail.pam330,
    "pam340": parasail.pam340,
    "pam350": parasail.pam350,
    "pam360": parasail.pam360,
    "pam370": parasail.pam370,
    "pam380": parasail.pam380,
    "pam390": parasail.pam390,
    "pam400": parasail.pam400,
    "pam410": parasail.pam410,
    "pam420": parasail.pam420,
    "pam430": parasail.pam430,
    "pam440": parasail.pam440,
    "pam450": parasail.pam450,
    "pam460": parasail.pam460,
    "pam470": parasail.pam470,
    "pam480": parasail.pam480,
    "pam490": parasail.pam490,
    "pam500": parasail.pam500
}
