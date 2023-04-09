from Bio import SeqUtils

PROTEIN_LETTERS = dict()
for k, v in SeqUtils.IUPACData.protein_letters_3to1_extended.items():
    PROTEIN_LETTERS[str.upper(k)] = v
PROTEIN_LETTERS["UNK"] = "X"
