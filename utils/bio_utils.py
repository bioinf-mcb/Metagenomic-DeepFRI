from Bio import SeqUtils

from utils.structure_files_parsers.parse_pdb import parse_pdb
from utils.structure_files_parsers.parse_mmcif import parse_mmcif

PROTEIN_LETTERS = dict()
for k, v in SeqUtils.IUPACData.protein_letters_3to1_extended.items():
    PROTEIN_LETTERS[str.upper(k)] = v
PROTEIN_LETTERS["UNK"] = "X"

STRUCTURE_FILES_PARSERS = {
    '.pdb': parse_pdb,
    '.pdb.gz': parse_pdb,
    '.cif': parse_mmcif,
    '.cif.gz': parse_mmcif,
    '.ent': parse_pdb,
    '.ent.gz': parse_pdb
}
