import pathlib

from Bio import SeqUtils


MAX_CHAIN_LENGTH = 2500

PROTEIN_LETTERS = dict()
for k, v in SeqUtils.IUPACData.protein_letters_3to1_extended.items():
    PROTEIN_LETTERS[str.upper(k)] = v
PROTEIN_LETTERS["UNK"] = "X"


STRUCTURE_FILES_PATTERNS = [
    '.pdb',
    '.pdb.gz',
    '.cif',
    '.cif.gz'
]

# folder structure
DATA_PATH = pathlib.Path.home() / "data"
STRUCTURE_FILES_PATH = DATA_PATH / "structure_files"
CONTACT_MAP_DATASET_PATH = DATA_PATH / "contact_map_dataset"
TARGET_MMSEQS_DATABASE_PATH = DATA_PATH / "mmseqs_db"
DEFAULT_MMSEQS_NAME = DATA_PATH / "targetDB"
TMP_FOLDER_PATH = DATA_PATH / "TMP"




