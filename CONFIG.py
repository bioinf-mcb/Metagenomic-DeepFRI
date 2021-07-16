import pathlib

from Bio import SeqUtils


MAX_CHAIN_LENGTH = 2500
ANGSTROM_CONTACT_THRESHOLD = 6

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
DATA_ROOT = pathlib.Path.home() / "data"
STRUCTURE_FILES_PATH = DATA_ROOT / "structure_files"
ATOMS_DATASET_PATH = DATA_ROOT / "atoms_dataset"
MMSEQS_DATABASES_PATH = DATA_ROOT / "mmseqs_db"
DEFAULT_MMSEQS_NAME = "targetDB"

WORK_PATH = DATA_ROOT / "workspace"
QUERY_FOLDER_PATH = DATA_ROOT / "query"
TMP_FOLDER_PATH = DATA_ROOT / "TMP"




