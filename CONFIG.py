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

FOLDSEEK_DATABASES_PATH = DATA_ROOT / "foldseek_db"
MMSEQS_DATABASES_PATH = DATA_ROOT / "mmseqs_db"

TARGET_DB_NAME = "targetDB"

WORK_PATH = DATA_ROOT / "workspace"
QUERY_PATH = DATA_ROOT / "query"
TMP_PATH = DATA_ROOT / "TMP"
FINISHED_PATH = DATA_ROOT / "finished"


FOLDSEEK_BIN_PATH = "/home/soliareofastora/foldseek/bin/"

