import multiprocessing
CPU_COUNT = multiprocessing.cpu_count()

# max chain length used in update_mmseqs_database.py and in filtering false positive sequences from queries
MAX_CHAIN_LENGTH = 2500

# DeepFRI was trained with default value of 10 as seen in DeepFRI/train_DeepFRI.py
ANGSTROM_CONTACT_THRESHOLD = 6
# GENERATE_CONTACTS are used to fill gaps in target sequence during contact map alignment
GENERATE_CONTACTS = 2

# parameters used to filter mmseqs2 search results before aligning
MMSEQS_MIN_BIT_SCORE = -99999
MMSEQS_MAX_EVAL = 99999

# parameters used to calculate alignment score
PAIRWISE_ALIGNMENT_IDENTICAL = 2
PAIRWISE_ALIGNMENT_MISSMATCH = -1
PAIRWISE_ALIGNMENT_GAP_OPEN = -0.5
PAIRWISE_ALIGNMENT_GAP_EXTENDING = -0.1

# parameter to filter alignments based on sequence identity
ALIGNMENT_MIN_SEQUENCE_IDENTITY = 0
