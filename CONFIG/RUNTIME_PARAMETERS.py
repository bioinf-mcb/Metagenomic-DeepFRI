import multiprocessing
CPU_COUNT = multiprocessing.cpu_count()

# MAX_TARGET_CHAIN_LENGTH is used in update_target_mmseqs_database.py:67 to truncate longer sequences
MAX_TARGET_CHAIN_LENGTH = 2500
# MAX_QUERY_CHAIN_LENGTH is used to filter out query sequences that are too long
MAX_QUERY_CHAIN_LENGTH = 2500

# DeepFRI was trained with default value of 10 as seen in DeepFRI/train_DeepFRI.py
ANGSTROM_CONTACT_THRESHOLD = 6
# GENERATE_CONTACTS are used to fill gaps in target sequence during contact map alignment
GENERATE_CONTACTS = 2

# parameters used to filter mmseqs2 search results before aligning
MMSEQS_MIN_BIT_SCORE = -99999
MMSEQS_MAX_EVAL = 99999

# parameters used to calculate alignment score
PAIRWISE_ALIGNMENT_MATCH = 2
PAIRWISE_ALIGNMENT_MISSMATCH = -1
PAIRWISE_ALIGNMENT_GAP_OPEN = -0.5
PAIRWISE_ALIGNMENT_GAP_CONTINUATION = -0.1

# parameter to filter alignments based on sequence identity
# ALIGNMENT_MIN_SEQUENCE_IDENTITY value should be between 0 and 1
ALIGNMENT_MIN_SEQUENCE_IDENTITY = 0.3

# this list determines types of outputs of metagenomic_deepfri
# mf = molecular_function
# bp = biological_process
# cc = cellular_component
# ec = enzyme_commission
# ['mf', 'bp', 'cc', 'ec']
DEEPFRI_PROCESSING_MODES = ['mf', 'bp', 'cc', 'ec']
