import pathlib

#####################################
# EDIT THIS - user's folder structure
#####################################
# DATA_ROOT - as name suggests, root directory for all the pipeline files.
# STRUCTURE_FILES_PATH - path in which one should put their protein structure files. In any structure.
# QUERY_PATH - root for all query .faa files. They should be placed in subdirectories named after project_name
# WORK_PATH - contains all currently running or interrupted tasks. Each project folder also contain runtime parameters in project_config.json
# FINISHED_PATH - root for all finished pipeline results. Will contain folders named after project_name

DATA_ROOT = pathlib.Path("/data")   # You probably want to edit only this
STRUCTURE_FILES_PATH = DATA_ROOT / "structure_files"
QUERY_PATH = DATA_ROOT / "query"
WORK_PATH = DATA_ROOT / "workspace"
FINISHED_PATH = DATA_ROOT / "finished"


##########################################
# Shouldn't edit anything below this point
##########################################
# pipeline folders
SEQ_ATOMS_DATASET_PATH = DATA_ROOT / "seq_atoms_dataset"
MMSEQS_DATABASES_PATH = DATA_ROOT / "mmseqs_db"

# folders and files names
DEFAULT_NAME = "default"
TARGET_MMSEQS_DB_NAME = "targetDB"
SEQUENCES = "seq"
ATOMS = "atom"

PROJECT_CONFIG = "project_config.json"
TASK_CONFIG = "task_config.json"
JOB_CONFIG = "job_config.json"
TARGET_DB_CONFIG = "target_db_config.json"

ALIGNMENTS = "alignments.json"
MERGED_SEQUENCES = 'merged_sequences.faa'
MMSEQS_SEARCH_RESULTS = 'mmseqs2_search_results.m8'

# DeepFri model weights.
DEEPFRI_TRAINED_MODELS_DOWNLOAD_URL = "https://users.flatironinstitute.org/vgligorijevic/public_www/DeepFRI_data/newest_trained_models.tar.gz"
DEEPFRI_MODEL_WEIGHTS_JSON_FILE = DATA_ROOT / "trained_models/model_config.json"
