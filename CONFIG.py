import pathlib

MAX_CHAIN_LENGTH = 2500
ANGSTROM_CONTACT_THRESHOLD = 6

# user's folder structure
DATA_ROOT = pathlib.Path("/data")
STRUCTURE_FILES_PATH = DATA_ROOT / "structure_files"
QUERY_PATH = DATA_ROOT / "query"
FINISHED_PATH = DATA_ROOT / "finished"

# pipeline folders
SEQ_ATOMS_DATASET_PATH = DATA_ROOT / "seq_atoms_dataset"
MMSEQS_DATABASES_PATH = DATA_ROOT / "mmseqs_db"
WORK_PATH = DATA_ROOT / "workspace"

TARGET_DB_NAME = "targetDB"
DEEPFRI_TRAINED_MODELS_DOWNLOAD_URL = "https://users.flatironinstitute.org/vgligorijevic/public_www/DeepFRI_data/newest_trained_models.tar.gz"
DEEPFRI_MODEL_WEIGHTS_JSON_PATH = DATA_ROOT / "trained_models/model_config.json"
