import pathlib

# user's folder structure
DATA_ROOT = pathlib.Path("/data")
STRUCTURE_FILES_PATH = DATA_ROOT / "structure_files"
QUERY_PATH = DATA_ROOT / "query"
FINISHED_PATH = DATA_ROOT / "finished"

# pipeline folders
WORK_PATH = DATA_ROOT / "workspace"
SEQ_ATOMS_DATASET_PATH = DATA_ROOT / "seq_atoms_dataset"
MMSEQS_DATABASES_PATH = DATA_ROOT / "mmseqs_db"

# names
DEFAULT_TARGET_DB_NAME = "default"
TARGET_MMSEQS_DB_NAME = "targetDB"
SEQUENCES = "seq"
ATOMS = "atom"

# DeepFri model weights. DO NOT CHANGE!
DEEPFRI_TRAINED_MODELS_DOWNLOAD_URL = "https://users.flatironinstitute.org/vgligorijevic/public_www/DeepFRI_data/newest_trained_models.tar.gz"
DEEPFRI_MODEL_WEIGHTS_JSON_FILE = DATA_ROOT / "trained_models/model_config.json"
