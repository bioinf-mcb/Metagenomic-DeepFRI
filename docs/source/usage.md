# How the pipeline works

The main feature of this project is the ability to generate a query contact map on the run
using results from `mmseqs2` target database search for similar protein sequences within known structures. In `metagenomic_deepfri.py` contact map alignment (implemented in CPP_lib/load_contact_maps.h) is performed to be later used as an input to DeepFRI GCN. 

A pipeline is built around following folder structure:

```{code-block} python
DATA_ROOT = data_root
# pipeline folders
STRUCTURE_FILES_PATH = data_root / "structure_files"
QUERY_PATH = data_root / "query"
WORK_PATH = data_root / "workspace"
FINISHED_PATH = data_root / "finished"

# database folders
SEQ_ATOMS_DATASET_PATH = data_root / "seq_atoms_dataset"
MMSEQS_DATABASES_PATH = data_root / "mmseqs_db"

# deepfri config file path
DEEPFRI_MODEL_WEIGHTS_JSON_FILE = data_root / "trained_models/model_config.json"
```

## Projects 
Project is a basic unit of workflow. It helps to avoid confusion from multiple pipeline runs.

Multiple users can have their separate **projects** - subdirectories inside `STRUCTURE_FILES_PATH`, `QUERY_PATH`, `WORK_PATH` and `FINISHED_PATH`. This allows to avoid set up of multiple databases locally, as they require a lot of physical memory.
You can execute `main_pipeline.py` or `update_target_mmseqs_database.py` with `--project_name` to easily control the files it works with.

**Without specifying `--project_name` pipeline will use `default` as a project name.**

## Tasks and timestamps

The task is a single run of `main_pipeline.py`. The task name is the timestamp at the beginning of the script run. Its path is `WORK_PATH / project_name / timestamp`. After completion, results will be stored in `FINISHED_PATH / project_name / timestamp`.

When running `main_pipeline.py` with a new project name, the current state of `CONFIG / RUNTIME_PARAMETERS.py` will be saved in `WORK_PATH / project_name / project_config.json` and will be used in all upcoming tasks in this project.

## Set up database

**Create a new database:**
1. Upload structure files to `STRUCTURE_FILES_PATH / your_project_name`.
2. Run `update_target_mmseqs_database.py` script.
   ```{code-block} bash
   python update_target_mmseqs_database.py --project_name your_project_name
   ```

The script will search for structure files, process them and store protein chain sequence and atoms positions inside `SEQ_ATOMS_DATASET_PATH / project_name`.
It will also create a `mmseqs2` database in `MMSEQS_DATABASES_PATH / project_name`. If the script is run for the second time it will append new structures to existing ones, creating a new folder named as a timestep. It would also store a `MAX_TARGET_CHAIN_LENGTH` parameter inside `target_db_config.json`.

Use `--input DIR_1 FILE_2 ...` argument list to parse structures from multiple sources. It accepts both absolute and relative `STRUCTURE_FILES_PATH` paths.
Use `--input .` to parse all structure files inside `STRUCTURE_FILES_PATH`.
Accepted formats are: `[".pdb", ".cif", ".ent"]` both raw and compressed `".gz"`

To add another structure file format edit `STRUCTURE_FILES_PARSERS` inside `update_target_mmseqs_database.py`

Protein ID is used as a filename. A new protein whose ID already exists in the database will be skipped.
Use `--overwrite` flag to overwrite existing sequences and atoms positions.
Use this argument if you want to apply changes to `MAX_TARGET_CHAIN_LENGTH` inside `target_db_config.json`


## Running the main pipeline

1. Upload `.faa` files into `QUERY_PATH / your_project_name` (default `project_name` is `default`)
2. Run main_pipeline.py
   ```{code-block} bash
   python main_pipeline.py --project_name your_project_name
   ```
3. Collect results from `FINISHED_PATH / your_project_name / timestamp`

The pipeline will attempt to use `project_name` target database name. If it's missing, the default target database will be used instead.
If you want to use another target database, use its name (project_name used during database creation) in `--target_db_name`.

Use `--input DIR_1 FILE_2 ...` argument list to process query `.faa` files from multiple sources.
Both absolute and relative to `QUERY_PATH`.
Use `--input .` to process all query `.faa` files inside `QUERY_PATH`.

`--delete_query` Use this flag so that source query files are deleted from input paths after being copied to the project workspace.

`--n_parallel_jobs` will divide query protein sequences evenly across all jobs.