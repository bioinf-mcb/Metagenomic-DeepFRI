# How it works

## Projects, tasks and timestamps
A pipeline is built around a folder structure described in `CONFIG / FOLDER_STRUCTURE.py`.

Multiple teams can have their separate **projects** - subdirectories inside `STRUCTURE_FILES_PATH`, `QUERY_PATH`, `WORK_PATH` and `FINISHED_PATH`.
You can execute `main_pipeline.py` or `update_target_mmseqs_database.py` with `--project_name` to easily control the files it touches.

**Without specifying `--project_name` pipeline will use `default` as a project name.**

The task is a single run of `main_pipeline.py`. The task name is the timestamp at the beginning of the script run. Its path is `WORK_PATH / project_name / timestamp`.
After completion, results will be stored in `FINISHED_PATH / project_name / timestamp`.

Target database creation `update_target_mmseqs_database.py` works in a similar fashion
appending new structures to `MMSEQS_DATABASES_PATH / project_name` creating a new timestamp folder.
The pipeline will use the database that was most recently created.

When running `main_pipeline.py` with a new project name,
the current state of `CONFIG / RUNTIME_PARAMETERS.py` will be saved in `WORK_PATH / project_name / project_config.json` and will be used in all upcoming tasks in this project.

Similarly `update_target_mmseqs_database.py`. It will store `MAX_TARGET_CHAIN_LENGTH` inside `target_db_config.json`.

## `mmseqs2` target database setup

1. Upload structure files to `STRUCTURE_FILES_PATH / your_project_name`.
2. Run `update_target_mmseqs_database.py` script.
   ```{code-block} bash
   python update_target_mmseqs_database.py --project_name your_project_name
   ```

The main feature of this project is its ability to generate a query contact map on the run
using results from `mmseqs2` target database search for similar protein sequences with known structures. Later in the `metagenomic_deepfri.py` contact map alignment is performed to use it as input to DeepFRI GCN. (implemented in CPP_lib/load_contact_maps.h)

`update_target_mmseqs_database.py` script will search for structure files,
process them and store protein chain sequence and atoms positions inside `SEQ_ATOMS_DATASET_PATH / project_name`.
It will also create a `mmseqs2` database in `MMSEQS_DATABASES_PATH / project_name`.
This operation will append new structures to existing ones.

Use `--input DIR_1 FILE_2 ...` argument list to parse structures from multiple sources. It accepts both absolute and relative `STRUCTURE_FILES_PATH` paths.
Use `--input .` to parse all structure files inside `STRUCTURE_FILES_PATH`.
Accepted formats are: `[".pdb", ".cif", ".ent"]` both raw and compressed `".gz"`

To add another structure file format edit `STRUCTURE_FILES_PARSERS` inside `update_target_mmseqs_database.py`

`target_db_config.json` contains `MAX_TARGET_CHAIN_LENGTH`.
This value is copied from `CONFIG / RUNTIME_PARAMETERS.py` while creating a new target database.

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