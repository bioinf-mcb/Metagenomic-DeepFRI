# Metagenomic-DeepFRI

## About The Project
Do you have **thousands of protein sequences** with **unknown structures**, but still want to know their
molecular function, biological process, cellular component and enzyme commission **predicted by DeepFRI Graph Convolutional Network?**

This is the right project for this task! Pipeline in a nutshell:
1. Search for similar target protein sequences using MMseqs2
2. Align target protein contact map to fit your query protein with unknown structure
3. Run predictions on query sequence combined with aligned target contact map or sequence alone if no alignment was found

### Built With

* [DeepFRI](https://github.com/SoliareofAstora/DeepFRI)
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [Boost.Python](https://www.boost.org/doc/libs/1_75_0/libs/python/doc/html/index.html)

# Installation
## Local
### Linux
1. Setup python environment
    ```
    pip install .
    ```
2. Install mmseqs2
    ```
    sudo apt install mmseqs2
   ```
3. Install boost libraries
    ```
    sudo apt-get install libboost-numpy1.71 libboost-python1.71
   ```
4. (optional) Edit `CONFIG/FOLDER_STRUCTURE.py` to customize your folder structure
   ```
   nano CONFIG/FOLDER_STRUCTURE.py
   ```
5. Run `post_setup.py` script to create folder structure according to `FOLDER_STRUCTURE.py` and to download and unzip DeepFRI model weights
   ```
   python post_setup.py
   ```
### MacOS
1. Setup conda environment
    ```
    conda upgrade
    ```
   and check what packages you need, for example:
    ```
    conda install pip
    ```
2. Setup python environment
    ```
    pip install .
    ```
3. Using Homebrew ($(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"), install mmseqs2:
   ```
    brew install mmseqs2
    ```
4. Using Homebrew install boost-python3:
   ```
    brew install boost-python3
    ```
5.  (optional) Edit `CONFIG/FOLDER_STRUCTURE.py` to customize your folder structure
   ```
   nano CONFIG/FOLDER_STRUCTURE.py
   ```
6. Run `post_setup.py` script to create folder structure according to `FOLDER_STRUCTURE.py` and to download and unzip DeepFRI model weights
   ```
   python post_setup.py
   ```
### Docker (covers steps 1-3 from Linux and 1-4 from MacOS Local installation)
1. Create `YOUR_DATA_ROOT` directory on your local machine
   ```
   mkdir /YOUR_DATA_ROOT
   ```
2. Docker run! `-u $(id -u):$(id -g)` is used to make sure all files created by pipeline are accessible for users
   ```
   docker run -it -u $(id -u):$(id -g) -v /YOUR_DATA_ROOT:/data soliareofastora/metagenomic-deepfri /bin/bash
   ```
3. Inside docker run `post_setup.py` script to create folder structure and unzip DeepFRI model weights
   ```
   python3 post_setup.py
   ```
### TL:DR! QUICK START
1. Upload structure files, for example from [PDB](https://www.rcsb.org/), to `STRUCTURE_FILES_PATH` (paths are defined in `CONFIG/FOLDER_STRUCTURE.py`)
3. Create target database
   ```
   python update_target_mmseqs_database.py --input all
   ```
4. Upload protein sequences `.faa` files into `QUERY_PATH`
5. Run main_pipeline.py.
   ```
   python main_pipeline.py --input all
   ```
6. Collect results from `FINISHED_PATH`

# How this pipeline works
## Projects, tasks and timestamps
Pipeline is build around folder structure described in `CONFIG / FOLDER_STRUCTURE.py`.

Multiple teams can have their separate **projects** - subdirectories inside `STRUCTURE_FILES_PATH`, `QUERY_PATH`,
`WORK_PATH` and `FINISHED_PATH`.
You can execute `main_pipeline.py` or `update_target_mmseqs_database.py` with `--project_name` to easily control the files it touches.

**Without specifying `--project_name` pipeline will use `default` as project name.**

Task is a single run of `main_pipeline.py`. The task name is the timestamp at the beginning of the scrip run. Its path is `WORK_PATH / project_name / timestamp`.
After completion, results will be stored in `FINISHED_PATH / project_name / timestamp`.

Target database creation `update_target_mmseqs_database.py` works in similar fashion
appending new structures to `MMSEQS_DATABASES_PATH / project_name` creating new timestamp folder.
Pipeline will use the database that was most recently created.

When running `main_pipeline.py` with a new project name,
current state of `CONFIG / RUNTIME_PARAMETERS.py` will be saved in `WORK_PATH / project_name / project_config.json`
and will be used in all upcoming tasks in this project.

Similarly `update_target_mmseqs_database.py`. It will store `MAX_TARGET_CHAIN_LENGTH` inside `target_db_config.json`.

## Mmseqs2 target database setup

1. Upload structure files to `STRUCTURE_FILES_PATH / your_project_name`.
2. Run `update_target_mmseqs_database.py` script.
   ```
   python update_target_mmseqs_database.py --project_name your_project_name
   ```

Main feature of this project is its ability to generate query contact map on the fly
using results from mmseqs2 target database search for similar protein sequences with known structures.
Later in the `metagenomic_deepfri.py` contact map alignment is performed to use it as input to DeepFRI GCN.
(implemented in CPP_lib/load_contact_maps.h)

`update_target_mmseqs_database.py` script will search for structure files,
process them and store protein chain sequence and atoms positions inside `SEQ_ATOMS_DATASET_PATH / project_name`.
It will also create a mmseqs2 database in `MMSEQS_DATABASES_PATH / project_name`.
This operation will append new structures to existing ones.

You can also use `--input DIR_1 FILE_2 ...` argument list to parse structures from multiple sources.
Both absolute and relative to `STRUCTURE_FILES_PATH`.
Use `--input .` to parse all structure files inside `STRUCTURE_FILES_PATH`.
Accepted formats are: `.pdb .cif .ent` both raw and compressed `.gz`

To add another structure file format edit `STRUCTURE_FILES_PARSERS` inside `update_target_mmseqs_database.py`

`target_db_config.json` contains `MAX_TARGET_CHAIN_LENGTH`.
This value is copied from `CONFIG / RUNTIME_PARAMETERS.py` while creating new target database.

Protein ID is used as a filename. A new protein whose ID already exists in the database will be skipped.
Use `--overwrite` flag to overwrite existing sequences and atoms positions.
Also use this argument if you want to apply changes to `MAX_TARGET_CHAIN_LENGTH` inside `target_db_config.json`

## Running main pipeline

1. Upload `.faa` files into `QUERY_PATH / your_project_name` (default `project_name` is `default`)
2. Run main_pipeline.py
   ```
   python main_pipeline.py --project_name your_project_name
   ```
3. Upon completion, collect results from `FINISHED_PATH / your_project_name / timestamp`

Pipeline will attempt to use `project_name` target database name. If it's missing, default target database will be used instead.
If you want to use other target database use its name (project_name used during database creation) in `--target_db_name`.

You can use `--input DIR_1 FILE_2 ...` argument list to process query `.faa` files from multiple sources.
Both absolute and relative to `QUERY_PATH`.
Use `--input .` to process all query `.faa` files inside `QUERY_PATH`.

`--delete_query` Use this flag so that source query files are deleted from input paths after being copied to project workspace.

`--n_parallel_jobs` will divide query protein sequences evenly across all jobs.

## Results
Finished folder `FINISHED_PATH / project_name / timestamp` will contain:
1. `query_files/*` - directory containing all input query files.
2. `mmseqs2_search_results.m8`
3. `alignments.json` - results of alignment search implemented in `utils.search_alignments.py`
4. `metadata*` - files with some useful info
5. `results*` - multiple files from DeepFRI. Organized by model type ['GCN' / 'CNN'] and its mode ['mf', 'bp', 'cc', 'ec'] for the total of 8 files.
Sometimes results from one model can be missing which means that all query proteins sequences were aligned correctly or none of them were aligned.
   ```
   mf = molecular_function
   bp = biological_process
   cc = cellular_component
   ec = enzyme_commission
   ```

## Contributing

If you have a suggestion that would make this project better, email me or fork the repo and create a pull request.

### TODO
1. `main_pipeline.py` add possibility to use specific target_database path and timestamp instead of name only
2. `utils/search_alignments.py` make some runtime tests, maybe chunkified sequences will perform better with pathos.multiprocessing
3. `update_target_mmseqs_database.py` add max_target_chain_length argument and inform user if there is difference between this arg and existing target_db_config.json
4. `update_target_mmseqs_database.py` when already processed structures to another project, check if they already exists somewhere

### Contact

Piotr Kucharski - soliareofastorauj@gmail.com
