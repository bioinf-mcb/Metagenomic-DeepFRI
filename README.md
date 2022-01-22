# Metagenomic-DeepFRI

## About The Project
Pipeline for searching and aligning contact maps for proteins, then running DeepFri's GCN.
### Built With

* [DeepFRI](https://github.com/flatironinstitute/DeepFRI)
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [Boost.Python](https://www.boost.org/doc/libs/1_75_0/libs/python/doc/html/index.html)

## Installation
### Local
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
5. Run `post_setup.py` script to create folder structure according to `CONFIG.py` and to download DeepFRI model weights
   ```
   python post_setup.py
   ```
### Docker
1. Create `DATA_ROOT` directory on your local machine
   ```
   mkdir /DATA_ROOT
   ```
2. Docker run! `-u $(id -u):$(id -g)` is used to make sure all files created by pipeline are accessible for users
   ```
   docker run -it -u $(id -u):$(id -g) -v /DATA_ROOT:/data soliareofastora/metagenomic-deepfri
   ```
3. Inside docker run `post_setup.py` script to create folder structure and to unzip DeepFRI model weights
   ```
   python post_setup.py
   ```

## Mmseqs2 target database setup 

Main feature of this project  is its ability to find similar protein chains 
using mmseqs2 with known structures to use aligned contact maps as input to GCN from DeepFRI.

`update_target_mmseqs_database.py` script will parse structure files and store protein chain sequence and atoms positions inside `SEQ_ATOMS_DATASET_PATH`.
It will also create a mmseqs2 database in `MMSEQS_DATABASES_PATH`. This operation will append new structures to existing ones.

Protein ID is used as a filename. A new protein whose ID already exists in the database will be skipped.
Use `--overwrite` flag to overwrite existing sequences and atoms positions.

1. Upload structure files to `STRUCTURE_FILES_PATH`. Accepted formats are: `.pdb .cif .ent` both raw and compressed `.gz`
2. Run `update_target_mmseqs_database.py` script. 
   ```
   python update_target_mmseqs_database.py
   ```
You can also use `--input PATH_1 PATH_2` argument to parse structures from multiple directories.

Use argument `--output DB_NAME` to specify database name that you can later use in `main.py --target_db DB_NAME`.
By default, pipeline will use `default` database name.

To add another structure file format edit `STRUCTURE_FILES_PARSERS` inside `update_target_mmseqs_database.py :26`

## Running DeepFri predictions

1. Upload `.faa` files into `QUERY_PATH/task_name` (default `task_name` is `default`)
2. Run `python3 main.py -task_name`
4. Collect results from `FINISHED_PATH/task_name`

## Useful docker commands
```
docker build -t soliareofastora/metagenomic-deepfri .
docker push soliareofastora/metagenomic-deepfri
docker pull soliareofastora/metagenomic-deepfri
docker run -it -u $(id -u):$(id -g) -v /data:/data soliareofastora/metagenomic-deepfri
```
## Contributing

If you have a suggestion that would make this project better, email me or fork the repo and create a pull request.

### Contact

Piotr Kucharski - soliareofastorauj@gmail.com
