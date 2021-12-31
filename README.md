# Metagenomic-DeepFRI

## About The Project

### Built With

* [DeepFRI](https://github.com/flatironinstitute/DeepFRI)
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [Docker](https://www.docker.com/)

## Installation
### Docker
1. Create `DATA_ROOT` directory locally on your local machine
   ```
   mkdir DATA_ROOT
   ```
2. Run docker! `-u $(id -u):$(id -g)` is used to make sure all files created by container are accessible for users
   ```
   docker run -it -u $(id -u):$(id -g) -v /DATA_ROOT:/data soliareofastora/metagenomic-deepfri
   ```
3. Inside docker run `post_setup.py` script to create folder structure and to download DeepFRI model weights
   ```
   python3 post_setup.py
   ```
 
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
4. Edit `CONFIG.py` to set up your folder structure
   ```
   nano CONFIG.py
   ```
5. Run `post_setup.py` script to create folder structure and to download DeepFRI model weights
   ```
   python post_setup.py
   ```

## Mmseqs2 database setup 

Main feature of this project  is its ability to find similar protein chains 
using mmseqs2 with known structures to use aligned contact maps as input to GCN from DeepFRI.

1. Upload structure files to `STRUCTURE_FILES_PATH`. Accepted formats are: .pdb .pdb.gz .cif .cif.gz
2. Run `update_mmaseqs_database.py` script. You can also use `-i YOUR_PATH` argument to omit first step.

This script will parse structure files and store protein chain sequence and atoms positions inside `SEQ_ATOMS_DATASET_PATH`.
It will also create a mmseqs2 database in `MMSEQS_DATABASES_PATH`. This operation will append new structures to existing ones.


Protein ID is used as a filename. A new protein whose ID already exists in the database will be skipped.
Use `--overwrite` flag to overwrite existing sequences and atoms positions.

## Running experiments

1. Upload `.faa` files into `QUERY_PATH`
2. Run `python3 main_pipeline.py`
4. Collect results from `FINISHED_PATH`

## Useful docker commands
```
docker build -t soliareofastora/metagenomic-deepfri .
docker push soliareofastora/metagenomic-deepfri
docker pull soliareofastora/metagenomic-deepfri
docker run -it -u $(id -u):$(id -g) -v /data:/data soliareofastora/metagenomic-deepfri
```
## Contributing

If you have a suggestion that would make this project better, email me or fork the repo and create a pull request.

## Contact

Piotr Kucharski - soliareofastorauj@gmail.com
