# Installation
## Local

### Quick start
1. Upload structure files, for example from [PDB](https://www.rcsb.org/), to `STRUCTURE_FILES_PATH` (paths are defined in `CONFIG/FOLDER_STRUCTURE.py`)
2. Create target database
```{code-block} bash
python update_target_mmseqs_database.py --input all
```
3. Upload protein sequences `.faa` files into `QUERY_PATH`
4. Run main_pipeline.py.
```{code-block} bash
python main_pipeline.py --input all
```
5. Collect results from `FINISHED_PATH`

### Linux

1. Setup python environment
```{code-block} bash
pip install .
```
2. Install `mmseqs2` and `boost` libraries
```{code-block} bash
sudo apt install mmseqs2 libboost-numpy1.71 libboost-python1.71
```
3. (optional) Edit `CONFIG/FOLDER_STRUCTURE.py` to customize your folder structure
```{code-block} bash
nano CONFIG/FOLDER_STRUCTURE.py
```
4. Run `post_setup.py` script to create folder structure according to `FOLDER_STRUCTURE.py` and to download and unzip DeepFRI model weights
```{code-block} bash
python post_setup.py
```

### MacOS

1. Setup conda environment
```{code-block} bash
conda upgrade
```
and check what packages you need, for example:
```{code-block} bash
conda install pip
```
2. Setup python environment
```{code-block} bash
pip install .
```
3. Using Homebrew install `mmseqs2` and `boost` libraries
```{code-block} bash
brew install mmseqs2 boost-python3
```
4.  (optional) Edit `CONFIG/FOLDER_STRUCTURE.py` to customize your folder structure
```{code-block} bash
nano CONFIG/FOLDER_STRUCTURE.py
```
5. Run `post_setup.py` script to create folder structure according to `FOLDER_STRUCTURE.py` and to download and unzip DeepFRI model weights
```{code-block} bash
python post_setup.py
```

### Docker 

1. Create `YOUR_DATA_ROOT` directory on your local machine
```{code-block} bash
mkdir /YOUR_DATA_ROOT
```
2. Docker run! `-u $(id -u):$(id -g)` is used to ensure all files created by a pipeline are accessible to users
```{code-block} bash
docker run -it -u $(id -u):$(id -g) -v /YOUR_DATA_ROOT:/data soliareofastora/metagenomic-deepfri /bin/bash
```
3. Inside docker run `post_setup.py` script to create a folder structure and unzip DeepFRI model weights
```{code-block} bash
python3 post_setup.py
```
