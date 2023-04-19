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

1. Clone repo locally
```{code-block} bash
git clone https://github.com/bioinf-mcb/Metagenomic-DeepFRI --recursive
cd Metagenomic-DeepFRI
git submodule init
git submodule update --recursive --remote
```
2. Setup conda environment
```{code-block} bash
conda env create --name deepfri --file deepfri.yaml
conda activate deepfri
```
3. Install `meta-DeepFRI`
```{code-block} bash
pip install .
```
4. Verify installation
```{code-block} bash
pytest
deepfri --help
```

## Retrieve `DeepFRI` model weights

- [CPU weights](https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz)
- [GPU weights](https://users.flatironinstitute.org/~renfrew/DeepFRI_data/trained_models.tar.gz)

**Attention:** Using DeepFRI with GPU requires extra installation steps, described in section 4. GPU Setup of [tensorfow installation guide](https://www.tensorflow.org/install/pip).
## Database setup

1. Upload structure files to a folder.
2. Run `build_database.py` script.
   ```
   deepfri_build_db --input path/to/folder/with/strucures --output path/to/database
   ```
**Tip:** building a database from AF2 predicted structures took ~30 min.

Use parameter `-max_len` to define maximal length of the protein. Due to initial DeepFRI training set
default value is set to `1000 aa`.

Main feature of this project is its ability to generate query contact map on the fly
using results from mmseqs2 target database search for similar protein sequences with known structures.
Later in the `metagenomic_deepfri.py` contact map alignment is performed to use it as input to DeepFRI GCN.
(implemented in CPP_lib/load_contact_maps.h)

`build_database.py` script will search for structure files,
process them and store protein sequence and atoms positions inside `database/SEQ_ATOMS_DATASET_PATH`.
It will also create a mmseqs2 database in `database/MMSEQS_DATABASES_PATH`.

You can also use `--input DIR_1 FILE_2 ...` argument list to parse structures from multiple sources.
Accepted formats are: `.pdb .cif .ent` both raw and compressed with `.gz`

Protein ID is used as a filename. A new protein whose ID already exists in the database will be skipped.
Use `--overwrite` flag to overwrite existing sequences and atoms positions.

## Running main pipeline

1. Run main_pipeline.py
   ```
   deepfri -i /path/to/protein/sequences -db /path/to/database/folder/from/previous/step -c /path/to/deepfri/confg -o /output_path
   ```
**Attention:** Single instance of DeepFRI on GPU requires 10GB VRAM.
Other available parameters can be found upon command `deepfri --help`.

## Results
Finished folder will contain:
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

### Contact

Piotr Kucharski - soliareofastorauj@gmail.com \
Valentyn Bezshapkin - valentyn.bezshapkin@micro.biol.ethz.ch
