# Metagenomic-DeepFRI

## About The Project
Do you have **thousands of protein sequences** with **unknown structures**, but still want to know their
molecular function, biological process, cellular component and enzyme commission **predicted by DeepFRI Graph Convolutional Network?**

This is the right project for this task! Pipeline in a nutshell:
1. Search for similar target protein sequences using MMseqs2.
2. Align target protein contact map to fit your query protein with unknown structure.
3. Run predictions on query sequence combined with aligned target contact map or sequence alone if no alignment was found.

### Built With

* [DeepFRI](https://github.com/SoliareofAstora/DeepFRI)
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [Boost.Python](https://www.boost.org/doc/libs/1_75_0/libs/python/doc/html/index.html)
* [Parasail](https://github.com/jeffdaily/parasail)

# Installation

## 1. Install environment and DeepFRI

1. Clone repo locally
```{code-block} bash
git clone https://github.com/bioinf-mcb/Metagenomic-DeepFRI
cd Metagenomic-DeepFRI
```
2. Setup conda environment
```{code-block} bash
conda env create --name deepfri --file environment.yml
conda activate deepfri
```
3. Install `mDeepFRI`
```{code-block} bash
pip install .
```
4. Verify installation
```{code-block} bash
pytest
mDeepFRI --help
```

# Usage
## 1. Download models
Run command:
```
mDeepFRI get-models --output path/to/weights/folder
```

## 2. Prepare database

1. Upload structure (`.pdb` or `.mmcif`) files to a folder in your system.
2. Run command:
```
mDeepFRI build-db --input path/to/folder/with/strucures --output path/to/database -t threads
```
**Tip:** building a database from AF2Swissprot (~550k predicted structures) on 32 CPU cores took ~30 min.

Use parameter `-max_len` to define maximal length of the protein. Due to initial DeepFRI training set default value is set to `1000`.

Main feature of this project is its ability to generate query contact map on the fly
using results from mmseqs2 target database search for similar protein sequences with known structures.
Later in the `metagenomic_deepfri.py` contact map alignment is performed to use it as input to DeepFRI GCN.
(implemented in CPP_lib/load_contact_maps.h)

The command will search for structure files,
process them and store protein sequence and atoms positions inside `database/seq_atom_db`.
It will also create a mmseqs2 database within `database/`.

You can also use `--input DIR_1 FILE_2 ...` argument list to parse structures from multiple sources.
Accepted formats are: `.pdb`, `.cif`, `.ent` both raw and compressed with `.gz`

Protein ID is used as a filename. A new protein whose ID already exists in the database will be skipped.
Use `--overwrite` flag to overwrite existing sequences and atoms positions.

## 3. Predict protein function
```
mDeepFRI predict-function -i /path/to/protein/sequences -d /path/to/database/folder/from/previous/step -w /path/to/deepfri/weights/folder -o /output_path
```
**Attention:** Single instance of DeepFRI on GPU requires 10GB VRAM.

Other available parameters can be found upon command `mDeepFRI --help`.

## Results
Finished folder will contain:
1. `query_files/*` - directory containing all input query files.
2. `mmseqs2_search_results.m8`
3. `alignments.json` - results of alignment search implemented in `utils.search_alignments.py`
4. `metadata*` - files with some useful info
5. `results*` - multiple files from DeepFRI. Organized by model type (`GCN` or `CNN`) and its mode (`mf`, `bp`, `cc`, `ec`) for the total of 8 files.
Sometimes results from one model can be missing which means that all query proteins sequences were aligned correctly or none of them were aligned.
   ```
   mf = molecular_function
   bp = biological_process
   cc = cellular_component
   ec = enzyme_commission
   ```

## GPU / CPU utilization
If CUDA is installed on your machine, `metaDeepFRI` will automatically use it for prediction, no additional installations are needed. If not, the model will use CPUs. If argument `threads` is provided, the prediction will run on multiple CPU cores.

## Citations
If you use this software please cite:
- Gligorijević et al. "Structure-based protein function prediction using graph convolutional networks" Nat. Comms. (2021). https://doi.org/10.1038/s41467-021-23303-9
- Steinegger & Söding "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets" Nat. Biotechnol. (2017) https://doi.org/10.1038/nbt.3988
- Maranga et al. "Comprehensive Functional Annotation of Metagenomes and Microbial Genomes Using a Deep Learning-Based Method" mSystems (2023) https://doi.org/10.1128/msystems.01178-22
- Daily "Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments" BMC Bioinform. (2016) https://doi.org/10.1186/s12859-016-0930-z

## Contributing

If you have a suggestion that would make this project better, please send an e-mail or fork the repo and create a pull request.

### Contact

Piotr Kucharski - soliareofastorauj@gmail.com \
Valentyn Bezshapkin - valentyn.bezshapkin@micro.biol.ethz.ch
