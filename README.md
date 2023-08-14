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
* [FoldComp](https://github.com/steineggerlab/foldcomp)
* [pyOpal](https://github.com/steineggerlab/foldcomp)
* [ONNX](https://github.com/onnx/onnx)

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
## Prepare structural database
Download the database from the [website](https://foldcomp.steineggerlab.workers.dev/). The app was tested with `afdb_swissprot_v4`. You can use different databases, but be mindful that computation time might increase exponentially with the size of the database and the format of protein names might differ and app will crash.
## 1. Download models
Run command:
```
mDeepFRI get-models --output path/to/weights/folder
```

## 2. Predict protein function & capture log
```
mDeepFRI predict-function -i /path/to/protein/sequences -d /path/to/foldcomp/database/ -w /path/to/deepfri/weights/folder -o /output_path 2> log.txt
```

Logging module writes output into `stderr`, so use `2>` to redirect it to file.
Other available parameters can be found upon command `mDeepFRI --help`.
**Attention:** currently predictions are accumulated in memory, so to preserve original DeepFRI code. Therefore the amount of memory needed for the execution of pipeline grows with input size. We advise to split input into 10k protein chunks, thus memory usage is roughly 40GB RAM per chunk.
## Results
Finished folder will contain:
1. `mmseqs2_search_results.m8`
2. `metadata_skipped_ids_due_to_length` - too long or too short queries
3. `queryDB` + index from MMSeqs2 search.
4. `results*.tsv` - multiple files from DeepFRI. Organized by model type (`GCN` or `CNN`) and its mode (`mf`, `bp`, `cc`, `ec`).
Sometimes results from one model can be missing which means that all query proteins sequences were aligned correctly or none of them were aligned.
   ```
   mf = molecular_function
   bp = biological_process
   cc = cellular_component
   ec = enzyme_commission
   ```
## Prediction modes
The GO ontology contains three subontologies, defined by their root nodes:
- Molecular Function (MF)
- Biological Process (BP)
- Cellular Component (CC)
Additionally, Metagenomic-DeepFRI is able to predict Enzyme Comission number (EC). By default, the tool makes prediction in all 4 categories. To select only few pass the parameter `-p` or `--processing-modes` few times, i.e.:
```
mDeepFRI predict-function -i /path/to/protein/sequences -d /path/to/foldcomp/database/ -w /path/to/deepfri/weights/folder -o /output_path -p mf -p bp
```

## Temporary files
The first run of `mDeepFRI` with the database will create temporary files, needed for the pipeline. If you don't want to keep them for the next run use
flag `--no-keep-temporary`.

## CPU / GPU utilization
If argument `threads` is provided, the app will parallelize certain steps (alignment, contact map alignment, inference).
If CUDA is installed on your machine, `mDeepFRI` will automatically use it for prediction. If not, the model will use CPUs.
**Attention:** Single instance of DeepFRI on GPU requires 2GB VRAM.

## Citations
If you use this software please cite:
- Gligorijević et al. "Structure-based protein function prediction using graph convolutional networks" Nat. Comms. (2021). https://doi.org/10.1038/s41467-021-23303-9
- Steinegger & Söding "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets" Nat. Biotechnol. (2017) https://doi.org/10.1038/nbt.3988
- Kim, Midrita & Steinegger "Foldcomp: a library and format for compressing and indexing large protein structure sets" Bioinformatics (2023) https://doi.org/10.1093/bioinformatics/btad153
- Maranga et al. "Comprehensive Functional Annotation of Metagenomes and Microbial Genomes Using a Deep Learning-Based Method" mSystems (2023) https://doi.org/10.1128/msystems.01178-22

## Contributing

If you have a suggestion that would make this project better, please send an e-mail or fork the repo and create a pull request.

### Contact

Valentyn Bezshapkin - valentyn.bezshapkin@micro.biol.ethz.ch \
Piotr Kucharski - soliareofastorauj@gmail.com
