Pipeline
========

.. currentmodule:: mDeepFRI.pipeline

.. autofunction:: predict_protein_function

Description
-----------

The pipeline module provides the main orchestration for protein function prediction.
It coordinates sequence search, alignment, contact map generation, and DeepFRI prediction.

The :func:`predict_protein_function` function is the primary entry point that:

1. Searches query sequences against a structure database using MMseqs2
2. Aligns query sequences to structure templates using PyOpal
3. Generates contact maps from aligned structures
4. Runs DeepFRI prediction using graph convolutional networks
5. Outputs functional annotations with GO terms and EC numbers

Parameters
----------

The function accepts various parameters to control the prediction process:

- **query_fasta**: Path to input FASTA file with query proteins
- **database**: Path to structure database (MMseqs2 format)
- **output_dir**: Directory for output files
- **models_path**: Path to pre-trained DeepFRI models
- **mmseqs_sensitivity**: MMseqs2 search sensitivity (1-7, default: 4)
- **mmseqs_evalue**: E-value threshold for database search
- **skip_matrix**: Skip writing large prediction matrix files
- **scoring_matrix**: Custom scoring matrix for alignment (e.g., BLOSUM62)
- **threads**: Number of parallel threads to use

Output Files
------------

The pipeline generates several output files in the specified output directory:

- **results.tsv**: Main results file with functional annotations
- **alignment_summary.tsv**: Alignment statistics summary
- **prediction_matrix_*.tsv**: Full prediction matrices (optional)
- **database_search/**: MMseqs2 search results

Example
-------

.. code-block:: python

    from mDeepFRI.pipeline import predict_protein_function

    # Run prediction pipeline
    predict_protein_function(
        query_fasta="proteins.fasta",
        database="pdb100.mmseqsDB",
        output_dir="results/",
        models_path="models/v1.1/",
        mmseqs_sensitivity=4,
        mmseqs_evalue=1e-3,
        skip_matrix=True,
        scoring_matrix="BLOSUM62",
        threads=8
    )
