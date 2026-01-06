Examples
========

Basic Usage
-----------

Predict protein function with default settings:

.. code-block:: bash

    mDeepFRI predict-function -i /path/to/protein/sequences \
    -d /path/to/foldcomp/database/ \
    -w /path/to/deepfri/weights/folder \
    -o /output_path > log.txt

MMseqs2 Search
--------------

Programmatic database search using the MMseqs module:

.. code-block:: python

    from mDeepFRI.mmseqs import QueryFile

    # Create a query file
    query = QueryFile("sequences.fasta")
    # Load sequences to manipulate on them
    query.load_sequences()
    # filter sequences under 30 amino acids
    query.filter_sequences(min_length = 30)
    # search against MMseqs2 database
    result = query.search("mmseqs.db", eval=10e-3, mmseqs_sensitivity=4, threads=8)
    # save results
    result.save("output.tsv")

Hierarchical Database Search
-----------------------------

Search multiple databases hierarchically (e.g., PDB first, then AlphaFold):

.. code-block:: bash

    mDeepFRI predict-function -i /path/to/protein/sequences \
      -d /path/to/alphafold/database/ -d /path/to/esmfold/database/ \
      -w /path/to/deepfri/weights/folder -o /output_path

Selecting Prediction Modes
---------------------------

Predict only specific ontology categories:

.. code-block:: bash

    mDeepFRI predict-function -i /path/to/protein/sequences \
      -d /path/to/foldcomp/database/ -w /path/to/deepfri/weights/folder \
      -o /output_path -p mf -p bp

Skipping Prediction Matrices
-----------------------------

Save disk space by skipping large matrix files:

.. code-block:: bash

    mDeepFRI predict-function -i /path/to/protein/sequences \
      -d /path/to/foldcomp/database/ -w /path/to/deepfri/weights/folder \
      -o /output_path --skip-matrix

Filtering Results
-----------------

Load and filter prediction results in Python:

.. code-block:: python

    import pandas as pd

    # Load results
    df = pd.read_csv("results.tsv", sep="\t")

    # Filter high-confidence structure-based predictions
    filtered = df[
        (df["score"] >= 0.5) &
        (df["network_type"] == "gcn") &
        (df["aligned"] == True)
    ]

    # Group by protein
    for protein_id, group in filtered.groupby("protein"):
        print(f"\n{protein_id}:")
        for _, row in group.iterrows():
            print(f"  {row['go_term']}: {row['go_name']} (score: {row['score']:.3f})")
