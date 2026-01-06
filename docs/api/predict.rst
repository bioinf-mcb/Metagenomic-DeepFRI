Predict
=======

.. currentmodule:: mDeepFRI.predict

Classes
-------

.. autoclass:: Predictor
   :members:
   :undoc-members:
   :show-inheritance:

Functions
---------

.. autofunction:: aligned_predict

Description
-----------

The predict module provides Cython-accelerated functionality for contact map
alignment and DeepFRI model prediction. This module contains performance-critical
code that has been optimized using Cython for speed.

Key Components
--------------

Predictor Class
^^^^^^^^^^^^^^^

The :class:`Predictor` class loads and manages DeepFRI ONNX models for protein
function prediction. It supports multiple prediction modes:

- **GCN (Graph Convolutional Network)**: Uses both sequence and structure
- **CNN (Convolutional Neural Network)**: Uses sequence only
- **Modes**: Molecular Function (mf), Biological Process (bp), Cellular Component (cc), EC numbers (ec)

Contact Map Utilities
^^^^^^^^^^^^^^^^^^^^^

The module also provides utilities from :mod:`mDeepFRI.contact_map_utils`:

- :func:`contact_map_from_pdb`: Generate contact maps from PDB coordinates
- :func:`chain_res_to_bins`: Convert residue coordinates to distance bins

Aligned Prediction
^^^^^^^^^^^^^^^^^^

The :func:`aligned_predict` function performs contact map alignment and prediction:

1. Generates contact map from template structure
2. Aligns contact map to query sequence
3. Runs DeepFRI prediction on aligned contact map
4. Returns GO term/EC number predictions with scores

Example
-------

.. code-block:: python

    from mDeepFRI.predict import Predictor, aligned_predict
    from mDeepFRI.alignment import AlignmentResult

    # Initialize predictor
    predictor = Predictor(
        model_path="models/v1.1/",
        gcn=True,  # Use graph convolutional network
        device="cpu"
    )

    # Perform aligned prediction
    predictions = aligned_predict(
        predictor=predictor,
        query_seq="MKTAYIAKQRQISFVK...",
        alignment_result=alignment_result,
        mode="mf"  # Molecular function
    )

    # Process predictions
    for go_term, score in predictions:
        if score > 0.5:
            print(f"{go_term}: {score:.3f}")

Performance Notes
-----------------

This module uses Cython for performance-critical operations:

- Contact map generation is ~10x faster than pure Python
- Batch prediction supports parallel processing
- ONNX Runtime provides optimized inference

Type Stubs
----------

Type annotations are provided in `.pyi` stub files for IDE support and type checking.
