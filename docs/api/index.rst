API Reference
=============

.. toctree::
   :hidden:

    mmseqs <mmseqs>
    pipeline <pipeline>
    database <database>
    alignment <alignment>
    predict <predict>
    bio_utils <bio_utils>
    utils <utils>

.. currentmodule:: mDeepFRI

.. automodule:: mDeepFRI

.. only:: html

    Pipeline
    --------

    Main prediction pipeline for protein function annotation.

    .. autosummary::
        :nosignatures:

        mDeepFRI.pipeline.predict_protein_function

    Database
    --------

    Structure database handling and management.

    .. autosummary::
        :nosignatures:

        mDeepFRI.database.Database

    Alignment
    ---------

    Sequence-structure alignment using PyOpal.

    .. autosummary::
        :nosignatures:

        mDeepFRI.alignment.AlignmentResult
        mDeepFRI.alignment.align_mmseqs_results
        mDeepFRI.alignment.pairwise_against_database
        mDeepFRI.alignment.align_pairwise

    MMSeqs
    ------

    MMseqs2 database search functionality.

    .. autosummary::
        :nosignatures:

        mDeepFRI.mmseqs.QueryFile
        mDeepFRI.mmseqs.MMSeqsSearchResult

    Prediction
    ----------

    Cython-accelerated contact map alignment and prediction.

    .. autosummary::
        :nosignatures:

        mDeepFRI.predict.Predictor
        mDeepFRI.predict.aligned_predict
        mDeepFRI.contact_map_utils.chain_res_to_bins
        mDeepFRI.contact_map_utils.contact_map_from_pdb

    Utilities
    ---------

    General utility functions.

    .. autosummary::
        :nosignatures:

        mDeepFRI.utils.download_file
        mDeepFRI.utils.download_trained_models

    Bio Utilities
    -------------

    Biological data processing utilities.

    .. autosummary::
        :nosignatures:

        mDeepFRI.bio_utils.build_align_contact_map
        mDeepFRI.bio_utils.calculate_contact_map
