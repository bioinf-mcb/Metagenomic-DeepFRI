Alignment
=========

.. currentmodule:: mDeepFRI.alignment

Classes
-------

.. autoclass:: AlignmentResult
   :members:
   :undoc-members:
   :show-inheritance:

Functions
---------

.. autofunction:: insert_gaps

.. autofunction:: best_hit_database

.. autofunction:: align_mmseqs_results

.. autofunction:: pairwise_against_database

.. autofunction:: align_pairwise

Description
-----------

The alignment module provides sequence-structure alignment functionality using PyOpal,
a fast SIMD-accelerated pairwise alignment library.

Key Features
------------

- **PyOpal Integration**: High-performance SIMD-accelerated alignment
- **Custom Scoring Matrices**: Support for BLOSUM and custom matrices
- **Batch Processing**: Efficiently align multiple query-target pairs
- **Detailed Statistics**: Provides identity, coverage, and alignment coordinates

Alignment Workflow
------------------

1. **Database Search**: MMseqs2 identifies candidate structure templates
2. **Sequence Alignment**: PyOpal performs pairwise alignment
3. **Statistics Calculation**: Computes identity, coverage, and quality metrics
4. **Coordinate Mapping**: Maps aligned residues to structure coordinates

Example
-------

.. code-block:: python

    from mDeepFRI.alignment import align_pairwise
    from mDeepFRI.database import Database

    # Initialize database
    db = Database("pdb100.mmseqsDB")

    # Perform alignment
    result = align_pairwise(
        query_seq="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL",
        target_seq="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL",
        target_coords=db.get_structure("1abc_A").coords,
        scoring_matrix="BLOSUM62"
    )

    print(f"Identity: {result.query_identity:.1f}%")
    print(f"Coverage: {result.query_coverage:.1f}%")
    print(f"Aligned coordinates: {len(result.coords)}")

Scoring Matrices
----------------

Supported scoring matrices include:

- **BLOSUM62** (default): Balanced for diverse sequences
- **BLOSUM45**: More permissive for distant homologs
- **BLOSUM80**: More stringent for close homologs
- **PAM250**: Alternative evolutionary model

Custom scoring matrices can be provided as dictionaries.
