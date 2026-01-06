Database
========

.. currentmodule:: mDeepFRI.database

.. autoclass:: Database
   :members:
   :undoc-members:
   :show-inheritance:

Description
-----------

The Database class handles structure database management for protein function prediction.
It supports both FoldComp-compressed structure databases and standard PDB/mmCIF formats.

Key Features
------------

- **FoldComp Integration**: Efficiently decompress and access protein structures
- **MMseqs2 Compatibility**: Works with MMseqs2 database search results
- **Structure Caching**: Optimizes repeated structure access
- **Multi-format Support**: Handles PDB, mmCIF, and FoldComp formats

Usage
-----

The Database class is typically used internally by the pipeline to retrieve
structure coordinates for alignment and contact map generation.

Example
-------

.. code-block:: python

    from mDeepFRI.database import Database

    # Initialize database
    db = Database("pdb100.mmseqsDB")

    # Access structure by ID
    structure = db.get_structure("1abc_A")
    coords = structure.get_coordinates()

    # Get structure information
    chain_info = db.get_chain_info("1abc_A")
