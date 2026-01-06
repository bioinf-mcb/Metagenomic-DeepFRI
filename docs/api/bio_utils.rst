Bio Utilities
=============

.. currentmodule:: mDeepFRI.bio_utils

Functions
---------

.. autofunction:: build_align_contact_map

.. autofunction:: decompress_and_decode

.. autofunction:: get_calpha_coordinates

.. autofunction:: construct_contact_map

.. autofunction:: align_coordinates

.. autofunction:: calculate_contact_map

Description
-----------

The bio_utils module provides biological utilities for structure processing and contact map generation.
Values computed from these utilities are critical inputs for structure-based protein function prediction using DeepFRI.

Key Features
------------

- **Structure Loading**: Extract structures from specific database formats like FoldComp or standard PDB/CIF files
- **Coordinate Extraction**: Isolate C-alpha coordinates crucial for backbone representation
- **Contact Map Generation**: create adjacency matrices representing protein residue interactions
- **Contact Alignment**: Map structural contact information onto query sequences guided by alignments

Usage
-----
This module is primarily used internally by the pipeline to process structural data, but functions can be used independently for structural analysis tasks.

Example
-------

.. code-block:: python

    from mDeepFRI.bio_utils import get_calpha_coordinates, calculate_contact_map

    # Assuming 'structure' is a loaded Biotite AtomArray
    coords = get_calpha_coordinates(structure)
    if coords is not None:
        # Calculate contact map with 6 Angstrom threshold
        cmap = calculate_contact_map(coords, threshold=6.0)
