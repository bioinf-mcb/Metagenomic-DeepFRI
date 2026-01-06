Utilities
=========

.. currentmodule:: mDeepFRI.utils

Functions
---------

.. autofunction:: download_file

.. autofunction:: download_trained_models

Description
-----------

The utils module provides general utility functions for file handling,
model management, and common operations.

Key Features
------------

- **File Download**: Download files with progress tracking
- **Model Management**: Download and manage pre-trained DeepFRI models
- **Path Handling**: Manage file paths and directories
- **Configuration**: Handle configuration files and settings

Model Download
--------------

The :func:`download_trained_models` function automatically downloads
pre-trained DeepFRI models from the official repository:

.. code-block:: python

    from mDeepFRI.utils import download_trained_models

    # Download models to default location
    download_trained_models()

    # Download to custom directory
    download_trained_models(target_dir="custom/models/")

File Download
-------------

The :func:`download_file` function provides robust file downloading
with progress tracking:

.. code-block:: python

    from mDeepFRI.utils import download_file

    # Download with progress bar
    download_file(
        url="https://example.com/database.tar.gz",
        output_path="database.tar.gz",
        show_progress=True
    )

Example Usage
-------------

.. code-block:: python

    from mDeepFRI.utils import download_trained_models, download_file
    import os

    # Setup project directory
    models_dir = "models/v1.1"
    os.makedirs(models_dir, exist_ok=True)

    # Download pre-trained models
    download_trained_models(target_dir=models_dir)

    # Download additional database
    download_file(
        url="https://example.com/pdb100.tar.gz",
        output_path="pdb100.tar.gz"
    )
