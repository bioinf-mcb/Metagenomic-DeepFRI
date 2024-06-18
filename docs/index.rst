Welcome to Metagenomic-DeepFRI's documentation! |Stars|
=======================================================

.. |Stars| image:: https://img.shields.io/github/stars/bioinf-MCB/Metagenomic-DeepFRI.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/stargazers


Overview
--------

Protein function is dependent upon 3-dimensional conformation, often likened to a lock and key relationship.
Metagenomic-DeepFRI is a pioneering tool that brings the power of structure-based functional annotation to metagenomic data, enabling researchers to uncover the functional potential of microbial communities.

Key Features
------------

1. **Scalable structure processing:** Leveraging the compression capabilities of `FoldComp <https://github.com/steineggerlab/foldcomp/>`_, Metagenomic-DeepFRI efficiently handles massive datasets of predicted protein structures.
2. **Structural template identification:** The tool utilizes `MMSeqs2 <https://github.com/soedinglab/MMseqs2>`_ to identify the best template protein for functional annotation.
3. **Functional prediction:** Building upon the `DeepFRI <https://github.com/flatironinstitute/DeepFRI>`_ framework, Metagenomic-DeepFRI predicts protein function with high accuracy, enabling inference of the functional capabilities of proteins in metagenomic datasets.

By integrating these features, Metagenomic-DeepFRI offers a powerful solution for functional annotation of metagenomic data, facilitating an understanding of microbial communities and their roles in various ecosystems.

**Metagenomic-DeepFRI** is a joint project between research institutions:

* `Structural and Functional Genomics Lab, Sano Centre for Computational Medicine, Krakow, Poland <https://www.tomaszlab.org>`_
* `Systems Biology Research Group, Center for Computational Biology, Flatiron Institute, New York, USA <https://www.simonsfoundation.org/flatiron/center-for-computational-biology/>`_

Setup
-----

1. Download environment YAML.

.. code:: console

    wget https://raw.githubusercontent.com/bioinf-mcb/Metagenomic-DeepFRI/main/environment.yml

2. Setup conda environment and activate it.

.. code:: console

    conda env create --name deepfri --file environment.yml
    conda activate deepfri
    # Optional cleanup
    rm environment.yml

Contents
--------

.. toctree::
   :maxdepth: 2

    Installation <install>
    Commands <cli>
    API Reference <api/index>
    Examples <examples/index>
    Result Interpretation <results>


.. list-table::

    * - .. image :: _images/sano_logo.png
          :align: center
          :width: 400
      - .. image :: _images/flatiron_logo.png
          :align: center
          :width: 400
    * - .. image :: _images/full_col2x-100.jpg
          :align: center
          :width: 400
