
Welcome to Metagenomic-DeepFRI's documentation! |Stars|
=======================================================

.. |Stars| image:: https://img.shields.io/github/stars/bioinf-MCB/Metagenomic-DeepFRI.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/stargazers

Overview
--------

Protein function is strongly dependent on 3-dimensional confirmation of the protein, making protein targets resemble a lock and the key.
**Metagenomic-DeepFRI** scales the first structure-based functional annotation method to metagenomic data.
It encorporates massive datasets of predicted structures conviniently compressed by `FoldComp <https://github.com/steineggerlab/foldcomp/>`_, searches for the best
template protein via `MMSeqs2 <https://github.com/soedinglab/MMseqs2>`_, identifies differences between template and target with the help of alignment and predicts protein function
version of `DeepFRI <https://github.com/flatironinstitute/DeepFRI>`_.
It is a joint project of the `Systems Biology Research Group, Center for Computational Biology, Flatiron Institute, New York, USA <https://www.simonsfoundation.org/flatiron/center-for-computational-biology/>`_
and `Structural and Functional Genomics Lab, Sano Centre for Computational Medicine, Krakow, Poland <https://www.tomaszlab.org>`_.

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

Library
-------

.. toctree::
   :maxdepth: 2

    Installation <install>
    Command Line Interface <cli>
    API Reference <api/index>
    Examples <examples/index>

.. list-table::

    * - .. image :: _images/sano_logo.png
          :align: center
      - .. image :: _images/flatiron_logo.png
          :align: center
    * - .. image :: _images/full_col2x-100.jpg
          :align: center
          :width: 70 %
      -
