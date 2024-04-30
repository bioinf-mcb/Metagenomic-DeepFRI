
Welcome to Metagenomic-DeepFRI's documentation!
===============================================
Do you have **thousands of protein sequences** with **unknown structures**, but still want to know their
molecular function, biological process, cellular component and enzyme commission **predicted by DeepFRI?**
Metagenomic-DeepFRI is a deep learning-based software for functional annotation of metagenomes.
It combines both protein sequence and structure information for function prediction. It is a joint project of the `Systems Biology Research Group, Center for Computational Biology, Flatiron Institute, New York, USA <https://www.simonsfoundation.org/flatiron/center-for-computational-biology/>`_
and `Structural and Functional Genomics Lab, Małopolska Center of Biotechnology UJ, Krakow, Poland <https://www.tomaszlab.org>`_.

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
    API Reference <api/index>
    Examples <examples/index>

.. list-table::

    * - .. image :: _images/mcb_logo_en_v2.jpg
          :align: center
      - .. image :: _images/flatiron_logo.png
          :align: center
    * - .. image :: _images/full_col2x-100.jpg
          :align: center
          :width: 70 %
      -
