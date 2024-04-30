Installation
============

.. hint::

    Windows OS is not supported. Consider using a Python
    inside `Windows Subsystem for Linux <https://learn.microsoft.com/en-us/windows/wsl/install>`_.

1. Download environment YAML.

.. code:: console

    wget https://raw.githubusercontent.com/bioinf-mcb/Metagenomic-DeepFRI/main/environment.yml

2. Setup conda environment and activate it.

.. code:: console

    conda env create --name deepfri --file environment.yml
    conda activate deepfri
    # Optional cleanup
    rm environment.yml

3. Run and view the help message.

.. code:: console

    mDeepFRI --help
