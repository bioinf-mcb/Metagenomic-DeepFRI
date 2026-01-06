"""
Metagenomic-DeepFRI: High-performance protein function annotation pipeline.

This module provides a pipeline for annotating protein sequences with Gene Ontology (GO)
terms using DeepFRI, a deep learning model for functional protein annotation. It combines
structure information from FoldComp databases, sequence-based predictions using DeepFRI's
neural networks, and fast searches with MMseqs2 for database alignment.

Key Features:
    - Structure information from FoldComp databases (AlphaFold, ESMFold, PDB, etc.)
    - Sequence-based predictions using DeepFRI's neural networks
    - Fast searches with MMseqs2 for database alignment
    - 2-12× speedup compared to standard DeepFRI

Attributes:
    DEEPFRI_MODES (dict): Dictionary mapping prediction modes to their descriptions.
        - "bp": Gene Ontology Biological Process
        - "cc": Gene Ontology Cellular Component
        - "mf": Gene Ontology Molecular Function
        - "ec": Enzyme Commission numbers

Example:
    Basic usage of the pipeline::

        from mDeepFRI.pipeline import hierarchical_database_search
        from mDeepFRI.mmseqs import QueryFile

        query_file = QueryFile(filepath="proteins.fasta")
        hierarchical_database_search(
            query_file=query_file,
            output_path="./results",
            databases=["path/to/database"],
            threads=4
        )

References:
    - DeepFRI: https://github.com/flatironinstitute/DeepFRI
    - FoldComp: https://github.com/steineggerlab/foldcomp
    - MMseqs2: https://github.com/soedinglab/MMseqs2
"""

import os
from importlib.metadata import version

from mDeepFRI.mmseqs import QueryFile

repo_url = "https://huggingface.co/valentynbez/mDeepFRI/resolve/main/"

DEEPFRI_MODES = {
    "bp": "GO Biological Process",
    "cc": "GO Cellular Component",
    "mf": "GO Molecular Function",
    "ec": "Enzyme Commission"
}


def make_links(repo_url, prefix, terms):
    return {
        term: {
            "model": os.path.join(repo_url, f"{prefix}_{term}.onnx"),
            "config": os.path.join(repo_url,
                                   f"{prefix}_{term}_model_params.json"),
        }
        for term in terms
    }


cnn_model_links = make_links(repo_url, "DeepCNN-MERGED", DEEPFRI_MODES.keys())

gcn_model_links = {
    "1.0":
    make_links(repo_url,
               "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0",
               DEEPFRI_MODES.keys()),
    "1.1":
    make_links(
        repo_url,
        "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc",
        ["bp", "cc", "mf"])
}
__all__ = [
    "QueryFile",
    "MMSeqsResult",
]

__author__ = "Valentyn Bezshapkin <valentyn.bezshapkin@micro.biol.ethz.ch>"

__version__ = version(__name__)
