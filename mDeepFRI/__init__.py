import os

from mDeepFRI.mmseqs import QueryFile

__version__ = "1.1.8"
__author__ = "Valentyn Bezshapkin <valentyn.bezshapkin@micro.biol.ethz.ch>"
__licencse__ = "BSD-3-Clause"

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


terms = ["bp", "cc", "mf", "ec"]

cnn_model_links = make_links(repo_url, "DeepCNN-MERGED", terms)

gcn_model_links = {
    "1.0":
    make_links(repo_url,
               "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0",
               terms),
    "1.1":
    make_links(
        repo_url,
        "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc",
        ["bp", "cc", "mf"])
}

BAR_FORMAT = "{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}], {rate_fmt}{postfix}"
OUTPUT_HEADER = [
    'Protein', 'GO_term/EC_number', 'Score', 'Annotation', 'Neural_net',
    'DeepFRI_mode', 'DB_hit', 'DB_name', 'Identity', 'Coverage'
]

__all__ = [
    "QueryFile",
    "MMSeqsResult",
]
