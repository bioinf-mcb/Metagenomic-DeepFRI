import os

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
