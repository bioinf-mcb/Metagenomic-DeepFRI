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

cnn_model_links = {
    "bp": {
        "model": os.path.join(repo_url, "DeepCNN-MERGED_bp.onnx"),
        "config": os.path.join(repo_url, "DeepCNN-MERGED_bp_model_params.json")
    },
    "cc": {
        "model": os.path.join(repo_url, "DeepCNN-MERGED_cc.onnx"),
        "config": os.path.join(repo_url, "DeepCNN-MERGED_cc_model_params.json")
    },
    "mf": {
        "model": os.path.join(repo_url, "DeepCNN-MERGED_mf.onnx"),
        "config": os.path.join(repo_url, "DeepCNN-MERGED_mf_model_params.json")
    },
    "ec": {
        "model": os.path.join(repo_url, "DeepCNN-MERGED_ec.onnx"),
        "config": os.path.join(repo_url, "DeepCNN-MERGED_ec_model_params.json")
    },
}

gcn_model_links = {
    "1.0": {
        "bp": {
            "model":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_bp.onnx"
            ),
            "config":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_bp_model_params.json"
            )
        },
        "cc": {
            "model":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_cc.onnx"
            ),
            "config":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_cc_model_params.json"
            )
        },
        "mf": {
            "model":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_mf.onnx"
            ),
            "config":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_mf_model_params.json"
            )
        },
        "ec": {
            "model":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ec.onnx"
            ),
            "config":
            os.path.join(
                repo_url,
                "DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ec_model_params.json"
            )
        }
    },
    "1.1": {
        "bp": {
            "model":
            os.path.join(
                repo_url,
                "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_bp.onnx"
            ),
            "config":
            os.path.join(
                repo_url,
                "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_bp_model_params.json"
            )
        },
        "cc": {
            "model":
            os.path.join(
                repo_url,
                "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_cc.onnx"
            ),
            "config":
            os.path.join(
                repo_url,
                "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_cc_model_params.json"
            )
        },
        "mf": {
            "model":
            os.path.join(
                repo_url,
                "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_mf.onnx"
            ),
            "config":
            os.path.join(
                repo_url,
                "DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_mf_model_params.json"
            )
        },
    }
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
