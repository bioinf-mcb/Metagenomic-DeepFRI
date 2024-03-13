import os

__version__ = "1.1.3"

repo_url = "https://huggingface.co/valentynbez/mDeepFRI/resolve/main/"

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
