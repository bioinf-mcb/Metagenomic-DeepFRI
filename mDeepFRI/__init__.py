__version__ = "1.1.1"
cnn_model_links = {
    "bp":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_bp.onnx",
    "cc":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_cc.onnx",
    "mf":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_mf.onnx",
    "ec":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_ec.onnx",
}

gcn_model_links = {
    "bp":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_bp.onnx",
    "cc":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_cc.onnx",
    "mf":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_mf.onnx",
    "ec":
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ec.onnx",
}

config_links = [
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_bp_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_cc_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_mf_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepCNN-MERGED_ec_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_bp_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_cc_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ec_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/DeepFRI-MERGED_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_mf_model_params.json",
    "https://huggingface.co/crusher083/mDeepFRI/resolve/main/model_config.json"
]

TARGET_MMSEQS_DB_NAME = "targetDB"

MERGED_SEQUENCES = "merged_sequences.faa"
MMSEQS_SEARCH_RESULTS = "mmseqs2_search_results.m8"
