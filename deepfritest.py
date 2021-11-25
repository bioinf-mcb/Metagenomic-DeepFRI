import json
import shutil
import time

from Bio import pairwise2
from Bio import SeqIO

from CONFIG import *
from CPP_lib.libAtomDistanceIO import load_contact_map, initialize
from DeepFRI.deepfrier import Predictor
from utils import add_path_to_env
from run_mmseqs_search import run_mmseqs_search
# chromwell_process_fasta.py looks like the type of script i need to create

import matplotlib.pyplot as plt


initialize()
add_path_to_env(FOLDSEEK_BIN_PATH)
pipeline_start = str(time.time())
faa_files = list(QUERY_PATH.glob("**/*.faa"))

work_path = (WORK_PATH / pipeline_start)
while work_path.exists():
    pipeline_start = str(int(time.time()))
    work_path = (WORK_PATH / pipeline_start)
work_path.mkdir()

query_data = dict()
for query_file in faa_files:
    query_data.update({record.id: record.seq for record in SeqIO.parse(query_file, "fasta")})

with open(work_path / 'merged_sequences.faa', 'wb') as writer:
    for seq_file in faa_files:
        with open(seq_file, 'rb') as reader:
            shutil.copyfileobj(reader, writer)

output_data = run_mmseqs_search(work_path / 'merged_sequences.faa', work_path)

# shutil.rmtree(work_path)

# with open("/data/trained_models/model_config.json") as json_file:
#     params = json.load(json_file)
# gcn_params = params["gcn"]["models"]
# gcn = Predictor.Predictor(gcn_params["mf"], gcn=True)

for i in range(len(output_data)):
    query_id = output_data["query"].iloc[i]
    query_sequence = query_data[query_id]
    target_id = output_data["target"].iloc[i]

    with open(ATOMS_DATASET_PATH / "seq" / (target_id + ".faa")) as f:
        for record in SeqIO.parse(f, "fasta"):
            target_sequence = record.seq

    target_cmap = load_contact_map(str(ATOMS_DATASET_PATH / "positions" / (target_id + ".bin")),
                                   ANGSTROM_CONTACT_THRESHOLD)

    query_cmap = load_contact_map(str(ATOMS_DATASET_PATH / "positions" / (query_id + ".bin")),
                                   ANGSTROM_CONTACT_THRESHOLD)

    # f = plt.figure()
    # f.add_subplot(1, 2, 1)
    # plt.imshow(target_cmap)
    # plt.title("target")
    # f.add_subplot(1, 2, 2)
    # plt.imshow(query_cmap)
    # plt.title("query\n" + query_id)
    # plt.show(block=True)

    print(query_id)
    alignment = pairwise2.align.globalxx(query_sequence, target_sequence)[0]
    from Bio.pairwise2 import format_alignment
    print(format_alignment(*alignment))


#         indexses = []
#         if alignment.score == len(query_sequence):
#             if len(query_sequence) == len(alignment.seqA):
#                 aligned_cmap = target_cmap
#             else:
#                 for j in range(len(alignment.seqA)):
#                     if alignment.seqA[j] != '-':
#                         indexses.append(j)
#                 aligned_cmap = target_cmap[indexses]
#                 aligned_cmap = aligned_cmap[:, indexses]
#         else:
#             # print(i)
#             # from Bio.pairwise2 import format_alignment
#             # print(format_alignment(*alignment))
#             continue
#
#         gcn.predict_with_cmap(query_sequence, aligned_cmap, query_id)
#
#     gcn.export_csv("somethinf.csv", True)
#     gcn.save_predictions("somethinf.json")
#
#
# if __name__ == '__main__':
#     main_pipeline()
#
#


# import random
# import numpy as np
#
#
# def random_cmap(seq):
#     size = np.product((len(seq), len(seq)))
#     board = np.zeros(size, dtype=np.double)
#     i = np.random.choice(np.arange(size),  len(seq) * 10)
#     board[i] = 1
#     return board.reshape((len(seq), len(seq)))
#
#
# def random_seq():
#     return "SHVEAEKQRREKLNHRFYALRAIVPKVSRMDKASLLSDAVSYIESLKSKIDDLETEIKKMKM"
#     #return "".join(random.choice("EALWDNVECNRHMLSRYINPAKLTPYLRQCKVIDEQDEDEVLNAPMLPSKINRAGRLLDILHTKGQRGYVVFLESLEFYYPELYKLVT") for _ in range(50))
#
#
#
# import json
#
#
# with open("/data/trained_models/model_config.json") as json_file:
#     params = json.load(json_file)["gcn"]['models']
#
#
# from DeepFRI.deepfrier import Predictor
#
# a = Predictor.Predictor(params["mf"], gcn=True)
#
# seq = random_seq()
# cmap = b
#
#
# a.predict_with_cmap(seq, cmap, "chain_id13")
#
# a.export_csv("somethinf.csv",True)
# a.save_predictions("somethinf.json")
#
#
#
#
#
# from DeepFRI.deepfrier.utils import seq2onehot
#
# import tensorflow as tf
# import matplotlib.pyplot as plt
# import numpy as np
#
# model = tf.keras.models.load_model("/data/trained_models/lstm_lm_tf.hdf5")
#
# seq = "SHVEAEKQRREKLNHRFYALRAIVPKVSRMDKASLLSDAVSYIESLKSKIDDLETEIKKMKM"
# seq = seq2onehot(seq)
# seq = seq.reshape(1, *seq.shape)
# res = model(seq)
# res = tf.squeeze(res)
#
#
# b = np.zeros((res.shape[0], res.shape[0]))
# for i in range(res.shape[0]):
#     for j in range(res.shape[0]):
#         b[i,j] = np.sqrt(np.sum(np.power([res[i] - res[j]],2)))
# plt.imshow(b)
# plt.title("Python loop - euclidean distance")
# plt.colorbar()
# plt.show()
#
#
# b = res @ tf.transpose(res)
# plt.imshow(b)
# plt.title("matrix @ matrix.T")
# plt.colorbar()
# plt.show()
#
# from CPP_lib.libAtomDistanceIO import load_contact_map, initialize
# initialize()
# b = load_contact_map('/home/soliareofastora/data/atoms_dataset/positions/249_310_5gnj.4.B_607fc1be6bed59a94fd9c321.bin', 6.0)
# plt.imshow(b)
# plt.title("real")
# plt.colorbar()
# plt.show()
#
