import random
import numpy as np


def random_cmap(seq):
    size = np.product((len(seq), len(seq)))
    board = np.zeros(size, dtype=np.double)
    i = np.random.choice(np.arange(size),  len(seq) * 10)
    board[i] = 1
    return board.reshape((len(seq), len(seq)))


def random_seq():
    return "SHVEAEKQRREKLNHRFYALRAIVPKVSRMDKASLLSDAVSYIESLKSKIDDLETEIKKMKM"
    #return "".join(random.choice("EALWDNVECNRHMLSRYINPAKLTPYLRQCKVIDEQDEDEVLNAPMLPSKINRAGRLLDILHTKGQRGYVVFLESLEFYYPELYKLVT") for _ in range(50))



import json


with open("/data/trained_models/model_config.json") as json_file:
    params = json.load(json_file)["gcn"]['models']


from DeepFRI.deepfrier import Predictor

a = Predictor.Predictor(params["mf"], gcn=True)



seq = random_seq()
cmap = b


a.predict_with_cmap(seq, cmap, "chain_id13")

a.export_csv("somethinf.csv",True)
a.save_predictions("somethinf.json")





from DeepFRI.deepfrier.utils import seq2onehot

import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np

model = tf.keras.models.load_model("/data/trained_models/lstm_lm_tf.hdf5")

seq = "SHVEAEKQRREKLNHRFYALRAIVPKVSRMDKASLLSDAVSYIESLKSKIDDLETEIKKMKM"
seq = seq2onehot(seq)
seq = seq.reshape(1, *seq.shape)
res = model(seq)
res = tf.squeeze(res)


b = np.zeros((res.shape[0], res.shape[0]))
for i in range(res.shape[0]):
    for j in range(res.shape[0]):
        b[i,j] = np.sqrt(np.sum(np.power([res[i] - res[j]],2)))
plt.imshow(b)
plt.title("Python loop - euclidean distance")
plt.colorbar()
plt.show()


b = res @ tf.transpose(res)
plt.imshow(b)
plt.title("matrix @ matrix.T")
plt.colorbar()
plt.show()

from CPP_lib.libAtomDistanceIO import load_contact_map, initialize
initialize()
b = load_contact_map('/home/soliareofastora/data/atoms_dataset/positions/249_310_5gnj.4.B_607fc1be6bed59a94fd9c321.bin', 6.0)
plt.imshow(b)
plt.title("real")
plt.colorbar()
plt.show()

