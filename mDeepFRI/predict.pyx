import csv
import json

import numpy as np

cimport numpy as np

np.import_array()

import onnxruntime as rt

from mDeepFRI.bio_utils import seq2onehot


## TODO: do not accumulate predictions, write them out to csv asap to avoid unneccessary RAM usage
cdef class Predictor(object):
    """
    Class for loading trained models and computing GO/EC predictions and class activation maps (CAMs).
    """

    cdef public bint gcn
    cdef public str model_path
    cdef public int threads
    cdef public dict prot2goterms
    cdef public dict goidx2chains
    cdef public np.ndarray gonames
    cdef public np.ndarray goterms
    cdef public np.ndarray thresh
    cdef public session
    cdef public np.ndarray Y_hat
    cdef public dict data
    cdef public list test_prot_list

    def __init__(self, model_path: str, threads: int = 0):
        self.model_path = model_path
        self.threads = threads

        self._load_model()

    def _load_model(self):
        session_options = rt.SessionOptions()
        session_options.intra_op_num_threads = self.threads

        self.session = rt.InferenceSession(
            self.model_path,
            providers=[('CUDAExecutionProvider', {"cudnn_conv_algo_search": "DEFAULT"}),
                        'CPUExecutionProvider'],
            sess_options=session_options,
        )

        # load parameters
        with open(self.model_path.rsplit(".", 1)[0] + "_model_params.json") as json_file:
            metadata = json.load(json_file)

        self.gonames = np.asarray(metadata['gonames'])
        self.goterms = np.asarray(metadata['goterms'])
        self.thresh = 0.1 * np.ones(len(self.goterms))

    def predict_function(
        self,
        seqres: str,
        cmap = None,
        chain: str = "",
    ):
        """
        Computes GO/EC predictions for a single protein chain from sequence and contact map.

        Args:
            seqres (str): protein sequence
            cmap (np.array): contact map
            chain (str): chain ID

        Returns:
            None
        """

        cdef np.ndarray A
        cdef np.ndarray prediction
        cdef np.ndarray y
        cdef list output_rows = []
        cdef str go_term
        cdef float score
        cdef str annotation

        self.Y_hat = np.zeros((1, len(self.goterms)), dtype=float)
        self.data = {}
        self.test_prot_list = [chain]
        S = seq2onehot(seqres)
        S = S.reshape(1, *S.shape)
        inputDetails = self.session.get_inputs()
        # if cmap present use GCN with 2 inputs - sequence + cmap
        if cmap is not None:
            A = cmap.reshape(1, *cmap.shape)
            prediction = self.session.run(
                None, {
                    inputDetails[0].name: A.astype(np.float32),
                    inputDetails[1].name: S.astype(np.float32)
                })[0]
            self.data[chain] = [[A, S], seqres]

        # if no cmap use CNN with 1 input
        else:
            prediction = self.session.run(
                None, {inputDetails[0].name: S.astype(np.float32)})[0]
            self.data[chain] = [[S], seqres]

        y = prediction[:, :, 0].reshape(-1)
        self.Y_hat[0] = y
        go_idx = np.where(y >= self.thresh)[0]
        for idx in go_idx:
            go_term = self.goterms[idx].item()
            score = float(y[idx])
            annotation = self.gonames[idx].item()
            output_rows.append([chain, go_term, score, annotation])

        return output_rows
