import json

import numpy as np

cimport numpy as np

np.import_array()

import operator

import cython
import onnxruntime as rt


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef np.ndarray[float, ndim=2] seq2onehot(str seq):
    """
    Converts a protein sequence to 26-dim one-hot encoding.

    Args:
        seq (str): protein sequence

    Returns:
        np.ndarray: one-hot encoding of the protein sequence
    """
    cdef bytes seq_bytes = seq.encode()

    cdef bytes chars = b"-DGULNTKHYWCPVSOIEFXQABZRM"
    cdef int vocab_size = len(chars)
    cdef int seq_len = len(seq)
    cdef int i, j
    cdef int[:, ::1] onehot_view = np.zeros((seq_len, vocab_size), dtype=np.int32)

    for i in range(seq_len):
        j = chars.find(seq_bytes[i])
        if j != -1:
            onehot_view[i, j] = 1
        else:
            raise ValueError(f"Invalid character in sequence: {seq[i]}")

    return np.asarray(onehot_view, dtype=np.float32)


cdef class Predictor(object):
    """
    Class for loading trained models and computing GO/EC predictions and class activation maps (CAMs).
    """

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

    def __init__(self, model_path: str, threads: int = 0, ):
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

    def forward_pass(self, seqres: str, cmap = None):

        cdef np.ndarray A
        cdef np.ndarray prediction
        cdef np.ndarray y

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

        # if no cmap use CNN with 1 input
        else:
            prediction = self.session.run(
                None, {inputDetails[0].name: S.astype(np.float32)})[0]

        y = prediction[:, :, 0].reshape(-1)

        return y

    def format_predictions(self, y, chain: str = ""):

        cdef list output_rows = []
        cdef str go_term
        cdef float score
        cdef str annotation

        go_idx = np.where(y >= self.thresh)[0]
        for idx in go_idx:
            go_term = self.goterms[idx].item()
            score = float(y[idx])
            annotation = self.gonames[idx].item()
            output_rows.append([chain, go_term, score, annotation])

        # Sort output_rows based on score in descending order
        output_rows.sort(key=operator.itemgetter(2), reverse=True)

        return output_rows


    def predict_function(
        self,
        seqres: str,
        cmap = None,
        chain: str = "",
        format: str = "tsv",
    ):
        """
        Computes GO/EC predictions for a single protein chain from sequence and contact map.

        Args:
            seqres (str): protein sequence.
            cmap (np.array): contact map.
            chain (str): protein ID.
            format (str): output format. Options: "tsv", "vector".

        Returns:
            list: list of GO/EC predictions.

        """

        y = self.forward_pass(seqres, cmap)
        output_rows = self.format_predictions(y, chain)

        return output_rows
