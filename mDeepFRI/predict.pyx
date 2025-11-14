import json

import numpy as np

cimport numpy as np

np.import_array()
import gzip
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
    cdef float[:, ::1] onehot_view = np.zeros((seq_len, vocab_size), dtype=np.float32)

    for i in range(seq_len):
        j = chars.find(seq_bytes[i])
        if j != -1:
            onehot_view[i, j] = 1
        else:
            raise ValueError(f"Invalid character in sequence: {seq[i]}")

    return np.asarray(onehot_view)


cdef class Predictor(object):
    """
    Class for loading trained models and computing GO/EC predictions
    and class activation maps (CAMs).
    """

    cdef public str model_path
    cdef public int threads
    cdef public dict prot2goterms
    cdef public dict goidx2chains
    cdef public np.ndarray gonames
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
        meta_path = self.model_path.rsplit(".", 1)[0] + "_model_params.json"
        # read first two bytes
        with open(meta_path, 'rb') as f:
            sig = f.read(2)

        if sig == b"\x1f\x8b":  # gzip magic number
            with gzip.open(meta_path, "rt", encoding="utf-8") as json_file:
                metadata = json.load(json_file)
        else:
            with open(meta_path, "r", encoding="utf-8") as json_file:
                metadata = json.load(json_file)

        self.gonames = np.asarray(metadata['gonames'])

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
