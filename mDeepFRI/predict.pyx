import json

import numpy as np

cimport numpy as np

np.import_array()

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

    cdef bytes seq_bytes = seq.encode('ascii')
    cdef const unsigned char[:] seq_view = seq_bytes
    cdef int seq_len = len(seq_bytes)

    cdef float[:, ::1] onehot_view = np.zeros((seq_len, 26), dtype=np.float32)

    cdef int[256] char_map
    cdef bytes chars = b"-DGULNTKHYWCPVSOIEFXQABZRM"
    cdef int i, code
    cdef int error_index = -1

    for i in range(256):
        char_map[i] = -1
    for i in range(26):
        char_map[<unsigned char>chars[i]] = i

    with cython.nogil:
        for i in range(seq_len):
            code = char_map[seq_view[i]]

            if code != -1:
                onehot_view[i, code] = 1.0
            else:
                error_index = i
                break

    if error_index != -1:
        raise ValueError(f"Invalid character in sequence: {seq[error_index]}")

    return np.asarray(onehot_view)

cdef class Predictor(object):
    """
    Class for loading trained models and computing GO/EC predictions
    and class activation maps (CAMs).
    """

    cdef public str model_path
    cdef public int threads
    cdef public object session

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

    def forward_pass(self, seqres: str, cmap = None):

        cdef np.ndarray A
        cdef np.ndarray prediction
        cdef np.ndarray y
        cdef dict inputs

        S = seq2onehot(seqres)
        S = S.reshape(1, *S.shape)
        inputs = self.session.get_inputs()
        # if cmap present use GCN with 2 inputs - sequence + cmap
        if cmap is not None:
            # GCN Branch
            A = cmap.reshape(1, cmap.shape[0], cmap.shape[1])
            inputs = {
                self.input_name_1: A.astype(np.float32),
                self.input_name_2: S.astype(np.float32)
            }
        else:
            # CNN Branch
            inputs = {
                self.input_name_1: S.astype(np.float32)
            }

        # Run Inference
        prediction = self.session.run(None, inputs)[0]

        y = prediction[:, :, 0].reshape(-1)

        return y
