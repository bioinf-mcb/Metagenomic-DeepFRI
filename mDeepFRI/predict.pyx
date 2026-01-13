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
    cdef public str model_path
    cdef public int threads
    cdef public object session
    cdef public list input_names

    def __init__(self, model_path: str, threads: int = 1, ):
        self.model_path = model_path
        self.threads = threads

        self._load_model()

    def _load_model(self):
        session_options = rt.SessionOptions()
        session_options.intra_op_num_threads = self.threads
        session_options.inter_op_num_threads = self.threads

        self.session = rt.InferenceSession(
            self.model_path,
            providers=[('CUDAExecutionProvider', {"cudnn_conv_algo_search": "DEFAULT"}),
                        'CPUExecutionProvider'],
            sess_options=session_options,
        )
        self.input_names = [node.name for node in self.session.get_inputs()]

    def forward_pass(self, seqres: str, cmap = None):
        cdef np.ndarray A
        cdef np.ndarray prediction
        cdef np.ndarray y
        cdef dict inputs

        S = seq2onehot(seqres)
        S = S.reshape(1, *S.shape)
        # if cmap present use GCN with 2 inputs - sequence + cmap
        if cmap is not None:
            # GCN Branch
            A = cmap.reshape(1, cmap.shape[0], cmap.shape[1])
            inputs = {
                self.input_names[0]: A.astype(np.float32),
                self.input_names[1]: S.astype(np.float32)
            }
        else:
            # CNN Branch
            inputs = {
                self.input_names[0]: S.astype(np.float32)
            }

        # Run Inference
        prediction = self.session.run(None, inputs)[0]

        y = prediction[:, :, 0].reshape(-1)

        return y
