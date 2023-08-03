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

        # Not clear how parameter influences GPU exec
        if rt.get_device() == 'CPU':
            self.threads = threads
        elif rt.get_device() == 'GPU':
            self.threads = 0

        self._load_model()
        self.prot2goterms = {}
        self.goidx2chains = {}

    def _load_model(self):
        session_options = rt.SessionOptions()
        session_options.intra_op_num_threads = self.threads
        self.session = rt.InferenceSession(
            self.model_path,
            providers=['CUDAExecutionProvider', 'CPUExecutionProvider'],
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
        self.prot2goterms[chain] = []
        go_idx = np.where(y >= self.thresh)[0]
        for idx in go_idx:
            if idx not in self.goidx2chains:
                self.goidx2chains[idx] = set()
            self.goidx2chains[idx].add(chain)
            self.prot2goterms[chain].append(
                (self.goterms[idx], self.gonames[idx], float(y[idx])))


    def __save_file(self,
                    output_fn: str,
                    delimiter: str,
                    quotechar: str = '"'):
        """
        Exports predictions to .csv or .tsv format

        Args:
            output_fn (str): output file name
            delimiter (str): delimiter for .csv or .tsv file
            quotechar (str): quotechar for .csv file

        Returns:
            None
        """

        with open(output_fn, 'w', encoding="utf-8") as f:
            writer = csv.writer(f, delimiter=delimiter, quotechar=quotechar)
            writer.writerow([
                'Protein', 'GO_term/EC_number', 'Score',
                'GO_term/EC_number name'
            ])
            for prot, goterms in self.prot2goterms.items():
                sorted_rows = sorted(goterms, key=lambda x: x[2], reverse=True)
                for row in sorted_rows:
                    writer.writerow(
                        [prot, row[0], '{:.5f}'.format(row[2]), row[1]])

    def export_tsv(self, output_fn: str):
        """
        Exports predictions to .tsv format
        """

        self.__save_file(output_fn, "\t")
