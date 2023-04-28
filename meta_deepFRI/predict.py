import csv
import json

import numpy as np
import onnxruntime as rt

from meta_deepFRI.utils.bio_utils import seq2onehot


## TODO: do not accumulate predictions, write them out to csv asap to avoid unneccessary RAM usage
class Predictor(object):
    """
    Class for loading trained models and computing GO/EC predictions and class activation maps (CAMs).
    """
    def __init__(self, model_prefix: str, gcn: bool = True):
        self.model_prefix = model_prefix
        self.gcn = gcn
        self._load_model()
        self.prot2goterms = {}
        self.goidx2chains = {}

    def _load_model(self):
        self.session = rt.InferenceSession(
            self.model_prefix + '.onnx',
            providers=['CUDAExecutionProvider', 'CPUExecutionProvider'],
        )

        # load parameters
        with open(self.model_prefix + "_model_params.json") as json_file:
            metadata = json.load(json_file)

        self.gonames = np.asarray(metadata['gonames'])
        self.goterms = np.asarray(metadata['goterms'])
        self.thresh = 0.1 * np.ones(len(self.goterms))

    def predict_function(
        self,
        seqres: np.ndarray,
        cmap: np.ndarray = None,
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
        self.Y_hat = np.zeros((1, len(self.goterms)), dtype=float)
        self.data = {}
        self.test_prot_list = [chain]
        S = seq2onehot(seqres)
        S = S.reshape(1, *S.shape)
        inputDetails = self.session.get_inputs()
        if self.gcn:
            A = cmap.reshape(1, *cmap.shape)
            prediction = self.session.run(
                None, {
                    inputDetails[0].name: A.astype(np.float32),
                    inputDetails[1].name: S.astype(np.float32)
                })[0]
            self.data[chain] = [[A, S], seqres]
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

    def export_json(self, output_fn: str):
        """
        Exports predictions to .json format

        Args:
            output_fn (str): output file name

        Returns:
            None
        """
        data = []
        for prot, goterms in self.prot2goterms.items():
            sorted_rows = sorted(goterms, key=lambda x: x[2], reverse=True)
            for row in sorted_rows:
                data.append({
                    'Protein': prot,
                    'GO_term/EC_number': row[0],
                    'Score': '{:.5f}'.format(row[2]),
                    'GO_term/EC_number name': row[1]
                })

        with open(output_fn, 'w', encoding="utf-8") as f:
            json.dump(data, f, indent=4)

    def save_file(self, output_fn: str, delimiter: str, quotechar: str = '"'):
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

    def export_csv(self, output_fn: str):
        self.save_file(output_fn, ",")

    def export_tsv(self, output_fn: str):
        self.save_file(output_fn, "\t")
