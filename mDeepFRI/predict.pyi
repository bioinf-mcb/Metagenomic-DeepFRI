"""Type stubs for predict module."""

from typing import Optional

import numpy as np
from numpy.typing import NDArray

def seq2onehot(seq: str) -> NDArray[np.float32]:
    """
    Converts a protein sequence to 26-dim one-hot encoding.

    Args:
        seq: Protein sequence (single-letter amino acid codes).

    Returns:
        One-hot encoding of the protein sequence with shape (L, 26)
        where L is the sequence length.

    Raises:
        ValueError: If sequence contains invalid characters.
    """
    ...

class Predictor:
    """
    ONNX-based predictor for DeepFRI protein function annotation.

    This class loads trained DeepFRI models in ONNX format and performs
    forward passes for protein function prediction. It supports both GCN-based
    predictions (with contact maps) and CNN-based predictions (sequence only).

    The predictor automatically uses CUDA acceleration when available,
    falling back to CPU execution otherwise.

    Attributes:
        model_path: Path to the ONNX model file.
        threads: Number of threads for CPU inference.
            0 means use default (recommended).
        session: ONNX Runtime inference session.
        input_names: Names of model input tensors.

    Example:
        >>> # Sequence-only prediction
        >>> predictor = Predictor("DeepCNN-MERGED_mf.onnx", threads=4)
        >>> scores = predictor.forward_pass("MSKGEELFT")
        >>>
        >>> # Structure-based prediction with contact map
        >>> predictor = Predictor("DeepFRI-MERGED_mf.onnx")
        >>> scores = predictor.forward_pass("MSKGEELFT", cmap=contact_map)

    Note:
        For best performance, set threads=0 to let ONNX Runtime decide
        optimal thread count. Only adjust if you encounter threading issues.
    """

    model_path: str
    threads: int
    session: object
    input_names: list[str]

    def __init__(self, model_path: str, threads: int = 1) -> None:
        """
        Initialize predictor with model file.

        Args:
            model_path: Path to ONNX model file.
            threads: Number of threads for CPU inference.
                Defaults to 1.

        Raises:
            FileNotFoundError: If model_path does not exist.
            RuntimeError: If model loading fails.
        """
        ...

    def _load_model(self) -> None:
        """Load ONNX model into inference session."""
        ...

    def forward_pass(self, seqres: str, cmap: Optional[NDArray[np.float32]] = None) -> NDArray[np.float32]:
        """
        Perform forward pass through DeepFRI model.

        Runs inference on a protein sequence with optional contact map to
        predict GO terms or EC numbers. The model uses GCN when a contact map
        is provided, otherwise uses CNN for sequence-only prediction.

        Args:
            seqres: Protein sequence (single-letter amino acid codes).
                Should only contain standard amino acids and gap character '-'.
            cmap: Protein contact map (L x L matrix) where L is sequence length.
                If None, sequence-only prediction is performed. Defaults to None.

        Returns:
            1D array of prediction scores for each GO term/EC number
            in the model's vocabulary. Shape: (num_classes,)

        Raises:
            ValueError: If sequence contains invalid characters.
            RuntimeError: If inference fails.

        Example:
            >>> predictor = Predictor("model.onnx")
            >>> # Sequence-only
            >>> scores = predictor.forward_pass("MSKGEELFT")
            >>> # With structure
            >>> scores = predictor.forward_pass("MSKGEELFT", cmap=contact_map)
            >>> top_idx = np.argmax(scores)
            >>> print(f"Top prediction score: {scores[top_idx]:.3f}")

        Note:
            - Contact maps should be binary or distance-based matrices
            - Prediction scores are not probabilities (no softmax applied)
            - Higher scores indicate higher confidence
        """
        ...
