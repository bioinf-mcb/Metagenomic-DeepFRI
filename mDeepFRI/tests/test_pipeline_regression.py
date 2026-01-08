import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from mDeepFRI.database import Database
from mDeepFRI.pipeline import QueryFile, predict_protein_function


class TestPipelineRegression(unittest.TestCase):
    def setUp(self):
        # Locate the small_query.faa file
        self.test_dir = Path(__file__).parent
        self.data_dir = self.test_dir / "data"
        self.query_file_path = self.data_dir / "small_query.faa"

        # Verify file exists
        if not self.query_file_path.exists():
            # Fallback for when running from a different root
            self.query_file_path = Path(
                "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/vbezshapkin/Metagenomic-DeepFRI/mDeepFRI/tests/data/small_query.faa"
            )

    @patch("mDeepFRI.pipeline.Pool")
    @patch("mDeepFRI.pipeline.load_deepfri_config")
    @patch("mDeepFRI.pipeline.Predictor")
    @patch("mDeepFRI.pipeline.align_mmseqs_results")
    @patch("mDeepFRI.pipeline.extract_calpha_coords")
    @patch("mDeepFRI.pipeline.build_align_contact_map")
    @patch("mDeepFRI.pipeline.get_json_values")
    def test_predict_protein_function(
        self,
        mock_get_json_values,
        mock_build_align_contact_map,
        mock_extract_calpha_coords,
        mock_align_mmseqs_results,
        mock_predictor_cls,
        mock_load_config,
        mock_pool,
    ):
        # --- Mock setup ---

        # Mock Pool to just execute the function immediately
        mock_pool_instance = mock_pool.return_value
        mock_pool_instance.__enter__.return_value = mock_pool_instance
        mock_pool_instance.map.side_effect = lambda func, iterable: [
            func(i) for i in iterable
        ]

        # 1. Mock Config
        mock_load_config.return_value = {
            "gcn": {
                "ec": "path/to/gcn_ec.onnx",
                "bp": "path/to/gcn_bp.onnx"
            },
            "cnn": {
                "ec": "path/to/cnn_ec.onnx",
                "bp": "path/to/cnn_bp.onnx"
            },
            "version": "1.1"  # Skips 'ec' mode
        }

        # 2. Mock JSON values (GO terms and names)
        # First call gets 'goterms', second call gets 'gonames' (simplified for loop)
        mock_get_json_values.side_effect = lambda path, key: (
            ["GO:001", "GO:002"] if key == "goterms" else ["Term1", "Term2"])

        # 3. Mock MMseqs Alignment Results
        # Simulating one aligned sequence and one unaligned
        mock_aln = MagicMock()
        mock_aln.query_name = "A0A3B4WVX2_Gasdermin_pore_forming"
        mock_aln.target_name = "Target1.1"
        mock_aln.db_name = "TestDB"
        mock_aln.query_identity = 0.9
        mock_aln.query_coverage = 0.8
        mock_aln.target_coverage = 0.8
        mock_aln.query_sequence = "MFSKATANFVRQIDPEGSLIHVSRVNDSQKLVPMALVVKRNRLWFWQRPKYHPTDF"  # Truncated

        mock_align_mmseqs_results.return_value = [mock_aln]

        # 4. Mock C-alpha Coords
        mock_coords = np.zeros((len(mock_aln.query_sequence), 3))
        mock_extract_calpha_coords.return_value = [mock_coords]

        # 5. Mock Contact Map Building
        mock_cmap = np.random.rand(len(mock_aln.query_sequence),
                                   len(mock_aln.query_sequence))
        # build_align_contact_map returns (AlignmentResult, ContactMap)
        mock_build_align_contact_map.return_value = (mock_aln, mock_cmap)

        # 6. Mock Predictor
        mock_predictor_instance = mock_predictor_cls.return_value
        # Prediction output is a probability vector matching GO terms length (2)
        mock_predictor_instance.forward_pass.return_value = np.array(
            [0.95, 0.05], dtype=np.float32)

        # --- Test execution ---

        # Prepare Inputs
        query_file = QueryFile(filepath=str(self.query_file_path))
        query_file.load_sequences()

        mock_db = MagicMock(spec=Database)
        mock_db.name = "TestDB"
        mock_db.mmseqs_result = "path/to/mmseqs_results.tsv"
        mock_db.sequence_db = "path/to/seq_db"

        with tempfile.TemporaryDirectory() as temp_out:
            output_path = Path(temp_out)

            # RUN PIPELINE
            predict_protein_function(
                query_file=query_file,
                databases=(mock_db, ),
                weights="path/to/weights",
                output_path=str(output_path),
                deepfri_processing_modes=[
                    "bp"
                ],  # 'ec' filtered out by 'version' 1.1 logic mocked above
                save_cmaps=True,
                save_structures=True)

            # --- Assertions ---

            # 1. Check Output Files Created
            self.assertTrue((output_path / "alignment_summary.tsv").exists())
            self.assertTrue((output_path / "results.tsv").exists())
            self.assertTrue((output_path / "contact_maps" /
                             "A0A3B4WVX2_Gasdermin_pore_forming.npy").exists())

            # 2. Check Result Content
            with open(output_path / "results.tsv", "r") as f:
                content = f.read()
                # Expect header and at least one prediction line
                self.assertIn(
                    "protein\tnetwork_type\tprediction_mode\tgo_term\tscore",
                    content)
                self.assertIn("A0A3B4WVX2_Gasdermin_pore_forming", content)
                # GCN prediction for the aligned sequence
                self.assertIn("gcn", content)

            # 3. Verify Mock Calls
            mock_align_mmseqs_results.assert_called_once()
            mock_build_align_contact_map.assert_called()
            mock_predictor_instance.forward_pass.assert_called()


if __name__ == '__main__':
    unittest.main()
