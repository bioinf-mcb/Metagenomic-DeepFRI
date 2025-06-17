import json
import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory

from mDeepFRI.utils import (generate_config_json, load_fasta_as_dict,
                            retrieve_fasta_entries_as_dict, run_command)


class TestCommands(unittest.TestCase):
    def test_run_command(self):
        command = "echo 'Hello World!'"
        result = run_command(command).strip()
        self.assertEqual(result, "Hello World!")


class TestFastaFunctions(unittest.TestCase):
    def setUp(self):
        self.fasta_content = ">seq1\nATGC\n>seq2\nATGCATGC\n>seq3\nATGCATGCATGC\n"
        self.fasta_file = NamedTemporaryFile(suffix=".fasta", mode="w")
        self.fasta_file.write(self.fasta_content)
        self.fasta_file.flush()

        self.gzip_fasta_content = self.fasta_content.encode()
        self.gzip_fasta_file = NamedTemporaryFile(suffix=".fasta.gz",
                                                  mode="wb")
        self.gzip_fasta_file.write(self.gzip_fasta_content)
        self.gzip_fasta_file.flush()

    def test_load_fasta_as_dict(self):
        fasta_dict = load_fasta_as_dict(self.fasta_file.name)
        self.assertEqual(len(fasta_dict), 3)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertIn("seq3", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")
        self.assertEqual(fasta_dict["seq3"], "ATGCATGCATGC")

        # test with gzip file
        fasta_dict = load_fasta_as_dict(self.gzip_fasta_file.name)
        self.assertEqual(len(fasta_dict), 3)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertIn("seq3", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")
        self.assertEqual(fasta_dict["seq3"], "ATGCATGCATGC")

    def test_retrieve_fasta_entries_as_dict(self):
        entries = ["seq1", "seq2"]
        fasta_dict = retrieve_fasta_entries_as_dict(self.fasta_file.name,
                                                    entries)
        self.assertEqual(len(fasta_dict), 2)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")

        # test with gzip file
        fasta_dict = retrieve_fasta_entries_as_dict(self.gzip_fasta_file.name,
                                                    entries)
        self.assertEqual(len(fasta_dict), 2)
        self.assertIn("seq1", fasta_dict)
        self.assertIn("seq2", fasta_dict)
        self.assertEqual(fasta_dict["seq1"], "ATGC")
        self.assertEqual(fasta_dict["seq2"], "ATGCATGC")

    def tearDown(self):
        self.fasta_file.close()
        self.gzip_fasta_file.close()


class TestGenerateConfigJson(unittest.TestCase):
    def setUp(self):
        self.tmpdir = TemporaryDirectory()
        self.weights_dir = Path(self.tmpdir.name)

    def tearDown(self):
        self.tmpdir.cleanup()

    def create_mock_weights(self, filenames):
        for name in filenames:
            (self.weights_dir / name).touch()

    def load_config(self):
        config_path = self.weights_dir / "model_config.json"
        with open(config_path, "r") as f:
            return json.load(f)

    def test_generate_config_v1_0_success(self):
        self.create_mock_weights([
            "CNN_bp.onnx", "CNN_cc.onnx", "CNN_mf.onnx", "CNN_ec.onnx",
            "GraphConv_bp.onnx", "GraphConv_cc.onnx", "GraphConv_mf.onnx",
            "GraphConv_ec.onnx"
        ])

        generate_config_json(str(self.weights_dir), version="1.0")
        config = self.load_config()

        self.assertEqual(config["version"], "1.0")
        for net in ["cnn", "gcn"]:
            for mode in ["bp", "cc", "mf", "ec"]:
                self.assertIsNotNone(config[net][mode])
                self.assertTrue(config[net][mode].endswith(
                    f"{'CNN' if net == 'cnn' else 'GraphConv'}_{mode}.onnx"))

    def test_generate_config_v1_1_success(self):
        self.create_mock_weights([
            "CNN_bp.onnx", "CNN_cc.onnx", "CNN_mf.onnx", "GraphConv_bp.onnx",
            "GraphConv_cc.onnx", "GraphConv_mf.onnx"
        ])

        generate_config_json(str(self.weights_dir), version="1.1")
        config = self.load_config()

        self.assertEqual(config["version"], "1.1")
        for net in ["cnn", "gcn"]:
            for mode in ["bp", "cc", "mf"]:
                self.assertIn(mode, config[net])
                self.assertTrue(config[net][mode].endswith(
                    f"{'CNN' if net == 'cnn' else 'GraphConv'}_{mode}.onnx"))
            self.assertNotIn("ec", config[net])

    def test_missing_model_file_raises_error(self):
        self.create_mock_weights([
            "CNN_bp.onnx", "CNN_cc.onnx", "CNN_mf.onnx", "GraphConv_bp.onnx",
            "GraphConv_cc.onnx", "GraphConv_mf.onnx"
        ])
        with self.assertRaises(ValueError) as context:
            generate_config_json(str(self.weights_dir), version="1.0")
        self.assertIn("Model weights for cnn ec not found",
                      str(context.exception))

    def test_extra_files_are_ignored(self):
        self.create_mock_weights([
            "CNN_bp.onnx", "CNN_cc.onnx", "CNN_mf.onnx", "CNN_ec.onnx",
            "GraphConv_bp.onnx", "GraphConv_cc.onnx", "GraphConv_mf.onnx",
            "GraphConv_ec.onnx", "README.txt", "invalid_model.onnx"
        ])

        generate_config_json(str(self.weights_dir), version="1.0")
        config = self.load_config()

        for net in ["cnn", "gcn"]:
            for mode in ["bp", "cc", "mf", "ec"]:
                self.assertIsNotNone(config[net][mode])


if __name__ == "__main__":
    unittest.main()
