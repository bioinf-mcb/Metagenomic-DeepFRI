import argparse
import json
import logging
from pathlib import Path

import tf2onnx
from tensorflow.keras.saving import load_model

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


# create argparse
def parse_args():

    parser = argparse.ArgumentParser(
        description="Convert tensorflow models to ONNX",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-w",
                        "--weights",
                        required=True,
                        default=None,
                        help="Path to folder containing TF models.")

    parser.add_argument("-o",
                        "--output",
                        required=False,
                        default=None,
                        help="Path to folder where ONNX model will be saved.")

    return parser.parse_args()


def convert_model(model_path, output_path):
    # IMPORTANT - IMPORT ALL CUSTOM LAYERS
    model = load_model(model_path + ".hdf5", custom_objects={})

    tf2onnx.convert.from_keras(model, output_path=str(output_path), opset=15)


# Newest CPU models

if __name__ == "__main__":

    args = parse_args()

    model_path = args.weights
    try:
        with open(model_path + "/model_config.json") as f:
            config = json.loads(f.read())
    except FileNotFoundError:
        raise FileNotFoundError(
            "File %s/model_config.json not found. Download models from"
            "https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz",
            model_path)

    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)

    for net_type in ["cnn", "gcn"]:
        for mode in ["ec", "mf", "bp", "cc"]:
            model_prefix = str(
                Path(model_path) /
                config[net_type]["models"][mode].split("/")[-1])
            model_name = Path(model_prefix).name
            out_path = Path(args.output) / f"{model_name}.onnx"
            logging.info("Converting model %s to %s", model_prefix + ".hdf5",
                         out_path)
            convert_model(model_prefix, out_path)
            logging.info(
                "Tensorflow to ONNX conversion of %s - %s is finished",
                net_type, mode)
