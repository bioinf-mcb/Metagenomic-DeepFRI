import argparse
import time

from update_contact_dataset import update_contact_dataset
from utils import run_command

from CONFIG import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=False, default=STRUCTURE_FILES_PATH)
    parser.add_argument("-c", "--contacts", required=False, default=CONTACT_MAP_DATASET_PATH)
    parser.add_argument("-o", "--output", required=False, default=TARGET_MMSEQS_DATABASE_PATH)
    parser.add_argument("--dont_update_contacts",action="store_true", required=False, help="dont update contact map dataset")
    parser.add_argument("--overwrite", action="store_true", required=False, help="Override existing contacts")
    return parser.parse_args()


def create_mmseq_db(contact_dir, output_dir):
    output_dir.mkdir(exist_ok=True, parents=True)
    creation_time = str(time.time())
    db_path = (output_dir / creation_time)

    while db_path.exists():
        creation_time = str(int(time.time()))
        db_path = (output_dir / creation_time)
    db_path.mkdir()

    # todo save contact_dir if more than one contact datasets are required
    run_command(f"mmseqs createdb {contact_dir / 'merged_sequences.faa'} {db_path / DEFAULT_MMSEQS_NAME} --dbtype 1")
    run_command(f"mmseqs createindex {db_path / DEFAULT_MMSEQS_NAME} {TMP_FOLDER_PATH}")


if __name__ == '__main__':
    args = parse_args()
    input_dir = pathlib.Path(args.input)
    contact_dir = pathlib.Path(args.contacts)
    output_dir = pathlib.Path(args.output)
    dont_update_contacts = args.dont_update_contacts
    overwrite = args.overwrite

    if not dont_update_contacts:
        print("Updating contact map dataset with files in", input_dir, ". Use --dont_update_contacts to skip this step")
        update_contact_dataset(input_dir, contact_dir, overwrite)

    create_mmseq_db(contact_dir, output_dir)
