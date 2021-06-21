import argparse
import time

from Bio import SeqUtils
import gzip
from itertools import repeat
import multiprocessing
import numpy as np
import pathlib
import traceback

from libContactMapper import contact_mapper
from mmcif_parser import parse_mmcif
from pdb_parser import parse_pdb

import shutil


INPUT_FILE_PATTERNS = [
    '**/*.pdb',
    '**/*.pdb.gz',
    '**/*.cif',
    '**/*.cif.gz'
]


PROTEIN_LETTERS = dict()
for k, v in SeqUtils.IUPACData.protein_letters_3to1_extended.items():
    PROTEIN_LETTERS[str.upper(k)] = v
PROTEIN_LETTERS["UNK"] = "X"

MAX_CHAIN_LENGTH = 2500


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--override", action="store_true", help="(not implemented) Override existing files in output directory")
    return parser.parse_args()


def create_chunks(lst, n):
    return [lst[i::n] for i in range(n)]


def extract_camp_and_seq(protein_structure_files, save_path):
    cm = contact_mapper()
    for file in protein_structure_files:
        try:
            save_name = file.name
            if save_name.endswith('.pdb'):
                save_name = save_name.replace('.pdb', '')
                with open(file, 'r') as f:
                    atom_amino_group, positions, groups = parse_pdb(f)
            elif save_name.endswith('.pdb.gz'):
                save_name = save_name.replace('.pdb.gz', '')
                with gzip.open(file, 'rt') as f:
                    atom_amino_group, positions, groups = parse_mmcif(f)
            elif save_name.endswith('.cif'):
                save_name = save_name.replace('.cif', '')
                with open(file, 'r') as f:
                    atom_amino_group, positions, groups = parse_mmcif(f)
            elif save_name.endswith('.cif.gz'):
                save_name = save_name.replace('.cif.gz', '')
                with gzip.open(file, 'rt') as f:
                    atom_amino_group, positions, groups = parse_mmcif(f)
            else:
                print("Unsupported file format of file " + str(file))
                continue
        except Exception:
            print("EXCEPTION DURING FILE READING ", str(file))
            trace = traceback.format_exc()
            print(trace)
            continue

        try:
            _, groups = np.unique(groups, return_index=True)
            groups.sort()

            if len(groups) < 9:
                # print("Proteins shorter than 9 amino acids are not supported", file)
                continue

            if len(groups) > MAX_CHAIN_LENGTH:
                group_indexes = groups[:MAX_CHAIN_LENGTH]
                group_indexes = np.append(group_indexes, groups[MAX_CHAIN_LENGTH]).astype(np.int32)
            else:
                group_indexes = np.append(groups, positions.shape[0]).astype(np.int32)

            sequence = ''.join([PROTEIN_LETTERS[atom_amino_group[i]] for i in group_indexes[:-1]])
            with open(save_path + "/seq/" + save_name + ".faa", "w") as f:
                f.write(">" + save_name + "\n" + sequence + "\n")

            cm.save_contact_map(positions, group_indexes, save_path + "/cmap/" + save_name + ".bin")

        except Exception:
            print("EXCEPTION DURING FILE PROCESSING ", str(file))
            trace = traceback.format_exc()
            print(trace)
            continue


def main():
    args = parse_args()
    input_dir = pathlib.Path(args.input)
    output_dir = pathlib.Path(args.output)

    output_dir.mkdir(exist_ok=True)
    (output_dir / 'seq').mkdir(exist_ok=True)
    (output_dir / 'cmap').mkdir(exist_ok=True)
    save_path = str(output_dir.absolute())
    print("Save path: ", save_path)

    for pattern in INPUT_FILE_PATTERNS:
        protein_structure_files = np.array(list(input_dir.glob(pattern)))
        if len(protein_structure_files) == 0:
            continue

        file_names = [x.name for x in protein_structure_files]
        _, index = np.unique(file_names, return_index=True)
        protein_structure_files = protein_structure_files[index]
        print("Processing ", len(protein_structure_files), pattern, " files.")

        x = int(len(protein_structure_files) / 1000)
        chunks = create_chunks(protein_structure_files, max(multiprocessing.cpu_count(), x))

        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as p:
            p.starmap(extract_camp_and_seq, zip(chunks,repeat(save_path)))

    sequence_files = list((output_dir / 'seq').glob("**/*.faa"))
    print("Merging " + str(len(sequence_files)) + " faa files.")
    with open(output_dir / 'merged_sequences.faa', 'wb') as writer:
        for seq_file in sequence_files:
            with open(seq_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)


if __name__ == '__main__':
    main()
