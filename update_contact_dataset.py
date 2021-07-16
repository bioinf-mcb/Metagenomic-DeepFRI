import argparse
import gzip
import multiprocessing
import shutil
import traceback

from itertools import repeat
import numpy as np

from CONFIG import *
from CPP_lib.libAtomDistanceIO import save_atoms
from structure_files_parsers.parse_mmcif import parse_mmcif
from structure_files_parsers.parse_pdb import parse_pdb
from utils import create_chunks


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=False, default=STRUCTURE_FILES_PATH)
    parser.add_argument("-o", "--output", required=False, default=CONTACT_MAP_DATASET_PATH)
    parser.add_argument("--overwrite", action="store_true",
                        help="Override existing")
    return parser.parse_args()


def process_files(protein_structure_files, save_path):
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

            save_atoms(positions, group_indexes, save_path + "/cmap/" + save_name + ".bin")

        except Exception:
            print("EXCEPTION DURING FILE PROCESSING ", str(file))
            trace = traceback.format_exc()
            print(trace)
            continue


def update_contact_dataset(input_dir, output_dir, overwrite=False):
    output_dir.mkdir(exist_ok=True, parents=True)
    (output_dir / 'seq').mkdir(exist_ok=True)
    (output_dir / 'cmap').mkdir(exist_ok=True)
    save_path = str(output_dir.absolute())
    print("output path: ", save_path)

    structure_files = dict()

    for pattern in STRUCTURE_FILES_PATTERNS:
        pattern_structures = np.array(list(input_dir.glob("**/*" + pattern)))
        pattern_ids = np.array([x.name[:-len(pattern)] for x in pattern_structures])
        if len(pattern_structures) == 0:
            continue

        _, index = np.unique(pattern_ids, return_index=True)
        pattern_structures = pattern_structures[index]
        pattern_ids = pattern_ids[index]
        print("Found", len(pattern_structures), pattern, "files")

        structure_files.update(zip(pattern_ids, pattern_structures))

    if len(structure_files) == 0:
        print("No structure files found")
        return

    if not overwrite:
        existing_file_ids = set([x.name[:-4] for x in (output_dir / 'cmap').glob("**/*.bin")])
        print("Existing files:", len(existing_file_ids))
        for id in list(structure_files.keys()):
            if id in existing_file_ids:
                structure_files.pop(id)

    x = int(len(structure_files) / 10000)
    chunks = create_chunks(list(structure_files.values()), max(multiprocessing.cpu_count(), x))

    print("Processing", len(structure_files), "files")
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as p:
        p.starmap(process_files, zip(chunks, repeat(save_path)))

    sequence_files = list((output_dir / 'seq').glob("**/*.faa"))
    print("Merging " + str(len(sequence_files)) + " sequence files.")
    with open(output_dir / 'merged_sequences.faa', 'wb') as writer:
        for seq_file in sequence_files:
            with open(seq_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)


if __name__ == '__main__':
    args = parse_args()
    input_dir = pathlib.Path(args.input)
    output_dir = pathlib.Path(args.output)
    overwrite = args.overwrite
    update_contact_dataset(input_dir, output_dir, overwrite)
