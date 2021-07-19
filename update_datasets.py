import argparse
import gzip
import multiprocessing
import shutil
import time
import traceback

from itertools import repeat
import numpy as np

from CONFIG import *
from CPP_lib.libAtomDistanceIO import save_atoms
from CPP_lib.libAtomDistanceIO import initialize as initialize_CPP_LIB
from structure_files_parsers.parse_mmcif import parse_mmcif
from structure_files_parsers.parse_pdb import parse_pdb
from utils import create_chunks, run_command


# Use this script to extract atoms positions from large text-based PDB/cif files
# and save them in space efficient binary file.
# for file specification please head over to CPP_lib
# Sequences will be extracted along the way and stored in mmseqs2 db

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=False, default=STRUCTURE_FILES_PATH)
    parser.add_argument("--atoms", required=False, default=ATOMS_DATASET_PATH)
    parser.add_argument("-o", "--output", required=False, default=MMSEQS_DATABASES_PATH)
    parser.add_argument("-t", "--temporary", required=False, default=TMP_FOLDER_PATH)
    parser.add_argument("--overwrite", action="store_true", help="Override existing")
    return parser.parse_args()


def load_file_extract_and_save_atoms(protein_structure_files, save_path):
    for file in protein_structure_files:
        try:
            file_id = file.name
            if file_id.endswith('.pdb'):
                file_id = file_id.replace('.pdb', '')
                with open(file, 'r') as f:
                    atom_amino_group, positions, groups = parse_pdb(f)
            elif file_id.endswith('.pdb.gz'):
                file_id = file_id.replace('.pdb.gz', '')
                with gzip.open(file, 'rt') as f:
                    atom_amino_group, positions, groups = parse_mmcif(f)
            elif file_id.endswith('.cif'):
                file_id = file_id.replace('.cif', '')
                with open(file, 'r') as f:
                    atom_amino_group, positions, groups = parse_mmcif(f)
            elif file_id.endswith('.cif.gz'):
                file_id = file_id.replace('.cif.gz', '')
                with gzip.open(file, 'rt') as f:
                    atom_amino_group, positions, groups = parse_mmcif(f)
            else:
                print("Unsupported file format of file " + str(file))
                continue
        except Exception:
            print("EXCEPTION WHILE READING FILE ", str(file))
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
            with open(save_path + "/seq/" + file_id + ".faa", "w") as f:
                f.write(">" + file_id + "\n" + sequence + "\n")
            save_atoms(positions, group_indexes, save_path + "/positions/" + file_id + ".bin")
        except Exception:
            print("EXCEPTION DURING FILE PROCESSING ", str(file))
            trace = traceback.format_exc()
            print(trace)
            continue


def update_atoms_dataset(input_path, atoms_path, output_path, tmp_path, overwrite):
    atoms_path.mkdir(exist_ok=True, parents=True)
    (atoms_path / 'seq').mkdir(exist_ok=True)
    (atoms_path / 'positions').mkdir(exist_ok=True)
    save_path = str(atoms_path.absolute())
    print("Atoms will be stored in: ", save_path)

    structure_files = dict()

    for pattern in STRUCTURE_FILES_PATTERNS:
        pattern_structures = np.array(list(input_path.glob("**/*" + pattern)))
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
        existing_file_ids = set([x.name[:-4] for x in (atoms_path / 'positions').glob("**/*.bin")])
        print("Existing files:", len(existing_file_ids))
        for id in list(structure_files.keys()):
            if id in existing_file_ids:
                structure_files.pop(id)

    x = int(len(structure_files) / 10000)
    chunks = create_chunks(list(structure_files.values()), max(multiprocessing.cpu_count(), x))

    initialize_CPP_LIB()
    print("Processing", len(structure_files), "files")
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as p:
        p.starmap(load_file_extract_and_save_atoms, zip(chunks, repeat(save_path)))

    sequence_files = list((atoms_path / 'seq').glob("**/*.faa"))
    print("Merging " + str(len(sequence_files)) + " sequence files.")
    with open(atoms_path / 'merged_sequences.faa', 'wb') as writer:
        for seq_file in sequence_files:
            with open(seq_file, 'rb') as reader:
                shutil.copyfileobj(reader, writer)

    # create new mmseqs2 database
    output_path.mkdir(exist_ok=True, parents=True)
    creation_time = str(time.time())
    db_path = (output_path / creation_time)
    while db_path.exists():
        creation_time = str(int(time.time()))
        db_path = (output_path / creation_time)
    db_path.mkdir()
    print("Creating new mmseqs2 database in " + str(db_path))

    # todo save contact_dir if more than one contact datasets are required
    run_command(f"mmseqs createdb {atoms_path / 'merged_sequences.faa'} {db_path / DEFAULT_MMSEQS_NAME} --dbtype 1")
    run_command(f"mmseqs createindex {db_path / DEFAULT_MMSEQS_NAME} {tmp_path}")


if __name__ == '__main__':
    args = parse_args()

    input_path = pathlib.Path(args.input)
    atoms_path = pathlib.Path(args.atoms)
    output_path = pathlib.Path(args.output)
    tmp_path = pathlib.Path(args.temporary)
    overwrite = args.overwrite

    update_atoms_dataset(input_path, atoms_path, output_path, tmp_path, overwrite)

