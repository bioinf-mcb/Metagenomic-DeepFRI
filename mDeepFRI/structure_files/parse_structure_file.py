import dataclasses

import numpy as np

from mDeepFRI.parsers import parse_mmcif, parse_pdb

# need to parse different type of files? Add a file name pattern with a parser in this dict.
# read structure_files_parsers/README.md for more information about how to create new parser.
PARSERS = {
    '.pdb': parse_pdb,
    '.pdb.gz': parse_pdb,
    '.cif': parse_mmcif,
    '.cif.gz': parse_mmcif,
    '.ent': parse_pdb,
    '.ent.gz': parse_pdb
}


@dataclasses.dataclass
class SeqAtoms:
    protein_id: str
    atom_amino_group: np.ndarray
    positions: np.ndarray
    groups: np.ndarray


def search_structure_files(input_paths: list):
    """

    :param input_paths:
    :return:
    """
    structure_files_paths = dict()
    # Iterate input_paths and check if file extension matches any key in structure_files.PARSERS
    for input_path in input_paths:
        print(f"{str(input_path)}")
        if input_path.is_file():
            for pattern in PARSERS:
                if str(input_path).endswith(pattern):
                    print(f"\tFile {input_path} is {pattern}")
                    structure_file_id = input_path.name[:-len(pattern)]
                    structure_files_paths[structure_file_id] = input_path
                    continue
        elif input_path.is_dir():
            files_inside_input_path = list(input_path.glob("**/*"))
            for pattern in PARSERS:
                structure_file_paths = list(
                    filter(lambda x: str(x).endswith(pattern),
                           files_inside_input_path))
                if len(structure_file_paths) == 0:
                    continue
                print(
                    f"\tFound {len(structure_file_paths)} {pattern} files in {input_path}"
                )
                structure_file_ids = [
                    x.name[:-len(pattern)] for x in structure_file_paths
                ]
                structure_files_paths.update(
                    zip(structure_file_ids, structure_file_paths))
        else:
            print(f"\tUnable to find {str(input_path)}")
    return structure_files_paths
