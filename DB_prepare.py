import pathlib
import numpy as np
from libContactMapper import contact_mapper
from Bio import SeqUtils
import gzip
from pdb_parser import parse_pdb
from mmcif_parser import parse_mmcif

output_dir = pathlib.Path("../genomics_data/database")
output_dir.mkdir(exist_ok=True)
(output_dir/'seq').mkdir(exist_ok=True)
(output_dir/'cmap').mkdir(exist_ok=True)
save_path = str(output_dir.absolute())

db_source = pathlib.Path('./SWISS-MODEL_Repository')
protein_structure_files = list(db_source.glob('**/*.pdb'))

cm = contact_mapper()

protein_letters = dict()
for k,v in SeqUtils.IUPACData.protein_letters_3to1_extended.items():
    protein_letters[str.upper(k)] = v

# file = pathlib.Path('/home/soliareofastora/Desktop/1a0a.cif')
with pathlib.Path('/home/soliareofastora/Desktop/1a0a.cif') as file:
# for file in protein_structure_files:

    save_name = file.name
    if save_name.endswith('.pdb'):
        save_name = save_name.replace('.pdb','')
        with open(file, 'r') as f:
            sequence, positions, groups = parse_pdb(f)
    elif save_name.endswith('.cif'):
        save_name = save_name.replace('.cif', '')
        with open(file, 'r') as f:
            sequence, positions, groups = parse_mmcif(f)
    elif save_name.endswith('.cif.gz'):
        save_name = save_name.replace('.cif.gz', '')
        with gzip.open(file,'rt') as f:
            sequence,positions,groups = parse_mmcif(f)
    else:
        print("Unsupported file format")
        continue

    _, groups_index = np.unique(groups, return_index=True)
    groups_index.sort()

    if len(groups_index) < 9:
        print("Proteins shorter than 9 amino groups are not supported")
        continue

    sequence = ''.join([protein_letters[sequence[i]] for i in groups_index])
    with open(save_path + "/seq/" + save_name + ".faa", "w") as f:
        f.write(">" + save_name + "\n" + sequence + "\n")

    groups_index = np.append(groups_index, positions.shape[0]).astype(np.int32)
    cm.save_contact_map(positions, groups_index, save_path + "/cmap/" + save_name + ".bin")
