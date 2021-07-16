import gzip
import time

import numpy as np

from CONFIG import *
from CPP_lib.libAtomDistanceIO import save_atoms, load_contact_map, initialize
from structure_files_parsers.parse_mmcif import parse_mmcif
from structure_files_parsers.parse_pdb import parse_pdb



file = pathlib.Path('/home/soliareofastora/genomics_data/PDB_Compressed/5ezz.cif.gz')

save_path = '/home/soliareofastora/genomics_data/tmp'

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



_, groupsx = np.unique(groups, return_index=True)
groupsx.sort()


group_indexes = np.append(groupsx, positions.shape[0]).astype(np.int32)
len(group_indexes)
sequence = ''.join([PROTEIN_LETTERS[atom_amino_group[i]] for i in group_indexes[:-1]])
with open(save_path + "/seq/" + save_name + ".faa", "w") as f:
    f.write(">" + save_name + "\n" + sequence + "\n")


start = time.time()
initialize()

print(time.time()-start)
start = time.time()
save_atoms(positions, group_indexes, save_path + "/cmap/" + save_name + ".bin")

print(time.time()-start)
start = time.time()
loaded = load_contact_map(save_path + "/cmap/" + save_name + ".bin", 6)

print(time.time()-start)
start = time.time()
import matplotlib.pyplot as plt

plt.imshow(loaded)
plt.show()