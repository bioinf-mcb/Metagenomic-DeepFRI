import gzip
import time

import numpy as np

from CONFIG import *
from CPP_lib.libAtomDistanceIO import save_atoms, load_contact_map, initialize
from structure_files_parsers.parse_mmcif import parse_mmcif
from structure_files_parsers.parse_pdb import parse_pdb



file = pathlib.Path('/home/soliareofastora/data/structure_files/Q8/0T/I0/swissmodel/87_199_5yqr.1.A_60823e8737a6a2a0adacc9f8.pdb')

save_path = '/home/soliareofastora/tmp'

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
# with open(save_path + "/seq/" + save_name + ".faa", "w") as f:
#     f.write(">" + save_name + "\n" + sequence + "\n")
print(len(sequence))

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