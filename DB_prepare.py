

from bitarray import bitarray
import pathlib
import numpy as np
import re
import matplotlib.pyplot as plt
import time

db_source = pathlib.Path('./SWISS-MODEL_Repository')
pdb_files = list(db_source.glob('**/*.pdb'))


def dst(a_atom, b_atom):
    a_position = np.array(a_atom[6:9], np.float)
    b_position = np.array(b_atom[6:9], np.float)
    return np.linalg.norm(a_position - b_position)

from CPP_lib.ContactMapper import contact_mapper
cm = contact_mapper()

# with pdb_files[2] as file:
for file in pdb_files[0:5]:
    atoms = []
    sequence = None

    with open(file, 'r') as f:
        line = f.readline()
        while line != "":
            if line.startswith("SEQRES"):
                print("SEQRES present! Use seqres", file)
                #  sequence = SEQRES
            if line.startswith("ATOM"):
                line = re.sub(' +', ' ', line)
                line = line.split(" ")
                if line[2] != "H":
                    atoms.append(line)
            line = f.readline()

    if sequence is None:
        sequence = []
        sequence.append(atoms[0][3])
        group_id = atoms[0][5]
        for i in range(1, len(atoms)):
            if atoms[i][5] == group_id:
                continue
            group_id = atoms[i][5]
            sequence.append(atoms[i][3])

    print("Sequence length", len(sequence))
    start = time.time()

    amino_i = 0
    amino_j = 0
    amino_i_group = atoms[0][5]
    amino_j_group = atoms[0][5]
    amino_contact_map = np.zeros((len(sequence), len(sequence)), dtype=np.bool)
    for i in range(0, len(atoms)):
        if atoms[i][5] != amino_i_group:
            amino_i_group = atoms[i][5]
            amino_i += 1
        amino_j = amino_i
        amino_j_group = amino_i_group
        for j in range(i+1, len(atoms)):
            if atoms[j][5] != amino_j_group:
                amino_j_group = atoms[j][5]
                amino_j += 1
            amino_contact_map[amino_i][amino_j] = \
                amino_contact_map[amino_i][amino_j] \
                or (dst(atoms[i], atoms[j]) < 6)

    print("Python generate contact map", ' %.4f' % (time.time() - start))

    start = time.time()

    save_path = pathlib.Path("./swiss_processed")/(file.stem + ".bin")
    distance = np.array([x[6:9] for x in atoms], dtype=np.float32)

    groups_count = []
    count = 0
    group_id = atoms[0][5]
    for x in atoms:
        if x[5] != group_id:
            groups_count.append(count)
            group_id = x[5]
        count += 1
    groups_count.append(count)
    groups = np.array(groups_count,dtype=np.int32)

    cm.fun(distance, groups, str(save_path.absolute()))

    print("C++ calculate and save cmap",' %.4f' % (time.time() - start))

    # plt.imshow(amino_contact_map, vmin=0, vmax=1)
    # plt.title("Generated contact map")
    # plt.show()
    #
    # n = len(sequence)
    # bits_size = int(n*(n-1) / 2)
    # bits = bitarray(bits_size)
    #
    # for i in range(n):
    #     for j in range(i+1, n):
    #         bit_index = int((n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1)
    #         bits[bit_index] = amino_contact_map[i][j]
    #
    # with open(save_path, 'wb') as f:
    #     bits.tofile(f, endian='little')

    loaded_bits = bitarray(endian="little")
    with open(save_path, 'rb') as f:
        loaded_bits.fromfile(f)

    pyamino_contact_map = amino_contact_map
    n = int((1+np.sqrt(8*len(loaded_bits)+1))/2)
    amino_contact_map = np.zeros((n, n), dtype=np.bool)
    for k in range(int(n*(n-1) / 2)):
        i = int(n - 2 - int(np.sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5))
        j = int(k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2)
        amino_contact_map[i][j] = loaded_bits[k]

    print("\n")

    plt.imshow(amino_contact_map, vmin=0, vmax=1)
    plt.title(f'C++ cmap Sequence length {len(sequence)}')
    plt.show()

    np.fill_diagonal(pyamino_contact_map, 0)
    plt.imshow(pyamino_contact_map, vmin=0, vmax=1)
    plt.title(f'Py cmap Sequence length {len(sequence)}')
    plt.show()
    #
    #
    # amino_contact_map += amino_contact_map.T
    # np.fill_diagonal(amino_contact_map, 1)
    # plt.imshow(amino_contact_map, vmin=0, vmax=1)
    # plt.title('Loaded and processed contact map')
    # plt.show()
