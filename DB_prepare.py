import numpy
from bitarray import bitarray
import pathlib
import numpy as np
import re
import matplotlib.pyplot as plt
import time
import tqdm
from libContactMapper import contact_mapper

db_source = pathlib.Path('./SWISS-MODEL_Repository')
pdb_files = list(db_source.glob('**/*.pdb'))

cm = contact_mapper()

# file = pdb_files[0]
# for file in tqdm.tqdm(pdb_files[0:10]):
for file in pdb_files[0:1000]:

    sequence = []
    positions = []
    groups = []
    save_path = pathlib.Path("./swiss_processed")/(file.stem + ".bin")

    with open(file, 'r') as f:
        line = f.readline()
        while line != "":
            if line.startswith("SEQRES"):
                print("SEQRES present! Use seqres", file)
            elif line.startswith("ATOM"):
                if line[77] != 'H':
                    sequence.append(line[17:20])
                    positions.append([line[30:38], line[38:46], line[46:54]])
                    groups.append(line[22:26])
            elif line.startswith("TER"):
                break
            line = f.readline()

    # start = time.time()

    groups = np.array(groups, dtype=np.int32)
    unique, index = np.unique(groups, return_index=True)
    sequence = np.array(sequence)[index]
    if len(sequence) < 9:
        print("Proteins shorter than 8 amino groups are not supported")

    positions = np.array(positions, dtype=np.float32)
    index = np.append(index, groups.shape[0]).astype(np.int32)
    cm.generate_contact_map(positions, index, str(save_path.absolute()))

    # print("C++ calculate and save cmap",' %.4f' % (time.time() - start))
    #
    #
    # loaded_bits = bitarray(endian="little")
    # with open(save_path, 'rb') as f:
    #     loaded_bits.fromfile(f)
    #
    # n = int((1+np.sqrt(8*len(loaded_bits)+1))/2)
    # amino_contact_map = np.zeros((n, n), dtype=np.bool)
    # for k in range(int(n*(n-1) / 2)):
    #     i = int(n - 2 - int(np.sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5))
    #     j = int(k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2)
    #     amino_contact_map[i][j] = loaded_bits[k]
    #
    # plt.imshow(amino_contact_map, vmin=0, vmax=1)
    # plt.title(f'C++ cmap Sequence length {len(sequence)}')
    # plt.show()


    # np.fill_diagonal(pyamino_contact_map, 0)
    # plt.imshow(pyamino_contact_map, vmin=0, vmax=1)
    # plt.title(f'Py cmap Sequence length {len(sequence)}')
    # plt.show()
    #
    #
    # amino_contact_map += amino_contact_map.T
    # np.fill_diagonal(amino_contact_map, 1)
    # plt.imshow(amino_contact_map, vmin=0, vmax=1)
    # plt.title('Loaded and processed contact map')
    # plt.show()


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
