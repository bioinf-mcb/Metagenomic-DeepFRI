import numpy as np


def parse_pdb(file):
    sequence = []
    positions = []
    groups = []

    line = file.readline()
    while line != "":
        if line.startswith("ATOM"):
            if line[76] != 'H' and line[17] != ' ':
                sequence.append(line[17:20])
                positions.append([line[30:38], line[38:46], line[46:54]])
                groups.append(line[21:26])
        line = file.readline()

    positions = np.array(positions, dtype=np.float32)
    return sequence, positions, groups
