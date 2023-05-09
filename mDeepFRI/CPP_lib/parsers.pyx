from libc.stdlib cimport strtof
from libcpp.string cimport string
from libcpp.vector cimport vector

from io import TextIOWrapper
from typing import Tuple

import numpy as np

cimport numpy as cnp

cnp.import_array()

cdef tuple parse_atom_line(char* line):
    cdef char* seq = &line[17]
    cdef char* group = &line[21]
    # assign
    cdef float x = strtof(&line[30], NULL)
    cdef float y = strtof(&line[38], NULL)
    cdef float z = strtof(&line[46], NULL)

    return (seq[0:3], group[0:5], x, y, z)


def parse_pdb(
        lines):
    cdef vector[string] sequence = []
    cdef vector[string] groups = []
    cdef vector[float] positions = []
    cdef char* line
    cdef float x, y, z
    cdef string seq, group

    for line in lines:
        if line == b"":
            break
        if line.startswith(b"TER"):
            break
        if line.startswith(b"ATOM"):
            if len(line) == 80:
                if line[76] != b'H' and line[17] != b' ':
                    seq, group, x, y, z = parse_atom_line(line)

                    positions.push_back(x)
                    positions.push_back(y)
                    positions.push_back(z)

                    sequence.push_back(seq)
                    groups.push_back(group)

    return sequence, positions, groups


def parse_mmcif(
        lines):
    """Parse a mmCIF file.

    Args:
        file (TextIOWrapper): Open file instance.

    Returns:
        sequence (np.array): Aminoacid sequence.
        positions (np.array): Positions of the atoms.
        groups (np.array): Groups of the atoms.
    """

    cdef vector[string] sequence = []
    cdef vector[string] groups = []
    cdef vector[float] positions = []
    cdef int n_atoms = -1
    cdef int label_counter = -1
    cdef bytes line
    cdef int index

    for line in lines:
        if line.startswith(b'_refine_hist.pdbx_number_atoms_protein'):
            n = line.split(b" ")[-1]
            if n != b'?' and n != b'0':
                try:
                    n_atoms = int(n)
                except ValueError:
                    continue

        if line.startswith(b'_atom_site.'):
            label_counter += 1
            if line.startswith(b'_atom_site.type_symbol '):
                atom_symbol = label_counter
            if line.startswith(b'_atom_site.label_asym_id '):
                assembly = label_counter
            if line.startswith(b'_atom_site.label_seq_id '):
                sequence_id = label_counter
            if line.startswith(b'_atom_site.label_comp_id '):
                residue = label_counter
            if line.startswith(b'_atom_site.Cartn_x '):
                x = label_counter
            if line.startswith(b'_atom_site.Cartn_y '):
                y = label_counter
            if line.startswith(b'_atom_site.Cartn_z '):
                z = label_counter

        if line.startswith(b"ATOM"):
            if label_counter > 0:
                break

    if n_atoms > 0:
        index = 0
        for line in lines:
            if line.startswith(b'loop_'):
                break
            if line.startswith(b'ATOM'):
                atom = line.split(b" ")
                if len(atom[residue]) == 3 and atom[atom_symbol] != b"H":

                    group = b''.join([atom[assembly], atom[sequence_id]])
                    sequence.push_back(atom[residue])
                    groups.push_back(group)
                    positions.push_back(float(atom[x]))
                    positions.push_back(float(atom[y]))
                    positions.push_back(float(atom[z]))

                    index += 1
                    if index == n_atoms:
                        break

    else:
        for line in lines:
            if line.startswith(b'loop_'):
                continue

            if line.startswith(b'ATOM'):
                atom = line.split()
                if len(atom[residue]) == 3 and atom[atom_symbol] != b"H":
                    sequence.push_back(atom[residue])
                    groups.push_back(b''.join([atom[assembly], atom[sequence_id]]))
                    positions.push_back(float(atom[x]))
                    positions.push_back(float(atom[y]))
                    positions.push_back(float(atom[z]))

    return sequence, positions, groups
