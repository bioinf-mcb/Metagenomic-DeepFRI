#ifndef ATOMS_FILE_IO
#define ATOMS_FILE_IO

#include <fstream>

void SaveAtomsFile(float *position_array, size_t atom_count, int *groups_array,
                   size_t chain_length, const std::string &save_path);

// static std::tuple<int, int *, float *>
// LoadAtomsFile(const std::string &file_path);

#endif
